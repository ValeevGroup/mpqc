//
// csgrad34qb.cc
// based on: csgrad.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Ida Nielsen <ibniels@ca.sandia.gov>
// Maintainer: LPS
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#ifdef __GNUC__
#pragma implementation
#endif

#include <math.h>

#include <util/misc/formio.h>
#include <chemistry/qc/mbpt/bzerofast.h>
#include <chemistry/qc/mbpt/csgrad34qb.h>
#include <chemistry/qc/mbpt/util.h>
#include <chemistry/qc/basis/distshpair.h>

using namespace std;
using namespace sc;

CSGrad34Qbtr::CSGrad34Qbtr(int mythread_a, int nthread_a,
                           int me_a, int nproc_a,
                           const Ref<MemoryGrp> &mem_a,
                           const Ref<MessageGrp> &msg_a,
                           const Ref<ThreadLock> &lock_a,
                           const Ref<GaussianBasisSet> &basis_a,
                           const Ref<TwoBodyInt> &tbint_a,
                           const Ref<TwoBodyDerivInt> &tbintder_a,
                           int nocc_a, int nfzc_a,
                           double **scf_vector_a,
                           double tol_a, int debug_a,
                           int dynamic_a, double print_percent_a,
                           DistShellPair::SharedData *shellpair_shared_data,
                           int dograd_a, int natom_a):
  shellpair_shared_data_(shellpair_shared_data)
{
  msg = msg_a;
  mythread = mythread_a;
  nthread = nthread_a;
  lock = lock_a;
  basis = basis_a;
  tbint = tbint_a;
  tbintder = tbintder_a;
  nocc = nocc_a;
  nfzc = nfzc_a;
  me = me_a;
  nproc = nproc_a;
  tol = tol_a;
  mem = mem_a;
  scf_vector = scf_vector_a;
  debug = debug_a;
  dynamic_ = dynamic_a;
  print_percent_ = print_percent_a;
  dograd = dograd_a;
  natom = natom_a;

  Lpi = 0;

  aointder_computed = 0;
  timer = new RegionTimer();

  ginter = new double*[natom];
  for (int i=0; i<natom; i++) {
    ginter[i] = new double[3];
    }
}

CSGrad34Qbtr::~CSGrad34Qbtr()
{
  delete[] Lpi;
  for (int i=0; i<natom; i++) {
    delete[] ginter[i];
    }
  delete[] ginter;
}

void
CSGrad34Qbtr::run()
{
  ////////////////////////////////////////////////////////////
  // Perform third and fourth quarter back-transformation and
  // compute contribution to gradient from non-sep 2PDM
  ////////////////////////////////////////////////////////////

  int P, Q, R, S;
  int p, q, r, s;
  int offset, p_offset, q_offset, r_offset, s_offset;
  int np, nq, nr, ns;
  int nbasis = basis->nbasis();
  int nshell = basis->nshell();
  int nfuncmax = basis->max_nfunction_in_shell();;
  int i, j;
  int jloop;
  int ij_proc, ij_index;
  int int_index;
  int index;
  int bf1, bf2, bf3, bf4;
  int nocc_act = nocc - nfzc;
  int qp, sr;
  int factor_pqrs;
  double c_rj;
  double pqrs;
  double tmpval;
  double dtol = 1.0e-10;
  double *c_sj, *c_pi, *c_qi;
  double *gammabuf;
  double *gamma_iqrs, *gamma_pqrs;
  double *gamma_iqjs_ptr, *gamma_irjq_ptr;
  double *gamma_iqrs_ptr, *gamma_iqsr_ptr, *gamma_iprs_ptr;
  double *gamma_pqrs_ptr;
  double *lpi_ptr, *lqi_ptr;
  double *grad_ptr1, *grad_ptr2;
  const double *intbuf = tbint->buffer();
  const double *intderbuf = tbintder->buffer(); // AO integral derivative buffer

  delete[] Lpi;
  Lpi = new double[nbasis*ni];
  
  // Initialize Lpi and ginter
  memset(Lpi, 0, sizeof(double)*basis->nbasis()*ni);
  for (i=0; i<natom; i++) memset(ginter[i], 0, sizeof(double)*3);
  
  MemoryGrpBuf<double> membuf_remote(mem);

  gamma_iqrs = new double[ni*nbasis*nfuncmax*nfuncmax];
  if (!gamma_iqrs) {
    ExEnv::errn() << "Could not allocate gamma_iqrs" << endl;
    abort();
    }

  gamma_pqrs = new double[nfuncmax*nfuncmax*nfuncmax*nfuncmax];
  if (!gamma_pqrs) {
    ExEnv::errn() << "Could not allocate gamma_pqrs" << endl;
    abort();
    }
  
  DerivCenters der_centers;

  DistShellPair shellpairs(msg,nthread,mythread,lock,basis,basis,dynamic_,
                           shellpair_shared_data_);
  shellpairs.set_print_percent(print_percent_);
  shellpairs.set_debug(debug);
  if (debug) shellpairs.set_print_percent(1);
  S = 0;
  R = 0;
  while (shellpairs.get_task(S,R)) {
    // If both PQRS and PQRS derivative are zero, skip this S,R pair
    // NB: The test is done after assigning an SR pair, and, when
    // using static load balancing, this may create some load imbalance
    // if more SR pairs are discarded in some threads than in others
    if (tbint->log2_shell_bound(R,S) < tol
        && (dograd && tbintder->log2_shell_bound(R,S) < tol)) continue;

    ns = basis->shell(S).nfunction();
    s_offset = basis->shell_to_function(S);
    nr = basis->shell(R).nfunction();
    r_offset = basis->shell_to_function(R);
    
    timer->enter("3. q.b.t.");
    // Begin third quarter back-transformation.

    bzerofast(gamma_iqrs,ni*nbasis*nfuncmax*nfuncmax);

    for (i=0; i<ni; i++) {
      for (jloop=me; jloop<me+nocc_act; jloop++) {
        // stagger j's to minimize contention
        j = jloop%nocc_act + nfzc;  // j runs from nfzc to nocc
        ij_proc =  (i*nocc + j)%nproc; // ij_proc has this ij pair
        ij_index = (i*nocc + j)/nproc;

        offset = s_offset*nbasis + ij_index*nbasis*nbasis;
        // Send for elements gamma_iqjs, if necessary
        gammabuf = (double*) membuf_remote.readonly_on_node(offset,
                                                            nbasis * ns,
                                                            ij_proc);
        for (bf1=0; bf1<nr; bf1++) {
          c_rj = scf_vector[bf1 + r_offset][j];
          gamma_iqjs_ptr = gammabuf;
          for (bf2=0; bf2<ns; bf2++) {
            gamma_iqrs_ptr = &gamma_iqrs[bf2 + ns*nbasis*(bf1 + nr*i)];
            for (q=0; q<nbasis; q++) {
              *gamma_iqrs_ptr += c_rj * *gamma_iqjs_ptr++;
              gamma_iqrs_ptr += ns;
              } // exit q loop
            }   // exit bf2 loop
          }     // exit bf1 loop

        membuf_remote.release();

        offset = r_offset*nbasis + ij_index*nbasis*nbasis;
        // Send for elements gamma_irjq, if necessary
        gammabuf = (double*) membuf_remote.readonly_on_node(offset,
                                                            nbasis*nr,
                                                            ij_proc);
        for (bf1=0; bf1<ns; bf1++) {
          s = bf1 + s_offset;
          c_sj = &scf_vector[s][j];
          gamma_irjq_ptr = gammabuf;
          for (bf2=0; bf2<nr; bf2++) {
            r = bf2 + r_offset;
            if (r != s) {
              gamma_iqsr_ptr = &gamma_iqrs[bf1 + ns*nbasis*(bf2 + nr*i)];
              for (q=0; q<nbasis; q++) {
                *gamma_iqsr_ptr += *c_sj * *gamma_irjq_ptr++;
                gamma_iqsr_ptr += ns;
                } // exit q loop
              }   // endif
            else gamma_irjq_ptr += nbasis;
            }     // exit bf2 loop
          }       // exit bf1 loop

        membuf_remote.release();

        }         // exit j loop
      }           // exit i loop

    // end of third quarter back-transformation
    // we now have gamma_iqrs (symmetrized)
    // for i-batch, all q, s in S, r in R
    timer->exit("3. q.b.t.");

    // only do this if integral is nonzero
    if (tbint->log2_shell_bound(R,S) >= tol) {

      // Compute contrib to Laj from (ov|vv) integrals
      // (done in AO basis to avoid generating (ov|vv);
      // here, generate Lpi for i-batch; later, transform
      // Lpi to get contribution to Laj
      timer->enter("(ov|vv) contrib to Laj");
      for (Q=0; Q<nshell; Q++) {
        nq = basis->shell(Q).nfunction();
        q_offset = basis->shell_to_function(Q);
        for (P=0; P<=Q; P++) {
          np = basis->shell(P).nfunction();
          p_offset = basis->shell_to_function(P);
       // if (scf_erep_bound(P,Q,R,S) < tol) {
       //   continue;  // skip ereps less than tol
       //   }
          if (tbint->log2_shell_bound(P,Q,R,S) < tol) {
            continue;  // skip ereps less than tol
            }
          timer->enter("erep");
          tbint->compute_shell(P,Q,R,S);
          timer->exit("erep");

          offset = nr*ns*nbasis;
          int_index = 0;

          for (bf1 = 0; bf1 < np; bf1++) {
            p = p_offset + bf1;
            for (bf2 = 0; bf2 < nq; bf2++) {
              q = q_offset + bf2;

              if (q < p) {
                int_index = ns*nr*(bf2+1 + nq*bf1);
                continue; // skip to next q value
                }

              for (bf3 = 0; bf3 < nr; bf3++) {
                r = r_offset + bf3;

                for (bf4 = 0; bf4 < ns; bf4++) {

                  if (fabs(intbuf[int_index]) > dtol) {
                    s = s_offset + bf4;

                    if (s < r) {
                      int_index++;
                      continue; // skip to next bf4 value
                      }

                    gamma_iqrs_ptr = &gamma_iqrs[bf4 + ns*(q + nbasis*bf3)];
                    gamma_iprs_ptr = &gamma_iqrs[bf4 + ns*(p + nbasis*bf3)];
                    pqrs = intbuf[int_index];

                    lpi_ptr = &Lpi[p*ni];
                    lqi_ptr = &Lpi[q*ni];

                    for (i=0; i<ni; i++) {
                      *lpi_ptr++ -= pqrs**gamma_iqrs_ptr;
                      if (p != q) {
                        *lqi_ptr++ -= pqrs**gamma_iprs_ptr;
                        }
                      gamma_iqrs_ptr += offset;
                      gamma_iprs_ptr += offset;
                      } // exit i loop
                    }   // endif

                  int_index++;
                  }     // exit bf4 loop
                }       // exit bf3 loop
              }         // exit bf2 loop
            }           // exit bf1 loop

          }             // exit P loop
        }               // exit Q loop
      timer->exit("(ov|vv) contrib to Laj");
      }                 // endif

    if (!dograd) continue;

    if (tbintder->log2_shell_bound(R,S) >= tol) {

      for (Q=0; Q<=S; Q++) {
        nq = basis->shell(Q).nfunction();
        q_offset = basis->shell_to_function(Q);

        for (P=0; P<=(Q==S ? R:Q); P++) {
          np = basis->shell(P).nfunction();
          p_offset = basis->shell_to_function(P);

          // If integral derivative is less than threshold skip to next P
          if (tbintder->log2_shell_bound(P,Q,R,S) < tol) continue;
          aointder_computed++;

          timer->enter("4. q.b.t.");
          bzerofast(gamma_pqrs,nfuncmax*nfuncmax*nfuncmax*nfuncmax);

          offset = nr*ns*nbasis;

          // Begin fourth quarter back-transformation
          gamma_pqrs_ptr = gamma_pqrs;
          for (bf1=0; bf1<np; bf1++) {
            p = bf1 + p_offset;
            for (bf2=0; bf2<nr; bf2++) {
              for (bf3=0; bf3<nq; bf3++) {
                q = bf3 + q_offset;
                for (bf4=0; bf4<ns; bf4++) {
                  c_pi = &scf_vector[p][i_offset];
                  c_qi = &scf_vector[q][i_offset];
                  gamma_iqrs_ptr = &gamma_iqrs[bf4 + ns*(q + nbasis*bf2)];
                  gamma_iprs_ptr = &gamma_iqrs[bf4 + ns*(p + nbasis*bf2)];
                  tmpval = 0.0;
                  for (i=0; i<ni; i++) {
                    tmpval += *c_pi * *gamma_iqrs_ptr;
                    if (p!=q) tmpval += *c_qi * *gamma_iprs_ptr;
                    c_pi++;
                    c_qi++;
                    gamma_iqrs_ptr += offset;
                    gamma_iprs_ptr += offset;
                    } // exit i loop
                  *gamma_pqrs_ptr += tmpval;
                  gamma_pqrs_ptr++;
                  }   // exit bf4 loop
                }     // exit bf3 loop
              }       // exit bf2 loop
            }         // exit bf1 loop
          // end of fourth quarter back-transformation
          timer->exit("4. q.b.t.");
          // (we now have the contribution from one i-batch to the
          // non-separable part of the 2PDM for one shell block PQRS)

          // Evaluate derivative integrals
          timer->enter("erep derivs");
          tbintder->compute_shell(P,Q,R,S,der_centers);
          timer->exit("erep derivs");

          // Compute contribution to gradient from non-sep 2PDM
          // (i.e., contract derivative integrals with gamma_pqrs)
          int_index = 0;
          timer->enter("non-sep 2PDM contrib.");
          for (int derset=0; derset<der_centers.n(); derset++) {
            for (int xyz=0; xyz<3; xyz++) {
              grad_ptr1 = &ginter[der_centers.atom(derset)][xyz];
              grad_ptr2 = &ginter[der_centers.omitted_atom()][xyz];
              for (bf1=0; bf1<np; bf1++) {
                p = bf1 + p_offset;
                for (bf2=0; bf2<nq; bf2++) {
                  q = bf2 + q_offset;
                  qp = q*(q+1)/2 + p;
                  for (bf3=0; bf3<nr; bf3++) {
                    r = bf3 + r_offset;
                    gamma_pqrs_ptr = &gamma_pqrs[ns*(bf2 + nq*(bf3 + nr*bf1))];
                    for (bf4=0; bf4<ns; bf4++) {
                      s = bf4 + s_offset;
                      sr = s*(s+1)/2 + r;
                      if (q == s && p == r) factor_pqrs = 1;
                      else factor_pqrs = 2;
                      tmpval = intderbuf[int_index]*factor_pqrs**gamma_pqrs_ptr;
                      gamma_pqrs_ptr++;
                      if (q>=p && s>=r && (P != R || Q != S || sr >= qp)) {
                        *grad_ptr1 += tmpval;
                        if (der_centers.has_omitted_center())
                          *grad_ptr2 -= tmpval;
                         }
                      int_index++;
                      } // exit bf4 loop
                    }   // exit bf3 loop
                  }     // exit bf2 loop
                }       // exit bf1 loop
              }         // exit xyz loop
            }           // exit derset loop
          timer->exit("non-sep 2PDM contrib.");

          } // exit P loop
        }   // exit Q loop
      }     // endif
    }       // end while

  delete[] gamma_iqrs;
  delete[] gamma_pqrs;

}

////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
