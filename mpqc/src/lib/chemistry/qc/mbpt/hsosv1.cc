//
// hsosv1.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Ida Nielsen <ida@kemi.aau.dk>
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

typedef int dmt_matrix;

#include <stdlib.h>
#include <math.h>

#include <util/misc/formio.h>
#include <util/misc/timer.h>
#include <util/class/class.h>
#include <util/state/state.h>
#include <util/group/message.h>
#include <math/scmat/matrix.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/scf/scf.h>
#include <chemistry/qc/mbpt/mbpt.h>
#include <chemistry/qc/mbpt/bzerofast.h>

void
MBPT2::compute_hsos_v1()
{
  int i, j, k;
  int s1, s2;
  int a, b;
  int isocc, asocc;   /* indices running over singly occupied orbitals */
  int nfuncmax = basis()->max_nfunction_in_shell();
  int nbasis = basis()->nbasis();
  int nvir;
  int nocc=0;
  int ndocc=0,nsocc=0;
  int i_offset; 
  int npass, pass;
  int ni;             /* batch size */
  int np, nq, nr, ns; 
  int P, Q, R, S;
  int p, q, r, s;
  int bf1, bf2, bf3, bf4;
  int index;
  int col_index;
  int docc_index, socc_index, vir_index;
  int me;
  int nproc;
  int rest;
  int a_rest;
  int a_number;              /* number of a-values processed by each node */
  int a_offset;
  int *a_vector;           /* each node's # of iajb integrals for one i,j */
  int shell_index;
  int compute_index;
  int tmp_index;
  int dim_ij;
  int aoint_computed = 0;
  int nshell;
  double A, B, C, ni_top, max, ni_double; /* Variables used to compute ni  */
  double *evals_open;    /* reordered scf eigenvalues                      */
  const double *intbuf;  /* 2-electron AO integral buffer                  */
  double *trans_int1;    /* partially transformed integrals                */
  double *trans_int2;    /* partially transformed integrals                */
  double *trans_int3;    /* partially transformed integrals                */
  double *trans_int4_node;/* each node's subset of fully transf. integrals */
  double *trans_int4;    /* fully transformed integrals                    */
  double *mo_int_do_so_vir; /*mo integral (is|sa); i:d.o.,s:s.o.,a:vir     */
  double *mo_int_tmp;    /* scratch array used in global summations        */
  double *socc_sum;      /* sum of 2-el integrals involving only s.o.'s    */
  double *iqrs, *iprs;
  double *iars_ptr, *iajs_ptr, *iajr_ptr;
  double iajr;
  double iars;
  double *iajb;
  double pqrs;
  double *c_qa;
  double *c_rb, *c_rj, *c_sj, *c_pi, *c_qi;
  double delta_ijab;
  double delta;
  double contrib1, contrib2;
  double ecorr_opt2=0,ecorr_opt1=0;
  double ecorr_zapt2;
  double ecorr_opt2_contrib=0, ecorr_zapt2_contrib=0;
  double escf;
  double eopt2,eopt1,ezapt2;
  double tol;     /* log2 of the erep tolerance (erep < 2^tol => discard) */

  me = msg_->me();
  
  cout << node0 << indent << "Just entered OPT2 program (opt2_v1)" << endl;

  tbint_ = integral()->electron_repulsion();
  intbuf = tbint_->buffer();

  tol = (int) (-10.0/log10(2.0));  /* discard ereps smaller than 10^-10 */

  nproc = msg_->n();
  cout << node0 << indent << "nproc = " << nproc << endl;

  ndocc = nsocc = 0;
  const double epsilon = 1.0e-4;
  for (i=0; i<nbasis; i++) {
    if      (reference_->occupation(i) >= 2.0 - epsilon) ndocc++;
    else if (reference_->occupation(i) >= 1.0 - epsilon) nsocc++;
    }

  /* do a few preliminary tests to make sure the desired calculation *
   * can be done (and appears to be meaningful!)                     */

  if (ndocc == 0 && nsocc == 0) {
    cerr << node0 << "There are no occupied orbitals; program exiting" << endl;
    abort();
    }

  if (nfzc > ndocc) {
    cerr << node0
         << "The number of frozen core orbitals exceeds the number" << endl
         << "of doubly occupied orbitals; program exiting" << endl;
    abort();
    }

  if (nfzv > nbasis - ndocc - nsocc) {
    cerr << node0
         << "The number of frozen virtual orbitals exceeds the number" << endl
         << "of unoccupied orbitals; program exiting" << endl;
    abort();
    }

  ndocc = ndocc - nfzc;
  /* nvir = # of unocc. orb. + # of s.o. orb. - # of frozen virt. orb. */
  nvir  = nbasis - ndocc - nfzc - nfzv; 
  /* nocc = # of d.o. orb. + # of s.o. orb - # of frozen d.o. orb. */
  nocc  = ndocc + nsocc;


  /* compute number of a-values (a_number) processed by each node */

  a_number = nvir/nproc; 
  a_rest = nvir%nproc;
  if (me < a_rest) a_number++;

  if (me == 0 && a_number < nsocc) { 
    cerr << "not enough memory allocated" << endl;
    /* must have all socc's on node 0 for computation of socc_sum*/
    abort();
    }

  if (me < a_rest) a_offset = me*a_number; /* a_offset for each node */
  else a_offset = a_rest*(a_number + 1) + (me - a_rest)*a_number;

  /* fill in elements of a_vector for gcollect */

  a_vector = (int*) malloc(nproc*sizeof(int));
  if (!a_vector) {
    cerr << "could not allocate storage for a_vector" << endl;
    abort();
    }
  for (i=0; i<nproc; i++) {
    a_vector[i] = nvir*(nvir/nproc)*sizeof(double);
    }
  for (i=0; i<a_rest; i++) {
    a_vector[i] += nvir*sizeof(double); /* first a_rest nodes hold an extra a */
    }


  /* compute batch size ni for opt2 loops                                 *
   * need to store the following arrays: trans_int1-4, trans_int4_node,   *
   * scf_vector, evals_open, socc_sum, mo_int_do_so_vir, mo_int_tmp and   *
   * a_vector;                                                            *
   * since a_number is not the same on all nodes, use node 0's a_number   *
   * (which is >= all other a_numbers) and broadcast ni afterwords        */

  if (me == 0) {
//    ni = (mem_alloc - sizeof(double)*(nvir*nvir + nbasis*(nvir+nocc)
//                                     + (nocc+nvir) + nsocc 
//                                     + 2*ndocc*nsocc*(nvir-nsocc) 
//                                     + nvir*a_number)
//                    - sizeof(int)*nproc)/
//          (sizeof(double)*(2*nfuncmax*nfuncmax*nbasis + nbasis*a_number*nocc));
    A = -0.5*sizeof(double)*nbasis*a_number;
    B = sizeof(double)*(2*nfuncmax*nfuncmax + (nocc+0.5)*a_number)*nbasis;
    C = sizeof(double)*(nvir*nvir + (nbasis+1)*(nocc+nvir) + nsocc 
                        + 2*ndocc*nsocc*(nvir-nsocc) + nvir*a_number)
        + sizeof(int)*nproc;
    ni_top = -B/(2*A);  /* ni value for which memory requirement is max. */
    max = A*ni_top*ni_top + B*ni_top + C;
    if (max <= mem_alloc) {
      ni = nocc;
      npass = 1;
      rest = 0;
      }
    else {
      ni_double = (-B + sqrt((double)(B*B - 4*A*(C-mem_alloc))))/(2*A);
      ni = (int) ni_double;
      if (ni > nocc) ni = nocc;
      rest = nocc%ni;
      npass = (nocc - rest)/ni + 1;
      if (rest == 0) npass--;
      }
    }
  msg_->bcast(ni);
  msg_->bcast(npass);
  msg_->bcast(rest);

  if (ni < nsocc) {
    cerr << "Not enough memory allocated" << endl;
    abort();
    }

  if (ni < 1) {     /* this applies only to a closed shell case */
    cerr << "Not enough memory allocated" << endl;
    abort();
    }

  cout << node0 << indent << "computed batchsize: " << ni << endl;

  nshell = basis()->nshell();
  if (me == 0) {
    cout << indent << " npass  rest  nbasis  nshell  nfuncmax"
         << "  ndocc  nsocc  nvir  nfzc  nfzv" << endl;
    cout << indent << scprintf("   %-4i   %-3i   %-5i   %-4i     %-3i"
                     "     %-3i     %-3i   %-3i    %-3i   %-3i",
                     npass,rest,nbasis,nshell,nfuncmax,ndocc,nsocc,nvir,nfzc,nfzv)
         << endl;
    cout << indent << scprintf("Using %i bytes of memory", mem_alloc) << endl;
    }

  /* the scf vector might be distributed between the nodes, but for OPT2 *
   * each node needs its own copy of the vector;                         *
   * therefore, put a copy of the scf vector on each node;               * 
   * while doing this, duplicate columns corresponding to singly         *
   * occupied orbitals and order columns as [socc docc socc unocc]       */
  /* also rearrange scf eigenvalues as [socc docc socc unocc]            *
   * want socc first to get the socc's in the first batch                *
   * (need socc's to compute energy denominators - see                   *
   * socc_sum comment below)                                             */

  evals_open = (double*) malloc((nbasis+nsocc-nfzc-nfzv)*sizeof(double));
  if (!evals_open) {
    cerr << "could not allocate storage for evals_open" << endl;
    abort();
    }

  RefDiagSCMatrix occ;
  RefDiagSCMatrix evals;
  RefSCMatrix Scf_Vec;
  eigen(evals, Scf_Vec, occ);

  if (debug_) {
    evals.print("eigenvalues");
    Scf_Vec.print("eigenvectors");
    }

  double *scf_vectort_dat = new double[nbasis*nbasis];
  Scf_Vec->convert(scf_vectort_dat);

  double** scf_vectort = new double*[nocc + nvir];

  int idoc = 0, ivir = 0, isoc = 0;
  for (i=nfzc; i<nbasis-nfzv; i++) {
    if (occ(i) >= 2.0 - epsilon) {
      evals_open[idoc+nsocc] = evals(i);
      scf_vectort[idoc+nsocc] = &scf_vectort_dat[i*nbasis];
      idoc++;
      }
    else if (occ(i) >= 1.0 - epsilon) {
      evals_open[isoc] = evals(i);
      scf_vectort[isoc] = &scf_vectort_dat[i*nbasis];
      evals_open[isoc+nocc] = evals(i);
      scf_vectort[isoc+nocc] = &scf_vectort_dat[i*nbasis];
      isoc++;
      }
    else {
      if (ivir < nvir) {
        evals_open[ivir+nocc+nsocc] = evals(i);
        scf_vectort[ivir+nocc+nsocc] = &scf_vectort_dat[i*nbasis];
        }
      ivir++;
      }
    }

  // need the transpose of the vector
  double **scf_vector = new double*[nbasis];
  double *scf_vector_dat = new double[(nocc+nvir)*nbasis];
  for (i=0; i<nbasis; i++) {
    scf_vector[i] = &scf_vector_dat[(nocc+nvir)*i];
    for (j=0; j<nocc+nvir; j++) {
      scf_vector[i][j] = scf_vectort[j][i];
      }
    }
  delete[] scf_vectort;
  delete[] scf_vectort_dat;

  if (debug_) {
    cout << node0 << indent << "Final eigenvalues and vectors" << endl;
    for (i=0; i<nocc+nvir; i++) {
      cout << node0 << indent << evals_open[i];
      for (j=0; j<nbasis; j++) {
        cout << node0 << " " << scf_vector[j][i];
        }
      cout << node0<< endl;
      }
    cout << node0 << endl;
    }

  /* allocate storage for integral arrays */

  dim_ij = nocc*ni - ni*(ni-1)/2;

  trans_int1 = (double*) malloc(nfuncmax*nfuncmax*nbasis*ni*sizeof(double));
  trans_int2 = (double*) malloc(nfuncmax*nfuncmax*nbasis*ni*sizeof(double));
  trans_int3 = (double*) malloc(nbasis*a_number*dim_ij*sizeof(double));
  trans_int4_node= (double*) malloc(nvir*a_number*sizeof(double));
  trans_int4 = (double*) malloc(nvir*nvir*sizeof(double));
  if (!(trans_int1 && trans_int2
        && (!a_number || trans_int3)
        && (!a_number || trans_int4_node) && trans_int4)){
    cerr << "could not allocate storage for integral arrays" << endl;
    abort();
    }
  if (nsocc) socc_sum  = (double*) malloc(nsocc*sizeof(double));
  if (nsocc) mo_int_do_so_vir = 
                   (double*) malloc(ndocc*nsocc*(nvir-nsocc)*sizeof(double));
  if (nsocc) mo_int_tmp = 
                   (double*) malloc(ndocc*nsocc*(nvir-nsocc)*sizeof(double));

  if (nsocc) bzerofast(mo_int_do_so_vir,ndocc*nsocc*(nvir-nsocc));

/**************************************************************************
*    begin opt2 loops                                                     *
***************************************************************************/

  for (pass=0; pass<npass; pass++) {
    i_offset= pass*ni;  
    if ((pass == npass - 1) && (rest != 0)) ni = rest;
    bzerofast(trans_int3,nbasis*a_number*dim_ij);

    shell_index = 0;

    tim_enter("RS loop");
    for (R = 0; R < basis()->nshell(); R++) {
      nr = basis()->shell(R).nfunction();

      for (S = 0; S <= R; S++) {
        ns = basis()->shell(S).nfunction();
        tim_enter("bzerofast trans_int1");
        bzerofast(trans_int1,nfuncmax*nfuncmax*nbasis*ni);
        tim_exit("bzerofast trans_int1");

        tim_enter("PQ loop");
        for (P = 0; P < basis()->nshell(); P++) {
          np = basis()->shell(P).nfunction();

          for (Q = 0; Q <= P; Q++) {
            shell_index++;
            if (shell_index%nproc != me) continue; 

            if (tbint_->log2_shell_bound(P,Q,R,S) < tol) {
              continue;                           /* skip ereps less than tol */
              }

            aoint_computed++;

            nq = basis()->shell(Q).nfunction();

            tim_enter("erep");
            tbint_->compute_shell(P,Q,R,S);
            tim_exit("erep");

            tim_enter("1. quart. tr."); 

            index = 0;

            for (bf1 = 0; bf1 < np; bf1++) {
              p = basis()->shell_to_function(P) + bf1;
 
              for (bf2 = 0; bf2 < nq; bf2++) {
                q = basis()->shell_to_function(Q) + bf2;
                if (q > p) {
                  /* if q > p: want to skip the loops over bf3-4  */
                  /* and larger bf2 values, so increment bf1 by 1 */
                  /* ("break") and adjust the value of index      */
                  index = (bf1 + 1) * nq * nr * ns;
                  break;
                  }

                for (bf3 = 0; bf3 < nr; bf3++) {

                  for (bf4 = 0; bf4 < ns; bf4++,index++) {
                    if (R==S && bf4>bf3) {
                      index = ((bf1*nq + bf2)*nr + (bf3+1))*ns;
                      break; 
                      }

                    if (fabs(intbuf[index])>1.0e-15) {
                      pqrs = intbuf[index];

                      iqrs = &trans_int1[((bf4*nr + bf3)*nbasis + q)*ni];
                      iprs = &trans_int1[((bf4*nr + bf3)*nbasis + p)*ni];
                    
                      if (p == q) pqrs *= 0.5;

                      col_index = i_offset;
                      c_pi = &scf_vector[p][col_index];
                      c_qi = &scf_vector[q][col_index];

                      for (i=ni; i; i--) {
                        *iqrs++ += pqrs * *c_pi++;
                        *iprs++ += pqrs * *c_qi++;
                        }
                      }
                    }   /* exit bf4 loop */
                  }     /* exit bf3 loop */
                }       /* exit bf2 loop */
              }         /* exit bf1 loop */
            tim_exit("1. quart. tr.");
            }           /* exit Q loop */
          }             /* exit P loop */

        tim_exit("PQ loop");

        tim_enter("sum int");
        msg_->sum(trans_int1,nr*ns*nbasis*ni,trans_int2);
        tim_exit("sum int");

        /* begin second quarter transformation */

        tim_enter("bzerofast trans_int2");
        bzerofast(trans_int2,nfuncmax*nfuncmax*nbasis*ni);
        tim_exit("bzerofast trans_int2");

        tim_enter("2. quart. tr.");

        for (bf3 = 0; bf3 < nr; bf3++) {

          for (bf4 = 0; bf4 < ns; bf4++) {
            if (R == S && bf4 > bf3) continue;

            for (q = 0; q < nbasis; q++) {
              c_qa = &scf_vector[q][nocc + a_offset];
              iqrs = &trans_int1[((bf4*nr + bf3)*nbasis + q)*ni];
              iars_ptr = &trans_int2[((bf4*nr + bf3)*a_number)*ni];

              for (a = 0; a < a_number; a++) {

                for (i=ni; i; i--) {
                  *iars_ptr++ += *c_qa * *iqrs++;
                  }

                iqrs -= ni;
                c_qa++;
                }
              }
            }
          }
        tim_exit("2. quart. tr.");

        /* begin third quarter transformation */
        tim_enter("3. quart. tr.");


        for (bf3 = 0; bf3<nr; bf3++) {
          r = basis()->shell_to_function(R) + bf3;

          for (bf4 = 0; bf4 <= (R == S ? bf3:(ns-1)); bf4++) {
            s = basis()->shell_to_function(S) + bf4;

            for (i=0; i<ni; i++) {
              tmp_index = i*(i+1)/2 + i*i_offset;

              for (a=0; a<a_number; a++) {
                iars = trans_int2[((bf4*nr + bf3)*a_number + a)*ni + i];
                if (r == s) iars *= 0.5;
                iajs_ptr = &trans_int3[tmp_index + dim_ij*(a + a_number*s)];
                iajr_ptr = &trans_int3[tmp_index + dim_ij*(a + a_number*r)];
                c_rj = scf_vector[r];
                c_sj = scf_vector[s];

                for (j=0; j<=i+i_offset; j++) {
                  *iajs_ptr++ += *c_rj++ * iars;
                  *iajr_ptr++ += *c_sj++ * iars;
                  }
                }
              }
            }     /* exit bf4 loop */
          }       /* exit bf3 loop */
        tim_exit("3. quart. tr.");
        }         /* exit S loop */
      }           /* exit R loop */
    tim_exit("RS loop");

    /* begin fourth quarter transformation;                                *
     * first tansform integrals with only s.o. indices;                    *
     * these integrals are needed to compute the denominators              *
     * in the various terms contributing to the correlation energy         *
     * and must all be computed in the first pass;                         *
     * the integrals are summed into the array socc_sum:                   *
     * socc_sum[isocc] = sum over asocc of (isocc asocc|asocc isocc)       *
     * (isocc, asocc = s.o. and the sum over asocc runs over all s.o.'s)   *
     * the individual integrals are not saved here, only the sums are kept */

    tim_enter("4. quart. tr.");
    if (pass == 0 && me == 0) {
      if (nsocc) bzerofast(socc_sum,nsocc);
      for (isocc=0; isocc<nsocc; isocc++) {

        for (r=0; r<nbasis; r++) {

          for (asocc=0; asocc<nsocc; asocc++) {
            socc_sum[isocc] += scf_vector[r][nocc+asocc]*
                                trans_int3[isocc*(isocc+1)/2 + isocc*i_offset 
                                          + isocc + dim_ij*(asocc + a_number*r)];
            }
          }
        }
      }

    tim_enter("bcast0 socc_sum");
    if (nsocc) msg_->bcast(socc_sum,nsocc);
    tim_exit("bcast0 socc_sum");

    tim_exit("4. quart. tr.");

    /* now we have all the sums of integrals involving s.o.'s (socc_sum);   *
     * begin fourth quarter transformation for all integrals (including     *
     * integrals with only s.o. indices); use restriction j <= (i_offset+i) *
     * to save flops                                                        */

    compute_index = 0;

    for (i=0; i<ni; i++) {

      for (j=0; j <= (i_offset+i); j++) {

       tim_enter("4. quart. tr.");

        bzerofast(trans_int4_node,nvir*a_number);

        for (r=0; r<nbasis; r++) {

          for (a=0; a<a_number; a++) {
            iajb = &trans_int4_node[a*nvir];
            c_rb = &scf_vector[r][nocc];
            iajr = trans_int3[i*(i+1)/2 + i*i_offset + j + dim_ij*(a+a_number*r)];

            for (b=0; b<nvir; b++) {
              *iajb++ += *c_rb++ * iajr;
              }
            }
          }

        tim_exit("4. quart. tr.");

        /* collect each node's part of fully transf. int. into trans_int4 */
        tim_enter("collect");
        msg_->collect(trans_int4_node,a_vector,trans_int4);
        tim_exit("collect");


        /* we now have the fully transformed integrals (ia|jb)      *
         * for one i, one j (j <= i_offset+i), and all a and b;     *
         * compute contribution to the OPT1 and OPT2 correlation    *
         * energies; use restriction b <= a to save flops           */

        tim_enter("compute ecorr");

        for (a=0; a<nvir; a++) {
          for (b=0; b<=a; b++) {
            compute_index++;
            if (compute_index%nproc != me) continue;

            docc_index = ((i_offset+i) >= nsocc && (i_offset+i) < nocc) 
                        + (j >= nsocc && j < nocc);
            socc_index = ((i_offset+i)<nsocc)+(j<nsocc)+(a<nsocc)+(b<nsocc);
            vir_index = (a >= nsocc) + (b >= nsocc);

            if (socc_index >= 3) continue; /* skip to next b value */
 
            delta_ijab = evals_open[i_offset+i] + evals_open[j] 
                       - evals_open[nocc+a] - evals_open[nocc+b];
            
            /* determine integral type and compute energy contribution */
            if (docc_index == 2 && vir_index == 2) {
              if (i_offset+i == j && a == b) {
                contrib1 = trans_int4[a*nvir + b]*trans_int4[a*nvir + b];
                ecorr_opt2 += contrib1/delta_ijab;
                ecorr_opt1 += contrib1/delta_ijab;
                }
              else if (i_offset+i == j || a == b) {
                contrib1 = trans_int4[a*nvir + b]*trans_int4[a*nvir + b];
                ecorr_opt2 += 2*contrib1/delta_ijab;
                ecorr_opt1 += 2*contrib1/delta_ijab;
                }
              else {
                contrib1 = trans_int4[a*nvir + b];
                contrib2 = trans_int4[b*nvir + a];
                ecorr_opt2 += 4*(contrib1*contrib1 + contrib2*contrib2
                             - contrib1*contrib2)/delta_ijab;
                ecorr_opt1 += 4*(contrib1*contrib1 + contrib2*contrib2
                             - contrib1*contrib2)/delta_ijab;
                }
              }
            else if (docc_index == 2 && socc_index == 2) {
              contrib1 = (trans_int4[a*nvir + b] - trans_int4[b*nvir + a])*
                         (trans_int4[a*nvir + b] - trans_int4[b*nvir + a]);
              ecorr_opt2 += contrib1/
                           (delta_ijab - 0.5*(socc_sum[a]+socc_sum[b]));
              ecorr_opt1 += contrib1/delta_ijab;
              }
            else if (socc_index == 2 && vir_index == 2) {
              contrib1 = (trans_int4[a*nvir + b] - trans_int4[b*nvir + a])*
                         (trans_int4[a*nvir + b] - trans_int4[b*nvir + a]);
              ecorr_opt2 += contrib1/
                          (delta_ijab - 0.5*(socc_sum[i_offset+i]+socc_sum[j]));
              ecorr_opt1 += contrib1/delta_ijab;
              }
            else if (docc_index == 2 && socc_index == 1 && vir_index == 1) {
              if (i_offset+i == j) {
                contrib1 = trans_int4[a*nvir + b]*trans_int4[a*nvir + b];
                ecorr_opt2 += contrib1/(delta_ijab - 0.5*socc_sum[b]);
                ecorr_opt1 += contrib1/delta_ijab;
                }
              else {
                contrib1 = trans_int4[a*nvir + b];
                contrib2 = trans_int4[b*nvir + a];
                ecorr_opt2 += 2*(contrib1*contrib1 + contrib2*contrib2 
                            - contrib1*contrib2)/(delta_ijab - 0.5*socc_sum[b]);
                ecorr_opt1 += 2*(contrib1*contrib1 + contrib2*contrib2 
                             - contrib1*contrib2)/delta_ijab;
                }
              }
            else if (docc_index == 1 && socc_index == 2 && vir_index == 1) {
              contrib1 = trans_int4[b*nvir+a]*trans_int4[b*nvir+a];
              if (j == b) {
                /* to compute the total energy contribution from an integral  *
                 * of the type (is1|s1a) (i=d.o., s1=s.o., a=unocc.), we need *
                 * the (is|sa) integrals for all s=s.o.; these integrals are  *
                 * therefore stored here in the array mo_int_do_so_vir, and   *
                 * the energy contribution is computed after exiting the loop *
                 * over i-batches (pass)                                      */
                mo_int_do_so_vir[a-nsocc + (nvir-nsocc)*
                                (i_offset+i-nsocc + ndocc*b)] =
                                trans_int4[b*nvir + a];
                ecorr_opt2_contrib += 1.5*contrib1/delta_ijab;
                ecorr_opt1         += 1.5*contrib1/delta_ijab;
                ecorr_zapt2_contrib += contrib1/
                               (delta_ijab - 0.5*(socc_sum[j]+socc_sum[b]))
                            + 0.5*contrib1/delta_ijab;
                }
              else {
                ecorr_opt2 += contrib1/
                             (delta_ijab - 0.5*(socc_sum[j] + socc_sum[b]));
                ecorr_opt1 += contrib1/delta_ijab;
                }
              }
            else if (docc_index == 1 && socc_index == 1 && vir_index == 2) {
              if (a == b) {
                contrib1 = trans_int4[a*nvir + b]*trans_int4[a*nvir + b];
                ecorr_opt2 += contrib1/(delta_ijab - 0.5*socc_sum[j]);
                ecorr_opt1 += contrib1/delta_ijab;
                }
              else {
                contrib1 = trans_int4[a*nvir + b];
                contrib2 = trans_int4[b*nvir + a];
                ecorr_opt2 += 2*(contrib1*contrib1 + contrib2*contrib2
                            - contrib1*contrib2)/(delta_ijab - 0.5*socc_sum[j]);
                ecorr_opt1 += 2*(contrib1*contrib1 + contrib2*contrib2
                             - contrib1*contrib2)/delta_ijab;
                }
              }
            }   /* exit b loop */
          }     /* exit a loop */
        tim_exit("compute ecorr");
        }       /* exit j loop */
      }         /* exit i loop */
    }           /* exit loop over i-batches (pass) */

  // don't need the AO integrals anymore
  tbint_ = 0;

  /* compute contribution from excitations of the type is1 -> s1a where   *
   * i=d.o., s1=s.o. and a=unocc; single excitations of the type i -> a,  *
   * where i and a have the same spin, contribute to this term;           *
   * (Brillouin's theorem not satisfied for ROHF wave functions);         */

  tim_enter("compute ecorr");

  if (nsocc > 0) {
    tim_enter("sum mo_int_do_so_vir");
    msg_->sum(mo_int_do_so_vir,ndocc*nsocc*(nvir-nsocc),mo_int_tmp);
    tim_exit("sum mo_int_do_so_vir");
    }

  /* add extra contribution for triplet and higher spin multiplicities *
   * contribution = sum over s1 and s2<s1 of (is1|s1a)*(is2|s2a)/delta */

  if (me == 0 && nsocc) {
    for (i=0; i<ndocc; i++) {

      for (a=0; a<nvir-nsocc; a++) {
        delta = evals_open[nsocc+i] - evals_open[nocc+nsocc+a];

        for (s1=0; s1<nsocc; s1++) {

          for (s2=0; s2<s1; s2++) {
            contrib1 = mo_int_do_so_vir[a + (nvir-nsocc)*(i + ndocc*s1)]*
                  mo_int_do_so_vir[a + (nvir-nsocc)*(i + ndocc*s2)]/delta;
            ecorr_opt2 += contrib1;
            ecorr_opt1 += contrib1;
            }
          }
        }     /* exit a loop */
      }       /* exit i loop */
    }

  tim_exit("compute ecorr");

  ecorr_zapt2 = ecorr_opt2 + ecorr_zapt2_contrib;
  ecorr_opt2 += ecorr_opt2_contrib;
  msg_->sum(ecorr_opt1);
  msg_->sum(ecorr_opt2);
  msg_->sum(ecorr_zapt2);
  msg_->sum(aoint_computed);

  escf = reference_->energy();
  hf_energy_ = escf;

  if (me == 0) {
    eopt2 = escf + ecorr_opt2;
    eopt1 = escf + ecorr_opt1;
    ezapt2 = escf + ecorr_zapt2;

    /* print out various energies etc.*/

    cout << indent
         << "Number of shell quartets for which AO integrals would" << endl
         << indent
         << "have been computed without bounds checking: "
         << npass*nshell*nshell*(nshell+1)*(nshell+1)/4 << endl;
    cout << indent
         << "Number of shell quartets for which AO integrals" << endl
         << indent << "were computed: " << aoint_computed << endl;
    cout << indent
         << scprintf("ROHF energy [au]:                  %17.12lf\n", escf);
    cout << indent
         << scprintf("OPT1 energy [au]:                  %17.12lf\n", eopt1);
    cout << indent
         << scprintf("OPT2 second order correction [au]: %17.12lf\n", ecorr_opt2);
    cout << indent
         << scprintf("OPT2 energy [au]:                  %17.12lf\n", eopt2);
    cout << indent
         << scprintf("ZAPT2 correlation energy [au]:     %17.12lf\n", ecorr_zapt2);
    cout << indent
         << scprintf("ZAPT2 energy [au]:                 %17.12lf\n", ezapt2);
    }
  msg_->bcast(eopt1);
  msg_->bcast(eopt2);
  msg_->bcast(ezapt2);

  if (method_ && !strcmp(method_,"opt1")) {
    set_energy(eopt1);
    set_actual_value_accuracy(reference_->actual_value_accuracy()
                              *ref_to_mp2_acc);
    }
  else if (method_ && !strcmp(method_,"opt2")) {
    set_energy(eopt2);
    set_actual_value_accuracy(reference_->actual_value_accuracy()
                              *ref_to_mp2_acc);
    }
  else {
    if (!(!method_ || !strcmp(method_,"zapt"))) {
      cout << node0 << indent
           << "MBPT2: bad method for closed shell case: " << method_
           << ", using zapt" << endl;
      }
    set_energy(ezapt2);
    set_actual_value_accuracy(reference_->actual_value_accuracy()
                              *ref_to_mp2_acc);
    }

  free(trans_int1);
  free(trans_int2);
  free(trans_int3);
  free(trans_int4_node);
  free(trans_int4);
  free(a_vector);
  if (nsocc) free(socc_sum);
  if (nsocc) free(mo_int_do_so_vir);
  if (nsocc) free(mo_int_tmp);
  free(evals_open);

  delete[] scf_vector;
  delete[] scf_vector_dat;
  }

////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
