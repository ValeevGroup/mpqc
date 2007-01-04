//
// csgrade12.cc
// based on: csgrad.cc
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

#ifdef __GNUC__
#pragma implementation
#endif

#include <math.h>

#include <util/misc/formio.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/mbpt/bzerofast.h>
#include <chemistry/qc/mbpt/csgrade12.h>
#include <chemistry/qc/basis/distshpair.h>

#include <chemistry/qc/mbpt/util.h>

using namespace std;
using namespace sc;

extern BiggestContribs biggest_ints_1;

#define PRINT1Q 0

CSGradErep12Qtr::CSGradErep12Qtr(int mythread_a, int nthread_a,
                                 int me_a, int nproc_a,
                                 const Ref<MemoryGrp> &mem_a,
                                 const Ref<MessageGrp> &msg_a,
                                 const Ref<ThreadLock> &lock_a,
                                 const Ref<GaussianBasisSet> &basis_a,
                                 const Ref<TwoBodyInt> &tbint_a,
                                 int nocc_a,
                                 double **scf_vector_a,
                                 double tol_a, int debug_a,
                                 int dynamic_a, double print_percent_a,
                                 DistShellPair::SharedData *shellpair_shared_data,
                                 int usep4):
  shellpair_shared_data_(shellpair_shared_data)
{
  msg = msg_a;
  mythread = mythread_a;
  nthread = nthread_a;
  lock = lock_a;
  basis = basis_a;
  tbint = tbint_a;
  nocc = nocc_a;
  me = me_a;
  nproc = nproc_a;
  tol = tol_a;
  mem = mem_a;
  scf_vector = scf_vector_a;
  debug = debug_a;
  dynamic_ = dynamic_a;
  print_percent_ = print_percent_a;
  usep4_ = usep4;

  aoint_computed = 0;
  timer = new RegionTimer();
}

CSGradErep12Qtr::~CSGradErep12Qtr()
{
}

void
CSGradErep12Qtr::run()
{
  int P,Q,R,S;
  int p,q,r,s;
  int np,nq,nr,ns;
  int bf1,bf2,bf3,bf4;
  int p_offset,q_offset,r_offset,s_offset;
  int offset;
  int nfuncmax = basis->max_nfunction_in_shell();
  int nshell = basis->nshell();
  int nbasis = basis->nbasis();
  double dtol = pow(2.0,tol);
  double *iqjs_ptr;
  double *iqrs_ptr, *iprs_ptr;
  double *c_pi, *c_qi;
  double tmpval;
  int i,j;
  double *iqjs_contrib;  // local contributions to integral_iqjs
  double *iqjr_contrib;  // local contributions to integral_iqjr

  const double *intbuf = tbint->buffer();

  iqjs_contrib  = mem->malloc_local_double(nbasis*nfuncmax);
  iqjr_contrib  = mem->malloc_local_double(nbasis*nfuncmax);

  double *integral_iqrs; // quarter transformed two-el integrals
  lock->lock();
  integral_iqrs = new double[ni*nbasis*nfuncmax*nfuncmax];
  lock->unlock();

  int work_per_thread = ((nshell*(nshell+1))/2)/(nproc*nthread);
  int print_interval = work_per_thread/100;
  int time_interval = work_per_thread/10;
  int print_index = 0;
  if (print_interval == 0) print_interval = 1;
  if (time_interval == 0) time_interval = 1;
  if (work_per_thread == 0) work_per_thread = 1;

  if (debug) {
    lock->lock();
    ExEnv::outn() << scprintf("%d:%d: starting get_task loop",me,mythread) << endl;
    lock->unlock();
    }

  // Use petite list for symmetry utilization
  Ref<PetiteList> p4list = tbint->integral()->petite_list();

  DistShellPair shellpairs(msg,nthread,mythread,lock,basis,basis,dynamic_,
                           shellpair_shared_data_);
  shellpairs.set_print_percent(print_percent_);
  shellpairs.set_debug(debug);
  if (debug) shellpairs.set_print_percent(1);
  S = 0;
  R = 0;
  Timer tim(timer);
  while (shellpairs.get_task(S,R)) {
    ns = basis->shell(S).nfunction();
    s_offset = basis->shell_to_function(S);

    nr = basis->shell(R).nfunction();
    r_offset = basis->shell_to_function(R);

    if (debug > 1 && (print_index++)%print_interval == 0) {
      lock->lock();
      ExEnv::outn() << scprintf("%d:%d: (PQ|%d %d) %d%%",
                       me,mythread,R,S,(100*print_index)/work_per_thread)
           << endl;
      lock->unlock();
      }
    if (debug > 1 && (print_index)%time_interval == 0) {
      lock->lock();
      ExEnv::outn() << scprintf("timer for %d:%d:",me,mythread) << endl;
      tim.print();
      lock->unlock();
      }

    bzerofast(integral_iqrs, ni*nbasis*nfuncmax*nfuncmax);

    for (Q=0; Q<nshell; Q++) {
      nq = basis->shell(Q).nfunction();
      q_offset = basis->shell_to_function(Q);
      for (P=0; P<=Q; P++) {
        np = basis->shell(P).nfunction();
        p_offset = basis->shell_to_function(P);

	// check if symmetry unique and compute degeneracy
        int deg;
        if (usep4_) deg = p4list->in_p4(P,Q,R,S);
        else deg = 1;
        double symfac = (double) deg;
        if (deg == 0)
          continue;

        if (tbint->log2_shell_bound(P,Q,R,S) < tol) {
          continue;  // skip ereps less than tol
          }

        aoint_computed++;

        tim.enter("erep");
        tbint->compute_shell(P,Q,R,S);
        tim.exit("erep");

        tim.enter("1. q.t.");
        // Begin first quarter transformation;
        // generate (iq|rs) for i active

        offset = nr*ns*nbasis;
        const double *pqrs_ptr = intbuf;
        for (bf1 = 0; bf1 < np; bf1++) {
          p = p_offset + bf1;
          for (bf2 = 0; bf2 < nq; bf2++) {
            q = q_offset + bf2;

            if (q < p) {
              pqrs_ptr = &intbuf[ns*nr*(bf2+1 + nq*bf1)];
              continue; // skip to next q value
              }

            for (bf3 = 0; bf3 < nr; bf3++) {
              r = r_offset + bf3;

              for (bf4 = 0; bf4 < ns; bf4++) {
                s = s_offset + bf4;

                if (s < r) {
                  pqrs_ptr++;
                  continue; // skip to next bf4 value
                  }

                if (fabs(*pqrs_ptr) > dtol) {
                  iprs_ptr = &integral_iqrs[bf4 + ns*(p + nbasis*bf3)];
                  iqrs_ptr = &integral_iqrs[bf4 + ns*(q + nbasis*bf3)];
                  c_qi = &scf_vector[q][i_offset];
                  c_pi = &scf_vector[p][i_offset];
                  tmpval = *pqrs_ptr;
		  // multiply each integral by its symmetry degeneracy factor
		  tmpval *= symfac;
                  for (i=0; i<ni; i++) {
                    *iprs_ptr += *c_qi++*tmpval;
                    iprs_ptr += offset;
                    if (p != q) {
                      *iqrs_ptr += *c_pi++*tmpval;
                      iqrs_ptr += offset;
                      }
                    } // exit i loop
                  }   // endif

                pqrs_ptr++;
                } // exit bf4 loop
              }   // exit bf3 loop
            }     // exit bf2 loop
          }       // exit bf1 loop
        // end of first quarter transformation
        tim.exit("1. q.t.");

        }           // exit P loop
      }             // exit Q loop

#if PRINT1Q
      {
      lock->lock();
      double *tmp = integral_iqrs;
      for (int i = 0; i<ni; i++) {
        for (int r = 0; r<nr; r++) {
          for (int q = 0; q<nbasis; q++) {
            for (int s = 0; s<ns; s++) {
              printf("1Q: (%d %d|%d %d) = %12.8f\n",
                     i,q,r+r_offset,s+s_offset,*tmp);
              tmp++;
              }
            }
          }
        }
      lock->unlock();
      }
#endif
#if PRINT_BIGGEST_INTS
      {
      lock->lock();
      double *tmp = integral_iqrs;
      for (int i = 0; i<ni; i++) {
        for (int r = 0; r<nr; r++) {
          for (int q = 0; q<nbasis; q++) {
            for (int s = 0; s<ns; s++) {
              if (i+i_offset==104) {
                biggest_ints_1.insert(*tmp,i+i_offset,q,r+r_offset,s+s_offset);
                }
              tmp++;
              }
            }
          }
        }
      lock->unlock();
      }
#endif

    tim.enter("2. q.t.");
    // Begin second quarter transformation;
    // generate (iq|jr) for i active and j active or frozen
    for (i=0; i<ni; i++) {
      for (j=0; j<nocc; j++) {

        bzerofast(iqjs_contrib, nbasis*nfuncmax);
        bzerofast(iqjr_contrib, nbasis*nfuncmax);

        for (bf1=0; bf1<ns; bf1++) {
          s = s_offset + bf1;
          double *c_sj = &scf_vector[s][j];
          double *iqjr_ptr = iqjr_contrib;
          for (bf2=0; bf2<nr; bf2++) {
            r = r_offset + bf2;
            if (r > s) {
              break; // skip to next bf1 value
              }
            double c_rj = scf_vector[r][j];
            iqjs_ptr = &iqjs_contrib[bf1*nbasis];
            iqrs_ptr = &integral_iqrs[bf1 + ns*nbasis*(bf2 + nr*i)];
            for (q=0; q<nbasis; q++) {
              *iqjs_ptr++ += c_rj * *iqrs_ptr;
              if (r != s) *iqjr_ptr += *c_sj * *iqrs_ptr;
              iqjr_ptr++;
              iqrs_ptr += ns;
              } // exit q loop
            }   // exit bf2 loop
          }     // exit bf1 loop

        // We now have contributions to iqjs and iqjr for one pair i,j,
        // all q, r in R and s in S; send iqjs and iqjr to the node
        // (ij_proc) which is going to have this ij pair
        int ij_proc =  (i*nocc + j)%nproc;
        int ij_index = (i*nocc + j)/nproc;

        // Sum the iqjs_contrib to the appropriate place
        size_t ij_offset = size_t(nbasis)*(s_offset + size_t(nbasis)*ij_index);
        mem->sum_reduction_on_node(iqjs_contrib,
                                   ij_offset, ns*nbasis, ij_proc);

        ij_offset = size_t(nbasis)*(r_offset + size_t(nbasis)*ij_index);
        mem->sum_reduction_on_node(iqjr_contrib,
                                   ij_offset, nr*nbasis, ij_proc);

        }     // exit j loop
      }       // exit i loop
    // end of second quarter transformation
    tim.exit("2. q.t.");

    }         // exit while get_task

  if (debug) {
    lock->lock();
    ExEnv::outn() << scprintf("%d:%d: done with get_task loop",me,mythread) << endl;
    lock->unlock();
    }

  lock->lock();
  delete[] integral_iqrs;
  mem->free_local_double(iqjs_contrib);
  mem->free_local_double(iqjr_contrib);
  lock->unlock();
}

////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
