//
// transform_13inds.cc
//
// Copyright (C) 2001 Edward Valeev
//
// Author: Edward Valeev <edward.valeev@chemistry.gatech.edu>
// Maintainer: EV
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
#include <stdexcept>

#include <util/misc/formio.h>
#include <util/misc/regtime.h>
#include <chemistry/qc/basis/gpetite.h>
#include <chemistry/qc/mbpt/bzerofast.h>
#include <chemistry/qc/mbpt/util.h>
#include <chemistry/qc/basis/distshpair.h>
#include <chemistry/qc/mbptr12/blas.h>
#include <chemistry/qc/mbptr12/transform_13inds.h>

using namespace std;
using namespace sc;

#define PRINT1Q 0
#define PRINT2Q 0
#define PRINT_NUM_TE_TYPES 1

// The FAST_BUT_WRONG flags is useful for exercising the communications
// layer.  It causes the first and second quarter transformation to be
// omitted, but all communication is still performed.  This permits
// problems in communications libraries to be more quickly identified.
#define FAST_BUT_WRONG 0

TwoBodyMOIntsTransform_13Inds::TwoBodyMOIntsTransform_13Inds(
  const Ref<TwoBodyMOIntsTransform>& tform, int mythread, int nthread,
  const Ref<ThreadLock>& lock, const Ref<TwoBodyInt> &tbint, double tol, int debug) :
  tform_(tform), mythread_(mythread), nthread_(nthread), lock_(lock), tbint_(tbint),
  tol_(tol), debug_(debug)
{
  timer_ = new RegionTimer();
  aoint_computed_ = 0;
  ni_ = tform_->batchsize();
  i_offset_ = 0;
}

TwoBodyMOIntsTransform_13Inds::~TwoBodyMOIntsTransform_13Inds()
{
  timer_ = 0;
}

/*
  Distribute work by SR
  
   for all PQ
    compute unique (PQ|RS)
    transform to (IM|RS) stored as rsim where M are all AOs for basis set 2
   end PQ

   transform (IM|RS) to (IM|JS) stored as ijsm and accumulate to the tasks that holds respective ij-pairs.

  end SR
*/

void
TwoBodyMOIntsTransform_13Inds::run()
{
  Timer tim(timer_);
  Ref<MemoryGrp> mem = tform_->mem();
  Ref<MessageGrp> msg = tform_->msg();
  Ref<R12IntsAcc> ints_acc = tform_->ints_acc();
  const int me = msg->me();
  const int nproc = msg->n();
  Ref<MOIndexSpace> space1 = tform_->space1();
  Ref<MOIndexSpace> space2 = tform_->space2();
  Ref<MOIndexSpace> space3 = tform_->space3();
  Ref<MOIndexSpace> space4 = tform_->space4();

  Ref<GaussianBasisSet> bs1 = space1->basis();
  Ref<GaussianBasisSet> bs2 = space2->basis();
  Ref<GaussianBasisSet> bs3 = space3->basis();
  Ref<GaussianBasisSet> bs4 = space4->basis();
  const bool bs1_eq_bs2 = (bs1 == bs2);
  const bool bs3_eq_bs4 = (bs3 == bs4);

  const bool dynamic = tform_->dynamic();
  const double print_percent = tform_->print_percent();

  const int ni = ni_;
  const int rank1 = space1->rank();
  const int rank3 = space3->rank();
  const int nfuncmax1 = bs1->max_nfunction_in_shell();
  const int nfuncmax2 = bs2->max_nfunction_in_shell();
  const int nfuncmax3 = bs3->max_nfunction_in_shell();
  const int nfuncmax4 = bs4->max_nfunction_in_shell();
  const int nsh1 = bs1->nshell();
  const int nsh2 = bs2->nshell();
  const int nsh3 = bs3->nshell();
  const int nsh4 = bs4->nshell();
  const int nbasis1 = bs1->nbasis();
  const int nbasis2 = bs2->nbasis();
  const int nbasis3 = bs3->nbasis();
  const int nbasis4 = bs4->nbasis();
  double dtol = pow(2.0,tol_);
  const size_t memgrp_blksize = tform_->memgrp_blksize()/sizeof(double);

  const int tbtype_anti1 = TwoBodyInt::r12t1;
  const int tbtype_anti2 = TwoBodyInt::r12t2;

  double** vector1 = new double*[nbasis1];
  double** vector3 = new double*[nbasis3];
  vector1[0] = new double[rank1*nbasis1];
  vector3[0] = new double[rank3*nbasis3];
  for(int i=1; i<nbasis1; i++) vector1[i] = vector1[i-1] + rank1;
  for(int i=1; i<nbasis3; i++) vector3[i] = vector3[i-1] + rank3;
  space1->coefs().convert(vector1);
  space3->coefs().convert(vector3);

  /*-------------------------------------------------------------
    Find integrals buffers to 1/r12, r12, and [r12,T1] integrals
   -------------------------------------------------------------*/
  const int num_te_types = tform_->num_te_types();
  const double *intbuf[TwoBodyInt::max_num_tbint_types];
  intbuf[TwoBodyInt::eri] = tbint_->buffer(TwoBodyInt::eri);
  intbuf[TwoBodyInt::r12] = tbint_->buffer(TwoBodyInt::r12);
  intbuf[TwoBodyInt::r12t1] = tbint_->buffer(TwoBodyInt::r12t1);
  intbuf[TwoBodyInt::r12t2] = tbint_->buffer(TwoBodyInt::r12t2);

  /*-----------------------------------------------------
    Allocate buffers for partially transformed integrals
   -----------------------------------------------------*/
  double *ijsq_contrib;  // local contributions to integral_ijsq
  double *ijrq_contrib;  // local contributions to integral_ijrq
  double **rsiq_ints = new double*[num_te_types];     // quarter-transformed integrals for each RS pair
  for(int te_type=0;te_type<num_te_types;te_type++) {
    rsiq_ints[te_type] = new double[ni*nbasis2*nfuncmax3*nfuncmax4];
  }
  ijsq_contrib  = mem->malloc_local_double(nbasis2*nfuncmax4);
  if (bs3_eq_bs4)
    ijrq_contrib  = mem->malloc_local_double(nbasis2*nfuncmax4);
  else
    ijrq_contrib  = NULL;
  
  /*-----------------------------
    Initialize work distribution
   -----------------------------*/
  sc::DistShellPair shellpairs(msg,nthread_,mythread_,lock_,bs4,bs3,dynamic,
                               tform_->shell_pair_data());
  shellpairs.set_debug(debug_);
  if (debug_) shellpairs.set_print_percent(print_percent/10.0);
  else shellpairs.set_print_percent(print_percent);
  int work_per_thread = bs3_eq_bs4 ? 
    ((nsh3*(nsh3+1))/2)/(nproc*nthread_) :
    (nsh3*nsh4)/(nproc*nthread_) ;
  int print_interval = work_per_thread/100;
  int time_interval = work_per_thread/10;
  int print_index = 0;
  if (print_interval == 0) print_interval = 1;
  if (time_interval == 0) time_interval = 1;
  if (work_per_thread == 0) work_per_thread = 1;

  if (debug_) {
    lock_->lock();
    ExEnv::outn() << scprintf("%d:%d: starting get_task loop",me,mythread_) << endl;
    lock_->unlock();
  }

  Ref<GenPetite4> p4list
    = construct_gpetite(bs1,bs2,bs3,bs4);

#if FAST_BUT_WRONG
  for(int te_type=0;te_type<num_te_types;te_type++) {
    bzerofast(rsiq_ints[te_type], ni*nbasis2*nfuncmax3*nfuncmax4);
    bzerofast(ijsq_contrib[te_type], nbasis2*nfuncmax4);
    if (bs3_eq_bs4)
      bzerofast(ijrq_contrib[te_type], nbasis2*nfuncmax4);
    }
#endif

  int R = 0;
  int S = 0;
  while (shellpairs.get_task(S,R)) {
    // if bs3_eq_bs4 then S >= R always (see sc::exp::DistShellPair)
    int nr = bs3->shell(R).nfunction();
    int r_offset = bs3->shell_to_function(R);
    
    int ns = bs4->shell(S).nfunction();
    int s_offset = bs4->shell_to_function(S);
    
    int nrs = nr*ns;

    if (debug_ > 1 && (print_index++)%print_interval == 0) {
      lock_->lock();
      ExEnv::outn() << scprintf("%d:%d: (PQ|%d %d) %d%%",
			       me,mythread_,R,S,(100*print_index)/work_per_thread)
		   << endl;
      lock_->unlock();
    }
    if (debug_ > 1 && (print_index)%time_interval == 0) {
      lock_->lock();
      ExEnv::outn() << scprintf("timer for %d:%d:",me,mythread_) << endl;
      tim.print();
      lock_->unlock();
    }

#if !FAST_BUT_WRONG
    // Zero out 1 q.t. storage
    for(int te_type=0;te_type<num_te_types;te_type++)
      bzerofast(rsiq_ints[te_type], nrs*ni*nbasis2);

    for (int P=0; P<nsh1; P++) {
      int np = bs1->shell(P).nfunction();
      int p_offset = bs1->shell_to_function(P);

      int Qmax = (bs1_eq_bs2 ? P : nsh2-1);
      for (int Q=0; Q<=Qmax; Q++) {
	int nq = bs2->shell(Q).nfunction();
	int q_offset = bs2->shell_to_function(Q);
        
	// check if symmetry unique and compute degeneracy
	int deg = p4list->in_p4(P,Q,R,S);
	if (deg == 0)
	  continue;
	double symfac = (double) deg;

        if (tbint_->log2_shell_bound(P,Q,R,S) < tol_) {
          continue;  // skip shell quartets less than tol
	}

        aoint_computed_++;

        tim.enter("AO integrals");
        tbint_->compute_shell(P,Q,R,S);
        tim.exit("AO integrals");

        tim.enter("1. q.t.");

        // Begin first quarter transformation;
        // generate (iq|rs) for i active
        // if bs1_eq_bs2 then (ip|rs) are also generated
        // store the integrals as rsiq
	for(int te_type=0; te_type<num_te_types; te_type++) {
	  const double *pqrs_ptr = intbuf[te_type];

	  for (int bf1 = 0; bf1 < np; bf1++) {
	    int p = p_offset + bf1;
            int qmax = (bs1_eq_bs2 && P == Q) ? bf1 : nq-1;

	    for (int bf2 = 0; bf2 <= qmax; bf2++) {
	      int q = q_offset + bf2;

	      for (int bf3 = 0; bf3 < nr; bf3++) {
                int smin = (bs3_eq_bs4 && R == S) ? bf3 : 0;
                pqrs_ptr += smin;

		for (int bf4 = smin; bf4 <ns; bf4++) {

                  // Only transform integrals larger than the threshold
		  if (fabs(*pqrs_ptr) > dtol) {

		    double* rsiq_ptr = &rsiq_ints[te_type][q + nbasis2*(0 + ni*(bf4 + ns*bf3))];
		    const double* c_pi = vector1[p] + i_offset_;

                    double* rsip_ptr;
		    const double* c_qi;
                    if (bs1_eq_bs2) {
		      rsip_ptr = &rsiq_ints[te_type][p + nbasis2*(0 + ni*(bf4 + ns*bf3))];
		      c_qi = vector1[q] + i_offset_;
                    }
                    
		    double rsiq_int_contrib = *pqrs_ptr;
		    // multiply each integral by its symmetry degeneracy factor
		    rsiq_int_contrib *= symfac;
                    
                    if (bs1_eq_bs2) {

                      double rsip_int_contrib = rsiq_int_contrib;
                      if (te_type == tbtype_anti1)
                        rsip_int_contrib = -1.0*rsiq_int_contrib;

                      if (p == q) {
                        for (int i=0; i<ni; i++) {
                          *rsiq_ptr += *c_pi++ * rsiq_int_contrib;
                          rsiq_ptr += nbasis2;
                        }
                      }
                      else {
                        // p != q
                        for (int i=0; i<ni; i++) {
                          *rsip_ptr += *c_qi++ * rsip_int_contrib;
                          rsip_ptr += nbasis2;
                          *rsiq_ptr += *c_pi++ * rsiq_int_contrib;
                          rsiq_ptr += nbasis2;
                        }
                      }
                        
                    }
                    else {

                      for (int i=0; i<ni; i++) {
                        *rsiq_ptr += *c_pi++ * rsiq_int_contrib;
                        rsiq_ptr += nbasis2;
                      }

                    } // endif bs1_eq_bs2
                  }   // endif dtol

		  pqrs_ptr++;
                } // exit bf4 loop
              }   // exit bf3 loop
            }     // exit bf2 loop
            pqrs_ptr += (nq - qmax - 1) * nrs;
          }       // exit bf1 loop
	  // end of first quarter transformation
	}
	tim.exit("1. q.t.");

        }           // exit P loop
      }             // exit Q loop
#endif // !FAST_BUT_WRONG

#if PRINT1Q
    {
        lock_->lock();
        for(int te_type=0; te_type<PRINT_NUM_TE_TYPES; te_type++) {
          for (int i = 0; i<ni; i++) {
            for (int q = 0; q<nbasis2; q++) {
              for (int r = 0; r<nr; r++) {
                int rr = r+r_offset;
                for (int s = 0; s<ns; s++) {
                  int ss = s+s_offset;
                  double value = rsiq_ints[te_type][q+nbasis2*(i+ni*(s+ns*r))];
                  printf("1Q: type = %d (%d %d|%d %d) = %12.8f\n",
                         te_type,i+i_offset_,q,rr,ss,value);
                }
              }
            }
          }
        }
        lock_->unlock();
    }
#endif

    const int niq = ni*nbasis2;

    tim.enter("2. q.t.");
    // Begin second quarter transformation;
    // generate (iq|js) stored as ijsq (also generate (iq|jr), if needed)

    for(int te_type=0; te_type<num_te_types; te_type++) {
      for (int i=0; i<ni; i++) {
        for (int j=0; j<rank3; j++) {

#if !FAST_BUT_WRONG
          bzerofast(ijsq_contrib, nbasis2*ns);
          if (bs3_eq_bs4)
            bzerofast(ijrq_contrib, nbasis2*nr);

          int ij_proc =  (i*rank3 + j)%nproc;
          int ij_index = (i*rank3 + j)/nproc;
          const size_t ijsq_start = (size_t)(num_te_types*ij_index + te_type) * memgrp_blksize;
          const double *rsiq_ptr = rsiq_ints[te_type] + i*nbasis2;

          if (bs3_eq_bs4) {

            const double ket_perm_pfac = (te_type == tbtype_anti2) ? -1.0 : 1.0;

            for (int bf3 = 0; bf3 < nr; bf3++) {
              int r = r_offset + bf3;
              int smin = (bs3_eq_bs4 && R == S) ? bf3 : 0;
              rsiq_ptr += smin*niq;

              for (int bf4 = smin; bf4 <ns; bf4++, rsiq_ptr+=niq) {
                int s = s_offset + bf4;

                // second quarter transform
                // rs = js
                // rs = jr
                
                double* ijsq_ptr = ijsq_contrib + bf4*nbasis2;
                double* ijrq_ptr = ijrq_contrib + bf3*nbasis2;
                const double* i_ptr = rsiq_ptr;

                const double c_rj = vector3[r][j];
                const double c_sj = vector3[s][j];

                if (r != s) {
                  for (int q=0; q<nbasis2; q++) {

                    double value = *i_ptr++;
                    *ijsq_ptr++ += c_rj * value;
                    *ijrq_ptr++ += c_sj * value * ket_perm_pfac;

                  }
                }
                else {
                  for (int q=0; q<nbasis2; q++) {

                    double value = *i_ptr++;
                    *ijsq_ptr++ += c_rj * value;

                  }
                }
              }
            }
          }
          else {

            for (int bf3 = 0; bf3 < nr; bf3++) {
              int r= r_offset + bf3;
              double* ijsq_ptr = ijsq_contrib;

              for (int bf4 = 0; bf4 <ns; bf4++, rsiq_ptr+=niq) {

                // second quarter transform
                // rs = js
                const double* i_ptr = rsiq_ptr;

                const double c_rj = vector3[r][j];

                for (int q=0; q<nbasis2; q++) {

                  double value = *i_ptr++;
                  *ijsq_ptr++ += c_rj * value;

                }
              }
            }
          }

          // We now have contributions to ijsq (and ijrq) for one pair i,j,
          // all q, and s in S (r in R); send ijsq (and ijrq) to the node
          // (ij_proc) which is going to have this ij pair
#endif // !FAST_BUT_WRONG

          // Sum the ijsq_contrib to the appropriate place
          size_t ij_offset = (size_t)nbasis2*s_offset + ijsq_start;
          mem->sum_reduction_on_node(ijsq_contrib,
                                     ij_offset, ns*nbasis2, ij_proc);

          if (bs3_eq_bs4) {
            size_t ij_offset = (size_t)nbasis2*r_offset + ijsq_start;
            mem->sum_reduction_on_node(ijrq_contrib,
                                       ij_offset, nr*nbasis2, ij_proc);
          }

        } // endif j
      } // endif i
    }  // endif te_type
    tim.exit("2. q.t.");
	  
  }         // exit while get_task

  if (debug_) {
    lock_->lock();
    ExEnv::outn() << scprintf("%d:%d: done with get_task loop",me,mythread_) << endl;
    lock_->unlock();
  }

  for(int te_type=0; te_type<num_te_types; te_type++) {
    delete[] rsiq_ints[te_type];
  }
  delete[] rsiq_ints;
  mem->free_local_double(ijsq_contrib);
  if (bs3_eq_bs4)
    mem->free_local_double(ijrq_contrib);
  delete[] vector1[0]; delete[] vector1;
  delete[] vector3[0]; delete[] vector3;
}

////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
