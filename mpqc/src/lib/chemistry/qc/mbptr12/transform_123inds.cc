//
// transform_123inds.cc
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
#include <util/misc/timer.h>
#include <chemistry/qc/basis/gpetite.h>
#include <chemistry/qc/mbpt/bzerofast.h>
#include <chemistry/qc/mbpt/util.h>
#include <chemistry/qc/basis/distshpair.h>
#include <chemistry/qc/mbptr12/transform_123inds.h>

using namespace std;
using namespace sc;

extern BiggestContribs biggest_ints_1;

#define PRINT1Q 0
#define PRINT_NUM_TE_TYPES 1
#define PRINT_BIGGEST_INTS_NUM_TE_TYPES 1

// The FAST_BUT_WRONG flags is useful for exercising the communications
// layer.  It causes the first and second quarter transformation to be
// omitted, but all communication is still performed.  This permits
// problems in communications libraries to be more quickly identified.
#define FAST_BUT_WRONG 0

TwoBodyMOIntsTransform_123Inds::TwoBodyMOIntsTransform_123Inds(
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

TwoBodyMOIntsTransform_123Inds::~TwoBodyMOIntsTransform_123Inds()
{
}

/*
  Distribute work by SR
  
   for all PQ
    compute unique (PQ|RS)
    transform to (RS|IM) where M are all AOs for basis set 2
   end PQ

   use BLAS to transform each rsIM to rsIX
   transform RSIX to IJXS and accumulate to the tasks that holds respective ij-pairs.

  end SR
*/

void
TwoBodyMOIntsTransform_123Inds::run()
{
  Ref<MemoryGrp> mem = tform_->mem();
  Ref<MessageGrp> msg = tform_->msg();
  int me = msg->me();
  int nproc = msg->n();
  Ref<MOIndexSpace> space1 = tform_->space1();
  Ref<MOIndexSpace> space2 = tform_->space2();
  Ref<MOIndexSpace> space3 = tform_->space3();
  Ref<MOIndexSpace> space4 = tform_->space4();

  Ref<GaussianBasisSet> bs1 = space1->basis();
  Ref<GaussianBasisSet> bs2 = space2->basis();
  Ref<GaussianBasisSet> bs3 = space3->basis();
  Ref<GaussianBasisSet> bs4 = space4->basis();
  bool bs1_eq_bs2 = (bs1 == bs2);
  if (!bs1_eq_bs2)
    throw std::runtime_error("TwoBodyMOIntsTransform_ixjy_12Qtr::run() -- bs1 != bs2");
  bool bs3_eq_bs4 = (bs3 == bs4);
  if (!bs3_eq_bs4)
    throw std::runtime_error("TwoBodyMOIntsTransform_ixjy_12Qtr::run() -- bs3 != bs4");

  bool dynamic = tform_->dynamic();
  double print_percent = tform_->print_percent();

  int nfuncmax1 = bs1->max_nfunction_in_shell();
  int nfuncmax2 = bs2->max_nfunction_in_shell();
  int nfuncmax3 = bs3->max_nfunction_in_shell();
  int nfuncmax4 = bs4->max_nfunction_in_shell();
  int nsh1 = bs1->nshell();
  int nsh2 = bs2->nshell();
  int nsh3 = bs3->nshell();
  int nsh4 = bs4->nshell();
  int nbasis1 = bs1->nbasis();
  int nbasis2 = bs2->nbasis();
  int nbasis3 = bs3->nbasis();
  int nbasis4 = bs4->nbasis();
  double dtol = pow(2.0,tol_);

  /*-------------------------------------------------------------
    Find integrals buffers to 1/r12, r12, and [r12,T1] integrals
   -------------------------------------------------------------*/
  int num_te_types = tform_->num_te_types();
  enum te_types {eri=0, r12=1, r12t1=2};
  const double *intbuf[num_te_types];
  intbuf[eri] = tbint->buffer(TwoBodyInt::eri);
  intbuf[r12] = tbint->buffer(TwoBodyInt::r12);
  intbuf[r12t1] = tbint->buffer(TwoBodyInt::r12t1);

  /*-----------------------------------------------------
    Allocate buffers for partially transformed integrals
   -----------------------------------------------------*/
  double *ijxs_contrib[num_te_types];  // local contributions to integral_ijxs
  double *ijxr_contrib[num_te_types];  // local contributions to integral_ijxr
  double *rsiq_ints[num_te_types];     // quarter-transformed integrals for each RS pair
  double *rsix_ints[num_te_types];     // 2 quarter-transformed integrals for each RS pair
  for(int te_type=0;te_type<num_te_types;te_type++) {
    ijxs_contrib[te_type]  = mem->malloc_local_double(nbasis2*nfuncmax4);
    if (bs3_eq_bs4)
      ijxr_contrib[te_type]  = mem->malloc_local_double(nbasis2*nfuncmax4);
    else
      ijxr_contrib[te_type]  = NULL;
    rsiq_ints[te_type] = new double[ni_*nbasis2*nfuncmax3*nfuncmax4];
    rsix_ints[te_type] = new double[ni_*nrank2*nfuncmax3*nfuncmax4];
  }

  /*-----------------------------
    Initialize work distribution
   -----------------------------*/
  sc::exp::DistShellPair shellpairs(msg,nthread,mythread,lock_,bs3,bs4,dynamic_);
  shellpairs.set_debug(debug);
  if (debug_) shellpairs.set_print_percent(print_percent/10.0);
  else shellpairs.set_print_percent(print_percent);
  int work_per_thread = bs3_eq_bs4 ? 
    ((nsh3*(nsh3+1))/2)/(nproc*nthread) :
    (nsh3*nsh4)/(nproc*nthread) ;
  int print_interval = work_per_thread/100;
  int time_interval = work_per_thread/10;
  int print_index = 0;
  if (print_interval == 0) print_interval = 1;
  if (time_interval == 0) time_interval = 1;
  if (work_per_thread == 0) work_per_thread = 1;

  if (debug) {
    lock_->lock();
    ExEnv::outn() << scprintf("%d:%d: starting get_task loop",me,mythread) << endl;
    lock_->unlock();
  }

  // Assuming all basis sets are the same (bs1_eq_bs2 and bs3_eq_bs4)
  // I don't know yet how to overcome the static type problem here
  canonical_aaaa c4(bs1,bs2,bs3,bs4);
  Ref<GPetite4<canonical_aaaa> > p4list
    = new GPetite4<canonical_aaaa>(bs1,bs2,bs3,bs4,c4);

#if FAST_BUT_WRONG
  for(int te_type=0;te_type<num_te_types;te_type++) {
    bzerofast(rsiq_ints[te_type], ni_*nbasis2*nfuncmax3*nfuncmax4);
    bzerofast(ijxs_contrib[te_type], nbasis2*nfuncmax4);
    if (bs3_eq_bs4)
      bzerofast(ijxr_contrib[te_type], nbasis2*nfuncmax4);
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

    if (debug > 1 && (print_index++)%print_interval == 0) {
      lock_->lock();
      ExEnv::outn() << scprintf("%d:%d: (PQ|%d %d) %d%%",
			       me,mythread,R,S,(100*print_index)/work_per_thread)
		   << endl;
      lock_->unlock();
    }
    if (debug > 1 && (print_index)%time_interval == 0) {
      lock_->lock();
      ExEnv::outn() << scprintf("timer for %d:%d:",me,mythread) << endl;
      timer_->print();
      lock_->unlock();
    }

#if !FAST_BUT_WRONG
    // Zero out 1 q.t. storage
    for(int te_type=0;te_type<num_te_types;te_type++)
      bzerofast(rsiq_ints[te_type], nrs*ni_*nbasis2);

    for (int P=0; P<nsh1; P++) {
      int np = bs1->shell(P).nfunction();
      int p_offset = bs1->shell_to_function(P);

      int Qmax = (bs1_eq_bs2 ? P : nsh2-1);
      for (int Q=0; Q<=Qmax; Q++) {
	int nq = bs2->shell(Q).nfunction();
	int q_offset = bs3->shell_to_function(Q);
        
	// check if symmetry unique and compute degeneracy
	int deg = p4list->in_p4(P,Q,R,S);
	if (deg == 0)
	  continue;
	double symfac = (double) deg;

        if (tbint->log2_shell_bound(P,Q,R,S) < tol_) {
          continue;  // skip shell quartets less than tol
	}

        aoint_computed++;

        timer_->enter("AO integrals");
        tbint->compute_shell(P,Q,R,S);
        timer_->exit("AO integrals");

        timer_->enter("1. q.t.");

        // Begin first quarter transformation;
        // generate (iq|rs) for i active
        // if bs1_eq_bs2 then (ip|rs) are also generated
        // store the integrals as rsiq
	for(int te_type=0; te_type<num_te_types; te_type++) {
	  const double *pqrs_ptr = intbuf[te_type];

	  for (int bf1 = 0; bf1 < np; bf1++) {
	    int p = p_offset + bf1;
            int qmax = (bs1_eq_bs2 && P == Q) ? p : nq-1;

	    for (int bf2 = 0; bf2 <= qmax; bf2++) {
	      int q = q_offset + bf2;

	      for (int bf3 = 0; bf3 < nr; bf3++) {
                int smin = (bs3_eq_bs4 && R == S) ? 0 : nr;
                pqrs_ptr += smin;

		for (int bf4 = smin; bf4 <nq; bf4++) {

                  // Only transform integrals larger than the threshold
		  if (fabs(*pqrs_ptr) > dtol) {

		    const double* rsiq_ptr = &rsiq_ints[te_type][bf2 + ni_*(bs4 + nq*bf3)];
		    const double* c_pi = vector1[p] + i_offset_;

                    const double* rsip_ptr;
		    const double* c_qi;
                    if (bs1_eq_bs2) {
		      rsip_ptr = &rsiq_ints[te_type][bf1 + ni_*(bs4 + nq*bf3)];
		      c_qi = vector1[q] + i_offset_;
                    }
                    
		    double rsiq_int_contrib = *pqrs_ptr;
		    // multiply each integral by its symmetry degeneracy factor
		    rsiq_int_contrib *= symfac;
                    
                    if (bs1_eq_bs2) {

                      double rsip_int_contrib = rsiq_int_contrib;
                      if (te_type == r12t1)
                      rsip_int_contrib = -1.0*rsiq_int_contrib;

                      if (p == q) {
                        for (i=0; i<ni_; i++) {
                          *rsiq_ptr += *c_pi++ * rsiq_int_contrib;
                          rsiq_ptr += nbasis2;
                        }
                      }
                      else {
                        // p != q
                        for (i=0; i<ni_; i++) {
                          *rsip_ptr += *c_qi++ * rsip_int_contrib;
                          rsip_ptr += nbasis2;
                          *rsiq_ptr += *c_pi++ * rsiq_int_contrib;
                          rsiq_ptr += nbasis2;
                        }
                      }
                        
                    }
                    else {

                      for (i=0; i<ni_; i++) {
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
	timer_->exit("1. q.t.");

        }           // exit P loop
      }             // exit Q loop
#endif // !FAST_BUT_WRONG

#if PRINT1Q
      {
      lock_->lock();
      for(int te_type=0; te_type<PRINT_NUM_TE_TYPES; te_type++) {
	double *tmp = rsiq_ints[te_type];
        for (int r = 0; r<nr; r++) {
          for (int s = 0; s<ns; s++) {
            for (int i = 0; i<ni_; i++) {
              for (int q = 0; q<nbasis2; q++) {
		printf("1Q: (%d %d|%d %d) = %12.8f\n",
		       i+i_offset_, q, bf3+r_offset, bf4+s_offset, *tmp);
		tmp++;
              }
            }
          }
        }
      }
      lock_->unlock();
      }
#endif
#if PRINT_BIGGEST_INTS
      {
      lock_->lock();
      for(te_type=0; te_type<PRINT_BIGGEST_INTS_NUM_TE_TYPES; te_type++) {
	double *tmp = integral_iqrs[te_type];
	for (int i = 0; i<ni_; i++) {
	  for (int r = 0; r<nr; r++) {
	    for (int q = 0; q<nbasis2; q++) {
	      for (int s = 0; s<ns; s++) {
		biggest_ints_1.insert(*tmp,i+i_offset,q,r+r_offset,s+s_offset);
		tmp++;
              }
            }
          }
        }
      }
      lock_->unlock();
      }
#endif

    timer_->enter("2. q.t.");
    // Begin second quarter transformation;
    // generate (ix|rs) stored as rsix

    for(int te_type=0; te_type<num_te_types; te_type++) {
      const double *rsiq_ptr = rsiq_ints[te_type];
      const double *rsix_ptr = rsix_ints[te_type];

      for (int bf3 = 0; bf3 < nr; bf3++) {
        int smin = (bs3_eq_bs4 && R == S) ? 0 : nr;
        rsiq_ptr += smin*ni_*nbasis2;
        rsix_ptr += smin*ni_*nrank2;

        for (int bf4 = smin; bf4 <nq; bf4++) {

          // second quarter transform
          // ix = iq * qx
          const char notransp = 'n';
          const double one = 1.0;
          const double zero = 0.0;
          F77_DGEMM(&notransp,&notransp,&rank2,&ni_,&nbasis2,&one,vector2[0],&nbasis2,
                    rsiq_ptr,&nbasis2,&zero,rsix_ptr,&rank2);

          rsiq_ptr += ni_*nbasis2;
          rsix_ptr += ni_*nrank2;

        }
      }
    }
    timer_->exit("2. q.t.");

    
    timer_->enter("3. q.t.");
    // Begin third quarter transformation;
    // generate (ix|js) stored as ijxs (also generate (ix|jr, if needed)

    for(int te_type=0; te_type<num_te_types; te_type++) {
      const double *rsix_ptr = rsix_ints[te_type];

      for (int i=0; i<ni_; i++) {
        for (int j=0; j<rank3; j++) {

#if !FAST_BUT_WRONG
          bzerofast(ijxs_contrib[te_type], rank2*ns);
          if (bs3_eq_bs4)
            bzerofast(ijxr_contrib[te_type], rank2*nr);

          int ij_proc =  (i*rank3 + j)%nproc;
          int ij_index = (i*rank3 + j)/nproc;
          const size_t ijxq_start = (size_t)(num_te_types*ij_index + te_type) * ints_acc->blocksize();
            
          if (bs3_eq_bs4) {

            for (int bf3 = 0; bf3 < nr; bf3++) {
              int smin = (bs3_eq_bs4 && R == S) ? 0 : nr;
              rsix_ptr += smin*ni_*nrank2;

              for (int bf4 = smin; bf4 <nq; bf4++) {

                // third quarter transform
                // rs = js
                // rs = jr
                
                const double* ijxs_ptr = ijxs_contrib[te_type] + bf4;
                const double* ijxr_ptr = ijxr_contrib[te_type] + bf3;
                const double* i_ptr = rsix_ptr + i*nrank2;

                const double c_rj = vector3[r][j];
                const double c_sj = vector3[s][j];

                if (r != s) {
                  for (x=0; x<rank2; x++) {

                    double value = *i_ptr++;
                    *ijxs_ptr += c_rj * value;
                    ijxs += ns;
                    *ijxr_ptr += c_sj * value;
                    ijxr += nr;

                  }
                }
                else {
                  for (x=0; x<rank2; x++) {

                    double value = *i_ptr++;
                    *ijxs_ptr += c_rj * value;
                    ijxs += ns;

                  }
                }
              }
            }
          }
          else {

            for (int bf3 = 0; bf3 < nr; bf3++) {
              for (int bf4 = 0; bf4 <nq; bf4++) {

                // third quarter transform
                // rs = js
                const double* ijxs_ptr = ijxs_contrib[te_type] + bf4;
                const double* i_ptr = rsix_ptr + i*nrank2;

                const double c_rj = vector3[r][j];

                for (x=0; x<rank2; x++) {

                  double value = *i_ptr++;
                  *ijxs_ptr += c_rj * value;
                  ijxs += ns;
                }
              }
            }
          }

          // We now have contributions to ijxs (and ijxr) for one pair i,j,
          // all x, and s in S (r in R); send ijxs (and ijxr) to the node
          // (ij_proc) which is going to have this ij pair
#endif // !FAST_BUT_WRONG

          // Sum the ijxs_contrib to the appropriate place
          size_t ij_offset = (size_t)rank2*s_offset + ijxq_start[te_type];
          mem->sum_reduction_on_node(ijxs_contrib[te_type],
                                     ij_offset, ns*rank2, ij_proc);

          if (bs3_eq_bs4) {
            size_t ij_offset = (size_t)rank2*r_offset + ijxq_start;
            mem->sum_reduction_on_node(ijxr_contrib[te_type],
                                       ij_offset, nr*rank2, ij_proc);
          }

        }
      }
    }
    timer_->exit("3. q.t.");
          
	  
  }         // exit while get_task

  if (debug) {
    lock_->lock();
    ExEnv::outn() << scprintf("%d:%d: done with get_task loop",me,mythread) << endl;
    lock_->unlock();
  }

  //  lock_->lock();
  for(int te_type=0; te_type<num_te_types; te_type++) {
    delete[] rsiq_ints[te_type];
    delete[] rsix_ints[te_type];
    mem->free_local_double(ijxs_contrib[te_type]);
    if (bs3_eq_bs4)
      mem->free_local_double(ijxr_contrib[te_type]);
  }
  //  lock_->unlock();
}

////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
