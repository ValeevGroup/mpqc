//
// trans12_grt.cc
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
#include <chemistry/qc/mbptr12/trans12_grt.h>
#include <chemistry/qc/basis/distshpair.h>

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

R12A_GRT_12Qtr::R12A_GRT_12Qtr(int mythread_a, int nthread_a,
                                 int me_a, int nproc_a,
                                 const Ref<MemoryGrp> &mem_a,
                                 const Ref<MessageGrp> &msg_a,
                                 const Ref<ThreadLock> &lock_a,
                                 const Ref<GaussianBasisSet> &basis_a,
                                 const Ref<GaussianBasisSet> &aux_basis_a,
                                 const Ref<TwoBodyInt> &tbint_a,
			         int nocc_a, int nocc_act_a,
                                 double **scf_vector_a,
                                 double tol_a, int debug_a,
                                 int dynamic_a, double print_percent_a)
{
  msg = msg_a;
  mythread = mythread_a;
  nthread = nthread_a;
  lock = lock_a;
  basis = basis_a;
  aux_basis = aux_basis_a;
  // aux_basis = basis_a;
  tbint = tbint_a;
  nocc = nocc_a;
  nocc_act = nocc_act_a;
  me = me_a;
  nproc = nproc_a;
  tol = tol_a;
  mem = mem_a;
  scf_vector = scf_vector_a;
  debug = debug_a;
  dynamic_ = dynamic_a;
  print_percent_ = print_percent_a;

  // Use the same basis for now
  bs1_ = basis;
  bs2_ = basis;
  bs3_ = basis;
  bs4_ = basis;

  aoint_computed = 0;
  timer = new RegionTimer();
}

R12A_GRT_12Qtr::~R12A_GRT_12Qtr()
{
}

void
R12A_GRT_12Qtr::run()
{
  bool bs1_eq_bs2 = (bs1_ == bs2_);
  if (!bs1_eq_bs2)
    throw std::runtime_error("R12A_GRT_12Qtr: bs1 != bs2");
  bool bs3_eq_bs4 = (bs3_ == bs4_);
  if (!bs3_eq_bs4)
    throw std::runtime_error("R12A_GRT_12Qtr: bs3 != bs4");

  int te_type;
  int P,Q,R,S;
  int p,q,r,s;
  int np,nq,nr,ns;
  int bf1,bf2,bf3,bf4;
  int p_offset,q_offset,r_offset,s_offset;
  int offset;
  int nfuncmax1 = bs1_->max_nfunction_in_shell();
  int nfuncmax2 = bs2_->max_nfunction_in_shell();
  int nfuncmax3 = bs3_->max_nfunction_in_shell();
  int nfuncmax4 = bs4_->max_nfunction_in_shell();
  int nsh1 = bs1_->nshell();
  int nsh2 = bs2_->nshell();
  int nsh3 = bs3_->nshell();
  int nsh4 = bs4_->nshell();
  int nbasis1 = bs1_->nbasis();
  int nbasis2 = bs2_->nbasis();
  int nbasis3 = bs3_->nbasis();
  int nbasis4 = bs4_->nbasis();
  double dtol = pow(2.0,tol);
  double *iqjs_ptr;
  double *iqrs_ptr, *iprs_ptr;
  double *c_pi, *c_qi;
  double tmpval;
  int i,j;

  /*-------------------------------------------------------------
    Find integrals buffers to 1/r12, r12, and [r12,T1] integrals
   -------------------------------------------------------------*/
  const int num_te_types = 3;
  enum te_types {eri=0, r12=1, r12t1=2};
  const double *intbuf[num_te_types];
  intbuf[eri] = tbint->buffer(TwoBodyInt::eri);
  intbuf[r12] = tbint->buffer(TwoBodyInt::r12);
  intbuf[r12t1] = tbint->buffer(TwoBodyInt::r12t1);

  /*-----------------------------------------------------
    Allocate buffers for partially transformed integrals
   -----------------------------------------------------*/
  double *iqjs_contrib[num_te_types];  // local contributions to integral_iqjs
  double *iqjr_contrib[num_te_types];  // local contributions to integral_iqjr
  double *integral_iqrs[num_te_types]; // quarter transformed two-el integrals
  for(te_type=0;te_type<num_te_types;te_type++) {
    iqjs_contrib[te_type]  = mem->malloc_local_double(nbasis2*nfuncmax4);
    // bs3_eq_bs4
    iqjr_contrib[te_type]  = mem->malloc_local_double(nbasis2*nfuncmax4);
    //  lock->lock();
    integral_iqrs[te_type] = new double[ni*nbasis2*nfuncmax3*nfuncmax4];
    //  lock->unlock();
  }

  /*-----------------------------
    Initialize work distribution
   -----------------------------*/
  DistShellPair shellpairs(msg,nthread,mythread,lock,bs3_,bs4_,dynamic_);
  shellpairs.set_debug(debug);
  if (debug) shellpairs.set_print_percent(print_percent_/10.0);
  else shellpairs.set_print_percent(print_percent_);
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
    lock->lock();
    ExEnv::outn() << scprintf("%d:%d: starting get_task loop",me,mythread) << endl;
    lock->unlock();
  }

  // Assuming all basis sets are the same (bs1_eq_bs2 and bs3_eq_bs4)
  Ref<GenPetite4> p4list = construct_gpetite(bs1_,bs2_,bs3_,bs4_);

#if FAST_BUT_WRONG
  for(te_type=0;te_type<num_te_types;te_type++) {
    bzerofast(integral_iqrs[te_type], ni*nbasis2*nfuncmax3*nfuncmax4);
    bzerofast(iqjs_contrib[te_type], nbasis2*nfuncmax4);
    bzerofast(iqjr_contrib[te_type], nbasis2*nfuncmax4);
    }
#endif

  R = 0;
  S = 0;
  while (shellpairs.get_task(S,R)) {
    nr = bs3_->shell(R).nfunction();
    r_offset = bs3_->shell_to_function(R);
    
    ns = bs4_->shell(S).nfunction();
    s_offset = bs4_->shell_to_function(S);

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
      timer->print();
      lock->unlock();
    }

#if !FAST_BUT_WRONG
    for(te_type=0;te_type<num_te_types;te_type++)
      bzerofast(integral_iqrs[te_type], ni*nbasis2*ns*nr);
//       bzerofast(integral_iqrs[te_type], ni*nbasis2*nfuncmax3*nfuncmax4);

    for (P=0; P<nsh1; P++) {
      np = bs1_->shell(P).nfunction();
      p_offset = bs1_->shell_to_function(P);

      int Qmax = (bs1_eq_bs2 ? P : nsh2-1);
      for (Q=0; Q<=Qmax; Q++) {
	nq = bs2_->shell(Q).nfunction();
	q_offset = bs3_->shell_to_function(Q);

	// check if symmetry unique and compute degeneracy
	int deg = p4list->in_p4(P,Q,R,S);
	double symfac = (double) deg;
	if (deg == 0)
	  continue;

        if (tbint->log2_shell_bound(P,Q,R,S) < tol) {
          continue;  // skip ereps less than tol
	}

        aoint_computed++;

        timer->enter("grt");
        tbint->compute_shell(P,Q,R,S);
        timer->exit("grt");

        timer->enter("1. q.t.");

        // Begin first quarter transformation;
        // generate (iq|rs) for i active
	for(te_type=0; te_type<num_te_types; te_type++) {
	  offset = nr*ns*nbasis2;
	  const double *pqrs_ptr = intbuf[te_type];
	  for (bf1 = 0; bf1 < np; bf1++) {
	    p = p_offset + bf1;
	    for (bf2 = 0; bf2 < nq; bf2++) {
	      q = q_offset + bf2;

	      // bs1_eq_bs2
	      if (p < q) {
		pqrs_ptr = &intbuf[te_type][ns*nr*(bf2+1 + nq*bf1)];
		continue; // skip to next q value
	      }

	      for (bf3 = 0; bf3 < nr; bf3++) {
		r = r_offset + bf3;

		for (bf4 = 0; bf4 < ns; bf4++) {
		  s = s_offset + bf4;

		  // bs3_eq_bs4
		  if (s < r) {
		    pqrs_ptr++;
		    continue; // skip to next bf4 value
                  }

		  if (fabs(*pqrs_ptr) > dtol) {
		    iprs_ptr = &integral_iqrs[te_type][bf4 + ns*(p + nbasis1*bf3)];
		    // nbasis1 == nbasis2
		    iqrs_ptr = &integral_iqrs[te_type][bf4 + ns*(q + nbasis1*bf3)];
		    c_qi = &scf_vector[q][i_offset];
		    c_pi = &scf_vector[p][i_offset];
		    tmpval = *pqrs_ptr;
		    // multiply each integral by its symmetry degeneracy factor
		    tmpval *= symfac;
#if 1 // this code has conditionals in the inner loop (apparently can be faster)
		    for (i=0; i<ni; i++) {
		      // bs1_eq_bs2
		      if (te_type!=2)
			*iprs_ptr += *c_qi++*tmpval;
		      else
			*iprs_ptr -= *c_qi++*tmpval;
		      iprs_ptr += offset;
		      // bs1_eq_bs2
		      if (p != q) {
			*iqrs_ptr += *c_pi++*tmpval;
			iqrs_ptr += offset;
                      }
                    } // exit i loop
#else // this code has conditionals moved out of the inner loop
                    if (te_type!=2) {
                      // bs1_eq_bs2
		      if (p == q) {
                        // te_type != 2, p == q
                        for (i=0; i<ni; i++) {
                          *iprs_ptr += *c_qi++*tmpval;
                          iprs_ptr += offset;
                          }
                        }
                      else {
                        // te_type != 2, p != q
                        for (i=0; i<ni; i++) {
                          *iprs_ptr += *c_qi++*tmpval;
                          iprs_ptr += offset;
                          *iqrs_ptr += *c_pi++*tmpval;
                          iqrs_ptr += offset;
                          }
                        }
                      }
                    else {
                      // bs1_eq_bs2
		      if (p == q) {
                        // te_type == 2, p == q
                        for (i=0; i<ni; i++) {
                          *iprs_ptr -= *c_qi++*tmpval;
                          iprs_ptr += offset;
                          }
                        }
                      else {
                        // te_type == 2, p != q
                        for (i=0; i<ni; i++) {
                          *iprs_ptr -= *c_qi++*tmpval;
                          iprs_ptr += offset;
                          *iqrs_ptr += *c_pi++*tmpval;
                          iqrs_ptr += offset;
                          }
                        }
                      }
#endif
                  }   // endif

		  pqrs_ptr++;
                } // exit bf4 loop
              }   // exit bf3 loop
            }     // exit bf2 loop
          }       // exit bf1 loop
	  // end of first quarter transformation
	}
	timer->exit("1. q.t.");

        }           // exit P loop
      }             // exit Q loop
#endif // !FAST_BUT_WRONG

#if PRINT1Q
      {
      lock->lock();
      for(te_type=0; te_type<PRINT_NUM_TE_TYPES; te_type++) {
	double *tmp = integral_iqrs[te_type];
	for (int i = 0; i<ni; i++) {
	  for (int r = 0; r<nr; r++) {
	    for (int q = 0; q<nbasis2; q++) {
	      for (int s = 0; s<ns; s++) {
		printf("1Q: (%d %d|%d %d) = %12.8f\n",
		       i+i_offset,q,r+r_offset,s+s_offset,*tmp);
		tmp++;
              }
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
      for(te_type=0; te_type<PRINT_BIGGEST_INTS_NUM_TE_TYPES; te_type++) {
	double *tmp = integral_iqrs[te_type];
	for (int i = 0; i<ni; i++) {
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
      lock->unlock();
      }
#endif

    timer->enter("2. q.t.");
    // Begin second quarter transformation;
    // generate (iq|jr) for i active and j active or frozen
    for (i=0; i<ni; i++) {
      for (j=0; j<nocc_act; j++) {
	int j_offset = nocc - nocc_act;
	int ij_proc =  (i*nocc_act + j)%nproc;
	int ij_index = (i*nocc_act + j)/nproc;
	int ijsq_start[num_te_types];
	ijsq_start[0] = num_te_types*nbasis2*nbasis4*ij_index;

	for(te_type=0; te_type<num_te_types; te_type++) {
	  if (te_type)
	    ijsq_start[te_type] = ijsq_start[te_type-1] + nbasis2*nbasis4;

#if !FAST_BUT_WRONG
	  bzerofast(iqjs_contrib[te_type], nbasis2*ns);
	  // bs3_eq_bs4
	  bzerofast(iqjr_contrib[te_type], nbasis2*nr);
// 	  bzerofast(iqjs_contrib[te_type], nbasis2*nfuncmax4);
// 	  // bs3_eq_bs4
// 	  bzerofast(iqjr_contrib[te_type], nbasis2*nfuncmax4);

	  for (bf1=0; bf1<ns; bf1++) {
	    s = s_offset + bf1;
	    double *c_sj = &scf_vector[s][j+j_offset];
	    // bs3_eq_bs4
	    double *iqjr_ptr = iqjr_contrib[te_type];
	    for (bf2=0; bf2<nr; bf2++) {
	      r = r_offset + bf2;
	      // bs3_eq_bs4
	      if (r > s) {
		break; // skip to next bf1 value
              }
	      // bs3_eq_bs4
	      double c_rj = scf_vector[r][j+j_offset];
	      iqjs_ptr = &iqjs_contrib[te_type][bf1*nbasis2];
	      iqrs_ptr = &integral_iqrs[te_type][bf1 + ns*nbasis2*(bf2 + nr*i)];
#if 1 // this code has conditionals in the inner loop (apparently can be faster)
	      for (q=0; q<nbasis2; q++) {
		*iqjs_ptr++ += c_rj * *iqrs_ptr;
		// bs3_eq_bs4
		if (r != s) *iqjr_ptr += *c_sj * *iqrs_ptr;
		iqjr_ptr++;
		iqrs_ptr += ns;
              } // exit q loop
#else // this code has conditionals removed from the inner loop
              // bs3_eq_bs4
              if (r == s) {
                for (q=0; q<nbasis2; q++) {
                  *iqjs_ptr++ += c_rj * *iqrs_ptr;
                  iqrs_ptr += ns;
                  }
                }
              else {
                double c_sj_val = *c_sj;
                for (q=0; q<nbasis2; q++) {
                  double iqrs_val = *iqrs_ptr;
                  *iqjs_ptr++ += c_rj * iqrs_val;
                  *iqjr_ptr += c_sj_val * iqrs_val;
                  iqjr_ptr++;
                  iqrs_ptr += ns;
                  }
                }
#endif
            }   // exit bf2 loop
          }     // exit bf1 loop
#endif // !FAST_BUT_WRONG
	  
	  // We now have contributions to iqjs and iqjr for one pair i,j,
	  // all q, r in R and s in S; send iqjs and iqjr to the node
	  // (ij_proc) which is going to have this ij pair

	  // Sum the iqjs_contrib to the appropriate place
	  int ij_offset = nbasis2*s_offset + ijsq_start[te_type];
	  mem->sum_reduction_on_node(iqjs_contrib[te_type],
				     ij_offset, ns*nbasis2, ij_proc);
	  
	  // bs3_eq_bs4
	  ij_offset = nbasis2*r_offset + ijsq_start[te_type];
	  mem->sum_reduction_on_node(iqjr_contrib[te_type],
				     ij_offset, nr*nbasis2, ij_proc);
	  
	}
      }     // exit j loop
    }       // exit i loop
    // end of second quarter transformation
    timer->exit("2. q.t.");
  }         // exit while get_task

  if (debug) {
    lock->lock();
    ExEnv::outn() << scprintf("%d:%d: done with get_task loop",me,mythread) << endl;
    lock->unlock();
  }

  //  lock->lock();
  for(te_type=0; te_type<num_te_types; te_type++) {
    delete[] integral_iqrs[te_type];
    mem->free_local_double(iqjs_contrib[te_type]);
    mem->free_local_double(iqjr_contrib[te_type]);
  }
  //  lock->unlock();
}

////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
