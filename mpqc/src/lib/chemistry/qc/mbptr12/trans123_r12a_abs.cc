//
// trans123_r12a_abs.cc
//
// Copyright (C) 2003 Edward Valeev
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
#include <chemistry/qc/mbptr12/trans123_r12a_abs.h>
#include <chemistry/qc/basis/distshpair.h>

using namespace std;
using namespace sc;

extern BiggestContribs biggest_ints_1;

#define USE_SYMMETRY 1
#define PRINT1Q 0
#define PRINT2Q 0
#define PRINT_NUM_TE_TYPES 1
#define PRINT_BIGGEST_INTS_NUM_TE_TYPES 1

R12A_ABS_123Qtr::R12A_ABS_123Qtr(int mythread_a, int nthread_a,
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

  bs1_ = basis;
  bs2_ = basis;
  bs3_ = basis;
  bs4_ = aux_basis;

  aoint_computed = 0;
  timer = new RegionTimer();
}

R12A_ABS_123Qtr::~R12A_ABS_123Qtr()
{
}

/*------------------------------------------------------------
  First 3 steps of the transformation are performed here
  pqrs -> iqrs
  iqrs -> ikrs
  ikrs -> ikjs (stored as ijks to send to the ij-destination)
  
  i - act
  k - occ
  j - act
 ------------------------------------------------------------*/

void
R12A_ABS_123Qtr::run()
{
  bool bs1_eq_bs2 = (bs1_ == bs2_);
  if (!bs1_eq_bs2)
    throw std::runtime_error("R12A_ABS_123Qtr: bs1 != bs2");
  bool bs3_eq_bs4 = (bs3_ == bs4_);

  int P,Q,R,S;
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

  /*-----------------------------------------------------------------------
    Find integrals buffers to 1/r12, r12, [r12,T1], and [r12,T2] integrals
   -----------------------------------------------------------------------*/
  const int num_te_types = 4;
  enum te_types {eri=0, r12=1, r12t1=2, r12t2=3};
  const double *intbuf[num_te_types];
  intbuf[eri] = tbint->buffer(TwoBodyInt::eri);
  intbuf[r12] = tbint->buffer(TwoBodyInt::r12);
  intbuf[r12t1] = tbint->buffer(TwoBodyInt::r12t1);
  intbuf[r12t2] = tbint->buffer(TwoBodyInt::r12t2);

  /*-----------------------------------------------------
    Allocate buffers for partially transformed integrals
   -----------------------------------------------------*/
  double *integral_iqrs[num_te_types]; // quarter transformed two-el integrals
  double *integral_ikrs[num_te_types]; // half-transformed two-el integrals (i act, k occ)
  double *integral_sk;                 // 3 q.t. two-el integrals (ikjs), stored as sk  for one ij at a time
  for(int te_type=0;te_type<num_te_types;te_type++) {
    integral_iqrs[te_type] = new double[ni*nbasis2*nfuncmax3*nfuncmax4];
    integral_ikrs[te_type] = new double[ni*nocc*nfuncmax3*nfuncmax4];
  }
  integral_sk = mem->malloc_local_double(nocc*nfuncmax4);

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

  // Assuming bs1_eq_bs2
  Ref<GenPetite4> p4list = construct_gpetite(bs1_,bs2_,bs3_,bs4_);

  R = 0;
  S = 0;
  while (shellpairs.get_task(R,S)) {
    int nr = bs3_->shell(R).nfunction();
    int r_offset = bs3_->shell_to_function(R);
    int ns = bs4_->shell(S).nfunction();
    int s_offset = bs4_->shell_to_function(S);
    int nrs = nr*ns;

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

    for(int te_type=0;te_type<num_te_types;te_type++) {
      bzerofast(integral_iqrs[te_type], ni*nbasis2*nfuncmax3*nfuncmax4);
      bzerofast(integral_ikrs[te_type], ni*nocc*nfuncmax3*nfuncmax4);
    }

    for (P=0; P<nsh1; P++) {
      int np = bs1_->shell(P).nfunction();
      int p_offset = bs1_->shell_to_function(P);

      int Qmax = (bs1_eq_bs2 ? P : nsh2-1);
      for (Q=0; Q<=Qmax; Q++) {
	int nq = bs2_->shell(Q).nfunction();
	int q_offset = bs2_->shell_to_function(Q);

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
	for(int te_type=0; te_type<num_te_types; te_type++) {
	  int qrs_offset = nrs*nbasis2;
	  const double *pqrs_ptr = intbuf[te_type];
	  for (int bf1 = 0; bf1 < np; bf1++) {
	    int p = p_offset + bf1;
	    for (int bf2 = 0; bf2 < nq; bf2++) {
	      int q = q_offset + bf2;

	      // bs1_eq_bs2
	      if (p < q) {
		pqrs_ptr = &intbuf[te_type][ns*nr*(bf2+1 + nq*bf1)];
		continue; // skip to next q value
	      }

	      int rs_offset = 0;
	      for (int bf3 = 0; bf3 < nr; bf3++) {
		int r = r_offset + bf3;

		for (int bf4 = 0; bf4 < ns; bf4++, rs_offset++) {
		  int s = s_offset + bf4;

		  if (fabs(*pqrs_ptr) > dtol) {
		    double *iprs_ptr = &integral_iqrs[te_type][rs_offset + nrs*p];
		    // bs1_eq_bs2
		    double *iqrs_ptr = &integral_iqrs[te_type][rs_offset + nrs*q];
		    double *c_qi = &scf_vector[q][i_offset];
		    double *c_pi = &scf_vector[p][i_offset];
		    double tmpval = *pqrs_ptr;
		    // multiply each integral by its symmetry degeneracy factor
		    tmpval *= symfac;
		    for (int i=0; i<ni; i++) {
		      // bs1_eq_bs2
		      if (te_type!=2)
			*iprs_ptr += *c_qi++ * tmpval;
		      else
			*iprs_ptr -= *c_qi++ * tmpval;
		      iprs_ptr += qrs_offset;
		      // bs1_eq_bs2
		      if (p != q) {
			*iqrs_ptr += *c_pi++ * tmpval;
			iqrs_ptr += qrs_offset;
                      }
                    } // exit i loop
                  }   // endif

		  pqrs_ptr++;
                } // exit bf4 loop
              }   // exit bf3 loop
            }     // exit bf2 loop
          }       // exit bf1 loop
	  // end of first quarter transformation
	}
	timer->exit("1. q.t.");

      }           // exit Q loop
    }             // exit P loop

#if PRINT1Q
      {
      lock->lock();
      for(int te_type=0; te_type<PRINT_NUM_TE_TYPES; te_type++) {
	double *tmp = integral_iqrs[te_type];
	for (int i = 0; i<ni; i++) {
	  for (int q = 0; q<nbasis2; q++) {
	    for (int r = 0; r<nr; r++) {
	      for (int s = 0; s<ns; s++) {
		printf("1Q: type = %d (%d %d|%d %d) = %12.8f\n",
		       te_type,i+i_offset,q,r+r_offset,s+s_offset,*tmp);
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
      for(int te_type=0; te_type<PRINT_BIGGEST_INTS_NUM_TE_TYPES; te_type++) {
	double *tmp = integral_iqrs[te_type];
	for (int i = 0; i<ni; i++) {
	  for (int q = 0; q<nbasis2; q++) {
	    for (int r = 0; r<nr; r++) {
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
    // generate (ik|rs) for i active and k occupied
    for (int te_type=0; te_type<num_te_types; ++te_type) {
      double *iqrs_ptr = integral_iqrs[te_type];
      for (int i=0; i<ni; i++) {
	for (int q=0; q<nbasis2; q++) {
	  
	  int rs_offset = 0;
	  for (int bf3=0; bf3<nr; bf3++) {
	    int r = r_offset + bf3;
	    for (int bf4=0; bf4<ns; bf4++, ++rs_offset) {
	      int s = s_offset + bf4;

	      double tmpval = *iqrs_ptr++;
	      double *ikrs_ptr = integral_ikrs[te_type] + rs_offset + nrs*nocc*i;
	      double *c_q = scf_vector[q];
	      for (int k=0; k<nocc; k++) {
		*ikrs_ptr += *c_q * tmpval;
		ikrs_ptr += nrs;
		++c_q;
	      } // exit k loop
	    } // exit s loop
	  }   // exit r loop
	}     // exit q loop
      } // exit i loop
    }
    // end of second quarter transformation
    timer->exit("2. q.t.");

#if PRINT2Q
      {
      lock->lock();
      for(int te_type=0; te_type<PRINT_NUM_TE_TYPES; te_type++) {
	double *tmp = integral_ikrs[te_type];
	for (int i = 0; i<ni; i++) {
	  for (int k = 0; k<nocc; k++) {
	    for (int r = 0; r<nr; r++) {
	      for (int s = 0; s<ns; s++) {
		printf("2Q: type = %d (%d %d|%d %d) = %12.8f\n",
		       te_type,i+i_offset,k,r+r_offset,s+s_offset,*tmp);
		tmp++;
              }
            }
          }
        }
      }
      lock->unlock();
      }
#endif

    timer->enter("3. q.t.");
    // Begin second quarter transformation;
    // generate (ik|js) for i active, k occupied, and j active
    for (int i=0; i<ni; i++) {
      for (int j=0; j<nocc_act; j++) {
	int j_offset = nocc - nocc_act;
	int ij_proc =  (i*nocc_act + j)%nproc;
	int ij_index = (i*nocc_act + j)/nproc;
	int ijsk_start[num_te_types];

	ijsk_start[0] = num_te_types*nocc*nbasis4*ij_index;
	for(int te_type=0; te_type<num_te_types; te_type++) {
	  if (te_type)
	    ijsk_start[te_type] = ijsk_start[te_type-1] + nocc*nbasis4;
	  
	  bzerofast(integral_sk, nfuncmax4*nocc);

	  double *ikrs_ptr = integral_ikrs[te_type] + i*nocc*nrs;
	  for (int k=0; k<nocc; k++) {
	    for (int bf3=0; bf3<nr; bf3++) {
	      int r = r_offset + bf3;
	      double c_rj = scf_vector[r][j+j_offset];
	      double *sk_ptr = integral_sk + k;
	      for (int bf4=0; bf4<ns; ++bf4) {
		*sk_ptr += c_rj * *ikrs_ptr;
		++ikrs_ptr;
		sk_ptr += nocc;

	      } // end of s loop
	    } // end of r loop
	  } // end of k loop


	  // We now have contributions to ikjs for from r in R and s in S
	  // send ikjs as ijsk to the node (ij_proc) which is going to have this ij pair
	  
	  int ijs_offset = s_offset*nocc + ijsk_start[te_type];
	  mem->sum_reduction_on_node(integral_sk, ijs_offset, ns*nocc, ij_proc);
	} // end of te_type loop
      } // end of j loop
    } // end of i loop

    timer->exit("3. q.t.");

  }         // exit while get_task

  if (debug) {
    lock->lock();
    ExEnv::outn() << scprintf("%d:%d: done with get_task loop",me,mythread) << endl;
    lock->unlock();
  }

  //  lock->lock();
  for(int te_type=0; te_type<num_te_types; te_type++) {
    delete[] integral_iqrs[te_type];
    delete[] integral_ikrs[te_type];
  }
  mem->free_local_double(integral_sk);
  //  lock->unlock();
}

////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
