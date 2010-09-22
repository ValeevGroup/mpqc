//
// transform_12inds.cc
//
// Copyright (C) 2001 Edward Valeev
//
// Author: Edward Valeev <evaleev@vt.edu>
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
#include <math/scmat/blas.h>
#include <chemistry/qc/mbptr12/transform_12inds.h>
#include <chemistry/qc/mbptr12/print.h>

using namespace std;
using namespace sc;

#define PRINT0Q 0
#define PRINT1Q 0
#define ALL_TASKS_ON_SAME_NODE 1
#define PRINT_NUM_TE_TYPES 1

// The FAST_BUT_WRONG flags is useful for exercising the communications
// layer.  It causes the first and second quarter transformation to be
// omitted, but all communication is still performed.  This permits
// problems in communications libraries to be more quickly identified.
#define FAST_BUT_WRONG 0

TwoBodyMOIntsTransform_12Inds::TwoBodyMOIntsTransform_12Inds(
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

TwoBodyMOIntsTransform_12Inds::~TwoBodyMOIntsTransform_12Inds()
{
  timer_ = 0;
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
TwoBodyMOIntsTransform_12Inds::run()
{
  const Ref<MemoryGrp>& mem = tform_->mem();
  const Ref<MessageGrp>& msg = tform_->msg();
  Ref<DistArray4> ints_acc = tform_->ints_acc();
  const int me = msg->me();
  const int nproc = msg->n();
  const Ref<OrbitalSpace>& space1 = tform_->space1();
  const Ref<OrbitalSpace>& space2 = tform_->space2();
  const Ref<OrbitalSpace>& space3 = tform_->space3();
  const Ref<OrbitalSpace>& space4 = tform_->space4();

  const Ref<GaussianBasisSet>& bs1 = space1->basis();
  const Ref<GaussianBasisSet>& bs2 = space2->basis();
  const Ref<GaussianBasisSet>& bs3 = space3->basis();
  const Ref<GaussianBasisSet>& bs4 = space4->basis();
  const bool bs1_eq_bs2 = (bs1 == bs2);
  const bool bs3_eq_bs4 = (bs3 == bs4);

  const bool dynamic = tform_->dynamic();
  const double print_percent = tform_->print_percent();

  const int ni = ni_;
  const int rank1 = space1->rank();
  const int rank2 = space2->rank();
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

  //find the type of integrals which is antisymmetric with respect to permuting functions of particle 1
  int tbtype_anti1 = -1;
  const Ref<TwoBodyOperSetDescr>& descr = tbint_->descr();
  const unsigned int ntypes = descr->size();
  for(unsigned int t=0; t<ntypes; ++t) {
      const TwoBodyOper::type ttype = descr->opertype(t);
      Ref<TwoBodyOperDescr> intdescr = TwoBodyOper::descr(ttype);
      if (intdescr->perm_symm(1) == -1) tbtype_anti1 = t;
  }

  double** vector1 = new double*[nbasis1];
  double** vector2 = new double*[nbasis2];
  vector1[0] = new double[rank1*nbasis1];
  vector2[0] = new double[rank2*nbasis2];
  for(int i=1; i<nbasis1; i++) vector1[i] = vector1[i-1] + rank1;
  for(int i=1; i<nbasis2; i++) vector2[i] = vector2[i-1] + rank2;
  space1->coefs().convert(vector1);
  space2->coefs().convert(vector2);

  /*-------------------------------------------------------------
    Get pointers to integral buffers
   -------------------------------------------------------------*/
  const int num_te_types = tform_->num_te_types();
  const double **intbuf = new const double*[num_te_types];
  for(int te_type=0; te_type<num_te_types; te_type++)
    intbuf[te_type] = tbint_->buffer( descr->opertype(te_type) );

  /*-----------------------------------------------------
    Allocate buffers for partially transformed integrals
   -----------------------------------------------------*/
  double *ijrs_contrib;  // local contributions to integral_ijrs
  double **rsiq_ints = new double*[num_te_types];     // quarter-transformed integrals for each RS pair
  for(int te_type=0;te_type<num_te_types;te_type++) {
    rsiq_ints[te_type] = new double[ni*nbasis2*nfuncmax3*nfuncmax4];
  }
  ijrs_contrib  = mem->malloc_local_double(ni*rank2*nfuncmax3*nfuncmax4);
  const int nij = ni*rank2;
  const int niq = ni*nbasis2;
  double* ij_ints = new double[nij];

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

  Ref<GPetiteList4> p4list
    = GPetiteListFactory::plist4(bs1,bs2,bs3,bs4);

#if FAST_BUT_WRONG
  for(int te_type=0;te_type<num_te_types;te_type++) {
    bzerofast(rsiq_ints[te_type], ni*nbasis2*nfuncmax3*nfuncmax4);
  }
  bzerofast(ijrs_contrib, ni*nbasis2*nfuncmax3*nfuncmax4);
#endif

  int R = 0;
  int S = 0;
  int RS_count = 0;
  while (shellpairs.get_task(S,R)) {
    // if bs3_eq_bs4 then S >= R always (see sc::exp::DistShellPair)
    int nr = bs3->shell(R).nfunction();
    int r_offset = bs3->shell_to_function(R);

    int ns = bs4->shell(S).nfunction();
    int s_offset = bs4->shell_to_function(S);

    const int nrs = nr*ns;

    if (debug_ >= DefaultPrintThresholds::fine && (print_index++)%print_interval == 0) {
      lock_->lock();
      ExEnv::outn() << scprintf("%d:%d: (PQ|%d %d) %d%%",
			       me,mythread_,R,S,(100*print_index)/work_per_thread)
		   << endl;
      lock_->unlock();
    }
    if (debug_ >= DefaultPrintThresholds::fine && (print_index)%time_interval == 0) {
      lock_->lock();
      ExEnv::outn() << scprintf("timer for %d:%d:",me,mythread_) << endl;
      timer_->print();
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
	int deg = p4list->in(P,Q,R,S);
	if (deg == 0)
	  continue;
	double symfac = (double) deg;

        if (tbint_->log2_shell_bound(P,Q,R,S) < tol_) {
          continue;  // skip shell quartets less than tol
	}

        aoint_computed_++;

        timer_->enter("AO integrals");
        tbint_->compute_shell(P,Q,R,S);
        timer_->exit("AO integrals");

#if PRINT0Q
    {
      // each task take its turn to write to the file, in case all tasks live on the same node
      for (int proc = 0; proc < nproc; ++proc) {
        if (me == proc) { // my turn to write

          lock_->lock();
          string filename = tform_->type() + "." + tform_->name()
                            + ".0q.dat";
          ios_base::openmode mode = ios_base::trunc;
          if (RS_count != 0 || ALL_TASKS_ON_SAME_NODE)
            mode = ios_base::app;
          ofstream ints_file(filename.c_str(), mode);

          const int ntetypes = std::min((int) num_te_types,
                                        PRINT_NUM_TE_TYPES);
          for (int te_type = 0; te_type < ntetypes; te_type++) {
            for (int p = 0; p < np; p++) {
              int pp = p + p_offset;
              for (int q = 0; q < nq; q++) {
                int qq = q + q_offset;
                for (int r = 0; r < nr; r++) {
                  int rr = r + r_offset;
                  for (int s = 0; s < ns; s++) {
                    int ss = s + s_offset;
                    double value = intbuf[te_type][s + ns * (r + nr * (q
                        + nq * p))];
                    ints_file << scprintf("0Q: type = %d |(%d %d|%d %d)| = %12.8f\n",
                                          te_type, pp, qq, rr, ss, fabs(value));
                  }
                }
              }
            }
          }

          ints_file.close();
          lock_->unlock();
        }
        msg->sync();
      }
    }
#endif

        timer_->enter("1. q.t.");

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
	timer_->exit("1. q.t.");

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

    timer_->enter("2. q.t.");
    // Begin second quarter transformation;
    // generate (ij|rs) stored as ijrs

    bzerofast(ijrs_contrib, ni*rank2*nrs);

    for(int te_type=0; te_type<num_te_types; te_type++) {
      const double *rsiq_ptr = rsiq_ints[te_type];
      double *ijrs_ptr = ijrs_contrib;

      for (int bf3 = 0; bf3 < nr; bf3++) {
        int smin = (bs3_eq_bs4 && R == S) ? bf3 : 0;
        rsiq_ptr += smin*niq;
        ijrs_ptr += smin;

        for (int bf4 = smin; bf4 <ns; bf4++) {

          // second quarter transform
          // ij = iq * qj
          const char notransp = 'n';
          const double one = 1.0;
          const double zero = 0.0;
          F77_DGEMM(&notransp,&notransp,&rank2,&ni,&nbasis2,&one,vector2[0],&rank2,
                    rsiq_ptr,&nbasis2,&zero,ij_ints,&rank2);

          // Copy rsij integrals into ijrs integrals
          const double* ij_src = ij_ints;
          double* ij_dst = ijrs_ptr;
          for(int i=0; i<ni; i++)
            for(int j=0; j<rank2; j++) {
              *ij_dst = *ij_src++;
              ij_dst += nrs;
            }

          rsiq_ptr += niq;
          ijrs_ptr++;

        }
      }

      // Send the integrals out
      double* ijr_ptr = ijrs_contrib;
      for (int i=0; i<ni; i++) {
        for (int j=0; j<rank2; j++) {

          int ij_proc =  (i*rank2 + j)%nproc;
          int ij_index = (i*rank2 + j)/nproc;
          const size_t ijrs_start = (size_t)(num_te_types*ij_index + te_type) * memgrp_blksize;
          size_t ijr_offset = (size_t)s_offset + r_offset*nbasis4 + ijrs_start;

          for (int bf3 = 0; bf3 < nr; bf3++) {
            // Sum the ijrs_contrib to the appropriate place
            mem->sum_reduction_on_node(ijr_ptr, ijr_offset, ns, ij_proc);
            ijr_offset += nbasis4;
            ijr_ptr += ns;
          }
        }
      }
    }

    timer_->exit("2. q.t.");

    ++RS_count;
  }         // exit while get_task

  if (debug_) {
    lock_->lock();
    ExEnv::outn() << scprintf("%d:%d: done with get_task loop",me,mythread_) << endl;
    lock_->unlock();
  }

  delete[] ij_ints;
  for(int te_type=0; te_type<num_te_types; te_type++) {
    delete[] rsiq_ints[te_type];
  }
  delete[] rsiq_ints;
  mem->free_local_double(ijrs_contrib);
  delete[] vector1[0]; delete[] vector1;
  delete[] vector2[0]; delete[] vector2;
  delete[] intbuf;
}

size_t
TwoBodyMOIntsTransform_12Inds::compute_required_dynamic_memory(const TwoBodyMOIntsTransform& tform,
                                                               int ibatchsize)
{
  const Ref<OrbitalSpace>& space1 = tform.space1();
  const Ref<OrbitalSpace>& space2 = tform.space2();
  const Ref<OrbitalSpace>& space3 = tform.space3();
  const Ref<OrbitalSpace>& space4 = tform.space4();

  const Ref<GaussianBasisSet>& bs1 = space1->basis();
  const Ref<GaussianBasisSet>& bs2 = space2->basis();
  const Ref<GaussianBasisSet>& bs3 = space3->basis();
  const Ref<GaussianBasisSet>& bs4 = space4->basis();
  const int rank1 = space1->rank();
  const int rank2 = space2->rank();
  const int nbasis1 = bs1->nbasis();
  const int nbasis2 = bs2->nbasis();
  const int nfuncmax3 = bs3->max_nfunction_in_shell();
  const int nfuncmax4 = bs4->max_nfunction_in_shell();
  const unsigned int num_te_types = tform.num_te_types();

  const size_t coefs1 = rank1*nbasis1;
  const size_t coefs2 = rank2*nbasis2;
  const size_t iqrs = num_te_types * ibatchsize * nbasis2 * nfuncmax3 * nfuncmax4;
  const size_t ijrs = ibatchsize * rank2 * nfuncmax3 * nfuncmax4;

  return (coefs1 + coefs2 + iqrs + ijrs) * sizeof(double);
}

////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
