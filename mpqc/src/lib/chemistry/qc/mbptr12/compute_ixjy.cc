//
// compute_ixjy.cc
//
// Copyright (C) 2004 Edward Valeev
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

#include <stdexcept>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

#include <scconfig.h>
#include <util/misc/formio.h>
#include <util/misc/timer.h>
#include <util/class/class.h>
#include <util/state/state.h>
#include <util/state/state_text.h>
#include <util/state/state_bin.h>
#include <math/scmat/matrix.h>
#include <chemistry/qc/mbptr12/transform_ixjy.h>
#include <chemistry/qc/mbptr12/transform_123inds.h>

using namespace std;
using namespace sc;

#define SINGLE_THREAD_E12   0

#define PRINT2Q 0
#define PRINT3Q 0
#define PRINT4Q 0
#define PRINT4Q_MP2 0
#define PRINT_NUM_TE_TYPES 3
#define PRINT_R12_INTERMED 0

#if PRINT_BIGGEST_INTS
BiggestContribs biggest_ints_1(4,40);
#endif

#define WRITE_DOUBLES 0

#if PRINT_CONTRIB
static void
sw(int&i,int&j)
{
  int tmp = i;
  i = j;
  j = tmp;
}

static void
print_contrib(double tmpval, int num, int onum,
              int P,int Q,int R,int S, int p,int q,int r,int s)
{

  printf("noncanon: z(%d)(%d %d %d %d)(%d %d %d %d) contrib = % 6.4f\n",
         num, P, Q, R, S, p, q, r, s, tmpval);
  printf("noncanon: z(%d)(%d %d %d %d)(%d %d %d %d) contrib = % 6.4f\n",
         onum, P, Q, R, S, p, q, r, s, -tmpval);

  if (p < q) {
      sw(p,q); sw(P,Q);
    }
  if (r < s) {
      sw(r,s); sw(R,S);
    }
  if (p < r || (p == r && q < s)) {
      sw(P,R); sw(p,r);
      sw(Q,S); sw(q,s);
    }

  printf("z(%d)(%d %d %d %d)(%d %d %d %d) contrib = % 6.4f\n",
         num, P, Q, R, S, p, q, r, s, tmpval);
  printf("z(%d)(%d %d %d %d)(%d %d %d %d) contrib = % 6.4f\n",
         onum, P, Q, R, S, p, q, r, s, -tmpval);
}
#endif

static inline void increment_ij(int& i, int& j, int n);


/*-------------------------------------
  Based on MBPT2::compute_mp2_energy()
 -------------------------------------*/
void
TwoBodyMOIntsTransform_ixjy::compute()
{
  init_acc();
  if (ints_acc_->is_committed())
    return;
  
  Ref<Integral> integral = factory_->integral();
  Ref<GaussianBasisSet> bs1 = space1_->basis();
  Ref<GaussianBasisSet> bs2 = space2_->basis();
  Ref<GaussianBasisSet> bs3 = space3_->basis();
  Ref<GaussianBasisSet> bs4 = space4_->basis();
  int rank1 = space1_->rank();
  int rank2 = space2_->rank();
  int rank3 = space3_->rank();
  int rank4 = space4_->rank();
  int nbasis1 = bs1->nbasis();
  int nbasis2 = bs2->nbasis();
  int nbasis3 = bs3->nbasis();
  int nbasis4 = bs4->nbasis();
  int nshell1 = bs1->nshell();
  int nshell2 = bs2->nshell();
  int nshell3 = bs3->nshell();
  int nshell4 = bs4->nshell();
  int nfuncmax1 = bs1->max_nfunction_in_shell();
  int nfuncmax2 = bs2->max_nfunction_in_shell();
  int nfuncmax3 = bs3->max_nfunction_in_shell();
  int nfuncmax4 = bs4->max_nfunction_in_shell();
  enum te_types {eri=0, r12=1, r12t1=2};

  // log2 of the erep tolerance
  // (erep < 2^tol => discard)
  const int tol = (int) (-10.0/log10(2.0));  // discard ints smaller than 10^-20

  int aoint_computed = 0; 

  BiggestContribs biggest_coefs(5,10);

#if PRINT_BIGGEST_INTS
  BiggestContribs biggest_ints_2(4,40);
  BiggestContribs biggest_ints_2s(4,40);
  BiggestContribs biggest_ints_3a(4,40);
  BiggestContribs biggest_ints_3(4,40);
#endif

  tim_enter("r12a-sbs-mem");

  int me = msg_->me();
  int nproc = msg_->n();

  if (debug_ >= 0)
    ExEnv::out0() << endl << indent
                  << "Entered (ix|jy) integrals evaluator" << endl;
  if (debug_ >= 1)
    ExEnv::out0() << indent << scprintf("nproc = %i", nproc) << endl;

  int restart_orbital = ints_acc_.null() ? 0 : ints_acc_->next_orbital();
  if (restart_orbital && debug_ >= 1) {
    ExEnv::out0() << indent
                  << scprintf("Restarting at orbital %d",
                              restart_orbital) << endl;
  }
  
  // Compute the storage remaining for the integral routines
  size_t dyn_mem = distsize_to_size(compute_transform_dynamic_memory_(batchsize_));
  
  if (debug_ >= 1) {
  ExEnv::out0() << indent
		<< "Memory available per node:      " << memory_ << " Bytes"
		<< endl;
  ExEnv::out0() << indent
		<< "Static memory used per node:    " << mem_static_ << " Bytes"
		<< endl;
  ExEnv::out0() << indent
		<< "Total memory used per node:     " << dyn_mem+mem_static_ << " Bytes"
		<< endl;
  ExEnv::out0() << indent
		<< "Memory required for one pass:   "
		<< compute_transform_dynamic_memory_(rank1)+mem_static_
		<< " Bytes"
		<< endl;
  ExEnv::out0() << indent
		<< "Minimum memory required:        "
		<< compute_transform_dynamic_memory_(1)+mem_static_
		<< " Bytes"
		<< endl;
  ExEnv::out0() << indent
		<< "Batch size:                     " << batchsize_
		<< endl;
  if (dynamic_)
    ExEnv::out0() << indent << "Using dynamic load balancing." << endl;

  /*
  ExEnv::out0() << indent
		<< scprintf(" npass  rest  nbasis  nshell  nfuncmax") << endl;
  ExEnv::out0() << indent
		<< scprintf("  %-4i   %-3i   %-5i    %-4i     %-3i",
			    npass,rest,nbasis,nshell,nfuncmax)
		<< endl;
  ExEnv::out0() << indent
		<< scprintf(" nocc   nvir   nfzc   nfzv") << endl;
  ExEnv::out0() << indent
		<< scprintf("  %-4i   %-4i   %-4i   %-4i",
			    nocc,nvir,nfzc,nfzv)
		<< endl;
  */
  }
  
  int nijmax = compute_nij(batchsize_,rank3,nproc,me);
  

  RefSCMatrix coefs1 = space1_->coefs();
  vector<int> mosym1 = space1_->mosym();
  RefSCMatrix coefs2 = space2_->coefs();
  vector<int> mosym2 = space2_->mosym();
  RefSCMatrix coefs3 = space3_->coefs();
  vector<int> mosym3 = space3_->mosym();
  RefSCMatrix coefs4 = space4_->coefs();
  vector<int> mosym4 = space4_->mosym();

  /////////////////////////////////////
  //  Begin transformation loops
  /////////////////////////////////////

  // debug print
  if (debug_ >= 2) {
    ExEnv::outn() << indent
		  << scprintf("node %i, begin loop over i-batches",me) << endl;
  }
  // end of debug print

  // Initialize the integrals
  integral->set_storage(memory_ - dyn_mem);
  Ref<TwoBodyInt>* tbints_ = new Ref<TwoBodyInt>[thr_->nthread()];
  for (int i=0; i<thr_->nthread(); i++) {
    tbints_[i] = integral->grt();
  }
  if (debug_ >= 1)
    ExEnv::out0() << indent << scprintf("Memory used for integral storage:       %i Bytes",
      integral->storage_used()) << endl;

  Ref<ThreadLock> lock = thr_->new_lock();
  TwoBodyMOIntsTransform_123Inds** e12thread = new TwoBodyMOIntsTransform_123Inds*[thr_->nthread()];
  for (int i=0; i<thr_->nthread(); i++) {
    e12thread[i] = new TwoBodyMOIntsTransform_123Inds(
    i, thr_->nthread(), me, nproc, mem_, msg_, lock, bs, bs, tbints_[i],
    nocc, nocc_act, scf_vector, tol, debug_, r12info()->dynamic(), r12info()->print_percent());
  }

  
  /*-----------------------------------

    Start the integrals transformation

   -----------------------------------*/
  tim_enter("mp2-r12/a passes");
  if (me == 0 && mole->if_to_checkpoint() && ints_acc_->can_restart()) {
    StateOutBin stateout(mole->checkpoint_file());
    SavableState::save_state(mole,stateout);
    ExEnv::out0() << indent << "Checkpointed the wave function" << endl;
  }

  for (int pass=0; pass<npass; pass++) {

    ExEnv::out0() << indent << "Beginning pass " << pass+1 << endl;

    int i_offset = restart_orbital + pass*ni + nfzc;
    if ((pass == npass - 1) && (rest != 0)) ni = rest;

    // Compute number of of i,j pairs on each node during current pass for
    // two-el integrals
    int nij = 0;
    {
      int index = 0;
      for (int i=0; i<ni; i++) {
	for (int j=0; j<nocc_act; j++) {
	  if (index++ % nproc == me) nij++;
	}
      }
    }

    // debug print
    if (debug_)
      ExEnv::outn() << indent << "node " << me << ", nij = " << nij << endl;
    // end of debug print

    // Allocate and initialize some arrays
    // (done here to avoid having these arrays
    // overlap with arrays allocated later)

    // Allocate (and initialize) some arrays

    double* integral_ijsq = (double*) mem->localdata();
    bzerofast(integral_ijsq, (num_te_types*nij*nbasis*nbasis));
    integral_ijsq = 0;
    mem->sync();
    ExEnv::out0() << indent
		  << scprintf("Begin loop over shells (grt, 1.+2. q.t.)") << endl;

    // Do the two electron integrals and the first two quarter transformations
    tim_enter("grt+1.qt+2.qt");
    for (int i=0; i<thr->nthread(); i++) {
      e12thread[i]->set_i_offset(i_offset);
      e12thread[i]->set_ni(ni);
      thr->add_thread(i,e12thread[i]);
#     if SINGLE_THREAD_E12
      e12thread[i]->run();
#     endif
    }
#   if !SINGLE_THREAD_E12
    thr->start_threads();
    thr->wait_threads();
#   endif
    tim_exit("grt+1.qt+2.qt");
    ExEnv::out0() << indent << "End of loop over shells" << endl;

    mem->sync();  // Make sure ijsq is complete on each node before continuing
    integral_ijsq = (double*) mem->localdata();

#if PRINT2Q
    if (me == 0) {
      int index = 0;
      int ij_index = 0;
      int j_offset = nfzc;
      for (int i = 0; i<ni; i++) {
	for (int j = 0; j<nocc_act; j++) {
	  if (index++ % nproc == me) {
	    double *integral_ij_offset = integral_ijsq + num_te_types*nbasis*nbasis*ij_index;
	    for(int te_type=0; te_type<PRINT_NUM_TE_TYPES; te_type++,integral_ij_offset+=nbasis*nbasis) {
	      for (int s = 0; s<nbasis; s++) {
		double *integral_ijsq_ptr = integral_ij_offset + s*nbasis;
		for (int q = 0; q<nbasis; q++) {
		  printf("2Q: (%d %d|%d %d) = %12.8f\n",
		       i+i_offset,q,j+j_offset,s,*integral_ijsq_ptr);
		integral_ijsq_ptr++;
		}
	      }
	    }
	    ij_index++;
	  }
	}
      }
    }
#endif

    // Allocate and initialize some arrays
    double* ijsx_tmp = new double[nbasis];

    ExEnv::out0() << indent << "Begin third q.t." << endl;
    tim_enter("3. q.t.");
    // Begin third quarter transformation;
    // generate (ix|js) for i act, j act, and x any MO
    int index = 0;
    int ij_index = 0;
    for (int i=0; i<ni; i++) {
      for (int j=0; j<nocc_act; j++) {
        if (index++ % nproc == me) {
	  double* integral_ij_offset = integral_ijsq + num_te_types*nbasis*nbasis*ij_index;
	  for(int te_type=0; te_type<num_te_types; te_type++,integral_ij_offset+=nbasis*nbasis) {
	    for (int s=0; s<nbasis; s++) {
	      double* integral_ijsq_ptr = integral_ij_offset + s*nbasis;
	      bzerofast(ijsx_tmp, nbasis);
	      for (int q=0; q<nbasis; q++,integral_ijsq_ptr++) {
		double* ijsx_ptr = ijsx_tmp;
		double* c_qx = scf_vector[q];
		double tmpval = *integral_ijsq_ptr;
#if PRINT_BIGGEST_INTS
		biggest_ints_2.insert(tmpval,i+i_offset,j,s,q);
		if ((i+i_offset==104 && j == 1)
		    ||(i+i_offset==104 && j == 2)) {
		  biggest_ints_2s.insert(tmpval,i+i_offset,j,s,q);
		}
#endif
		for (int x=0; x<noso; x++) {
		  *ijsx_ptr++ += *c_qx++ * tmpval;
		}
              }   // exit q loop

	      // Put ixjs into integral_ijsq, while overwriting what was there;
	      // i.e., integral_ijsq will now contain three-quarter transformed
	      // integrals ijsx
	      integral_ijsq_ptr = integral_ij_offset + s*nbasis;
	      double* ijsx_ptr = ijsx_tmp;
	      for (int x=0; x<noso; x++) {
#if PRINT_BIGGEST_INTS
		if (x>=nocc) {
		  biggest_ints_3a.insert(*ijsx_ptr,i+i_offset,j,s,x-nocc);
                }
#endif
		*integral_ijsq_ptr++ = *ijsx_ptr++;
              }
            }   // exit s loop
	  }
	  ij_index++;
	}     // endif
      }       // exit j loop
    }         // exit i loop
    // end of third quarter transformation
    tim_exit("3. q.t.");
    ExEnv::out0() << indent << "End of third q.t." << endl;


    delete[] ijsx_tmp;

    // The array of half-transformed integrals integral_ijsq has now
    // been overwritten by three-quarter transformed integrals ixjs;
    // rename the array integral_ixjs, where x = any MO
    double* integral_ijsx = integral_ijsq;

#if PRINT3Q
    if (me == 0) {
      int index = 0;
      int ij_index = 0;
      int j_offset = nfzc;
      for (int i = 0; i<ni; i++) {
	for (int j = 0; j<nocc_act; j++) {
	  if (index++ % nproc == me) {
	    double *integral_ij_offset = integral_ijsx + num_te_types*nbasis*nbasis*ij_index;
	    for(int te_type=0; te_type<PRINT_NUM_TE_TYPES; te_type++,integral_ij_offset+=nbasis*nbasis) {
	      for (int s = 0; s<nbasis; s++) {
		double *integral_ijsx_ptr = integral_ij_offset + s*nbasis;
		for (int x = 0; x<noso; x++) {
		  printf("3Q: (%d %d|%d %d) = %12.8f\n",
		       i+i_offset,x,j+j_offset,s,*integral_ijsx_ptr);
		integral_ijsx_ptr++;
		}
	      }
	    }
	    ij_index++;
	  }
	}
      }
    }
#endif

    double* ijyx_tmp = new double[noso];
    // in ijyx: i act; x any MO, j act; y any MO.

    // Begin fourth quarter transformation
    ExEnv::out0() << indent << "Begin fourth q.t." << endl;
    tim_enter("4. q.t.");
    index = 0;
    ij_index = 0;
    for (int i=0; i<ni; i++) {
      for (int j=0; j<nocc_act; j++) {
        if (index++ % nproc == me) {
	  double *integral_ij_offset = integral_ijsx + num_te_types*nbasis*nbasis*ij_index;
	  for(int te_type=0; te_type<num_te_types; te_type++,integral_ij_offset+=nbasis*nbasis) {
	    for (int x=0; x<noso; x++) {
	      bzerofast(ijyx_tmp, noso);
	      double* ijsx_ptr = integral_ij_offset + x;
	      for (int s=0; s<nbasis; s++) {
		double* c_sy = scf_vector[s];
		double* ijyx_ptr = ijyx_tmp;
		double tmpval = *ijsx_ptr;
#if PRINT_BIGGEST_INTS
		biggest_ints_3.insert(tmpval,i+i_offset,j,s,x);
		if ((i+i_offset==105 && j == 2 && s == 170 && x == 3)
		    ||(i+i_offset==102 && j == 2 && s == 170 && x == 2)) {
		  ExEnv::outn() << scprintf("3/4: %3d %3d %3d %3d: %16.10f",
					   i+i_offset, j, s, x-nocc)
			       << endl;
                }
#endif
		for (int y=0; y<noso; y++) {
		  *ijyx_ptr++ += *c_sy++ * tmpval;
		} // exit y loop
		ijsx_ptr += nbasis;
	      }   // exit s loop
	      // Put ijyx_tmp into ijsx for one i,x,j while
	      // overwriting elements of ijsx
	      ijsx_ptr = integral_ij_offset + x;
	      double* ijyx_ptr = ijyx_tmp;
	      for (int y=0; y<noso; y++) {
		*ijsx_ptr = *ijyx_ptr++;
		ijsx_ptr += nbasis;
	      } // exit y loop
	    }   // exit x loop
	  }
          ij_index++;
	}   // endif
      }     // exit j loop
    }       // exit i loop
    // end of fourth quarter transformation
    tim_exit("4. q.t.");
    ExEnv::out0() << indent << "End of fourth q.t." << endl;
    
    // The array integral_ijsx has now been overwritten by MO integrals ijyx
    // rename the array mo_int
    double* mo_int = integral_ijsx;
    delete[] ijyx_tmp;

    // Zero out nonsymmetric integrals
    {
    int index = 0;
    int ij_index = 0;
    int j_offset = nfzc;
    for (int i = 0; i<ni; i++) {
      for (int j = 0; j<nocc_act; j++) {
	if (index++ % nproc == me) {
	  double *integral_ij_offset = mo_int + num_te_types*nbasis*nbasis*ij_index;
	  for(int te_type=0; te_type<PRINT_NUM_TE_TYPES; te_type++,integral_ij_offset+=nbasis*nbasis) {
	    for (int y = 0; y < noso; y++) {
	      double *integral_ijyx_ptr = integral_ij_offset + y*nbasis;
	      for (int x = 0; x<noso; x++) {
		if (( mo_irrep[i+i_offset] ^
		      mo_irrep[j+j_offset] ^
		      mo_irrep[x] ^
		      mo_irrep[y]) ) {
		  *integral_ijyx_ptr = 0.0;
		}
		integral_ijyx_ptr++;
	      }
	    }
	  }
	  ij_index++;
	}
      }
    }
    }

#if PRINT4Q
    if (me == 0) {
      int index = 0;
      int ij_index = 0;
      int j_offset = nfzc;
      for (int i = 0; i<ni; i++) {
	for (int j = 0; j<nocc_act; j++) {
	  if (index++ % nproc == me) {
	    double *integral_ij_offset = mo_int + num_te_types*nbasis*nbasis*ij_index;
	    for(int te_type=0; te_type<PRINT_NUM_TE_TYPES; te_type++,integral_ij_offset+=nbasis*nbasis) {
	      for (int y = 0; y < noso; y++) {
		double *integral_ijyx_ptr = integral_ij_offset + y*nbasis;
		for (int x = 0; x<noso; x++) {
		  printf("4Q: type = %d (%d %d|%d %d) = %12.8f\n",
		       te_type,i+i_offset,x,j+j_offset,y,*integral_ijyx_ptr);
		integral_ijyx_ptr++;
		}
	      }
	    }
	    ij_index++;
	  }
	}
      }
    }
#endif

    // For now compute MP2 energy to verify the transformed ERIs
    tim_enter("compute emp2");

    index = 0;
    ij_index = 0;
    for (int i=0; i<ni; i++) {
      int ii = i + i_offset - nfzc;
      for (int j=0; j<nocc_act; j++) {
	int jj = j;
        double ecorr_ij = 0.0;
	// alpha-alpha pair energy is 0 when i == j - thus we only need strictly i > j
	int ij_aa = (ii > jj) ? ii*(ii-1)/2 + jj : jj*(jj-1)/2 + ii;
	int ij_ab = ii*nocc_act + jj;
	double eaa = 0.0;
	double eab = 0.0;

        if (index++ % nproc == me) {
	  for (int b=0; b<nvir; b++) {
	    double* iajb_ptr = &mo_int[nocc + nbasis*(b+nocc + nbasis*num_te_types*ij_index)];
	    double* ibja_ptr = &mo_int[b+nocc + nbasis*(nocc + nbasis*num_te_types*ij_index)];
	    for (int a=0; a<nvir; a++) {
#if PRINT4Q_MP2
	      printf("4Q: (%d %d|%d %d) = %12.8f\n",
		     i,a+nocc,j+nfzc,b+nocc,*iajb_ptr);
#endif
	      double delta_ijab = evals[i_offset+i]+evals[j+nfzc]-evals[nocc+a]-evals[nocc+b];
	      // only include determinants with unique coefficients
	      if (a>=b && i_offset+i>=j+nfzc) {
		if (a>b && i_offset+i>j+nfzc) {
		  // aaaa or bbbb
		  biggest_coefs.insert(*iajb_ptr - *ibja_ptr,
				       i_offset+i,j,a,b,1111);
		  // aabb or bbaa or abba or baab
		  biggest_coefs.insert(*ibja_ptr,i_offset+i,j,b,a,1212);
		} // endif
		// aabb or bbaa or abba or baab
		biggest_coefs.insert(*iajb_ptr,i_offset+i,j,a,b,1212);
	      } // endif
	      double tmpval;
	      if (ii != jj) {
		tmpval = (*iajb_ptr - *ibja_ptr)*(*iajb_ptr - *ibja_ptr)/delta_ijab;
		eaa += tmpval;
	      }
	      tmpval = 0.5*(*iajb_ptr * *iajb_ptr + *ibja_ptr * *ibja_ptr)/delta_ijab;
	      eab += tmpval;
	      ecorr_ij += *iajb_ptr*(2.0 * *iajb_ptr - *ibja_ptr)/delta_ijab;
	      iajb_ptr++;
	      ibja_ptr += nbasis;
	    } // exit a loop
	  }   // exit b loop
	  ij_index++;
	  if (ii > jj)
	    emp2_aa.set_element(ij_aa,eaa);
	  emp2_ab.set_element(ij_ab,eab);
	}     // endif
	if (debug_) {
	  msg_->sum(ecorr_ij);
	  ExEnv::out0() << indent
			<< scprintf("correlation energy for pair %3d %3d = %16.12f",
				    i+i_offset, j+nfzc, ecorr_ij)
			<< endl;
	}
      }         // exit j loop
    }           // exit i loop
    tim_exit("compute emp2");
    
    // debug print
    if (debug_) {
      ExEnv::out0() << indent << "End of ecorr" << endl;
    }
    // end of debug print
    
    integral_ijsq = 0;
    mem->sync(); // Make sure MO integrals are complete on all nodes before continuing

    // Push locally stored integrals to an accumulator
    // This could involve storing the data to disk or simply remembering the pointer
    tim_enter("MO ints store");
    ints_acc_->store_memorygrp(mem,ni);
    tim_exit("MO ints store");
    mem->sync();

    if (me == 0 && mole->if_to_checkpoint() && ints_acc_->can_restart()) {
      StateOutBin stateout(mole->checkpoint_file());
      SavableState::save_state(mole,stateout);
      ExEnv::out0() << indent << "Checkpointed the wave function" << endl;
    }

  } // end of loop over passes
  tim_exit("mp2-r12/a passes");
  if (debug_)
    ExEnv::out0() << indent << "End of mp2-r12/a transformation" << endl;
  // Done storing integrals - commit the content
  // WARNING: it is not safe to use mem until deactivate has been called on the accumulator
  //          After that deactivate the size of mem will be 0 [mem->set_localsize(0)]
  ints_acc_->commit();


  
  for (int i=0; i<thr->nthread(); i++) {
    delete e12thread[i];
  }
  delete[] e12thread;
  
  
  delete[] tbints_; tbints_ = 0;
  delete[] scf_vector;
  delete[] scf_vector_dat;
  delete[] evals;
  tim_exit("r12a-sbs-mem");

  evaluated_ = true;

  if (me == 0 && mole->if_to_checkpoint()) {
    StateOutBin stateout(mole->checkpoint_file());
    SavableState::save_state(mole,stateout);
    ExEnv::out0() << indent << "Checkpointed the wave function" << endl;
  }
  
  return;
}


// This function increments a pair of indices i and j subject to a condition i<n, j<n, i>=j
// If i==n and j==n then it sets i=0 and j=0
static inline void increment_ij(int& i, int& j, int n)
{
  if (i == j) {
    i = (++i)%n;
    j = 0;
  }
  else
    j++;
}

////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
