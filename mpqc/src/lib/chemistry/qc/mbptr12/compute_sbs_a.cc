//
// compute_sbs_a.cc
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
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/obint.h>
#include <chemistry/qc/basis/symmint.h>
#include <chemistry/qc/mbpt/util.h>
#include <chemistry/qc/mbpt/bzerofast.h>
#include <chemistry/qc/mbptr12/mbptr12.h>
#include <chemistry/qc/mbptr12/trans12_grt.h>
#include <chemistry/qc/mbptr12/r12ia.h>
#include <chemistry/qc/mbptr12/r12ia_memgrp.h>
#include <chemistry/qc/mbptr12/r12ia_node0file.h>
#ifdef HAVE_MPIIO
  #include <chemistry/qc/mbptr12/r12ia_mpiiofile.h>
#endif
#include <chemistry/qc/mbptr12/vxb_eval_info.h>
#include <chemistry/qc/mbptr12/vxb_eval_sbs_a.h>

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
R12IntEval_sbs_A::compute(RefSCMatrix& Vaa, RefSCMatrix& Xaa, RefSCMatrix& Baa,
			  RefSCMatrix& Vab, RefSCMatrix& Xab, RefSCMatrix& Bab,
			  RefSCVector& emp2_aa, RefSCVector& emp2_ab)
{
  if (evaluated_)
    return;
  
  int debug_ = r12info()->debug_level();

  Wavefunction* wfn = r12info()->wfn();
  Ref<Integral> integral = r12info()->integral();
  Ref<GaussianBasisSet> bs = r12info()->basis();
  bool two_basis_form = (bs != r12info()->basis_ri());
  LinearR12::ABSMethod abs_method = r12info()->abs_method();
  Ref<MessageGrp> msg = r12info()->msg();
  Ref<MemoryGrp> mem = r12info()->mem();
  Ref<ThreadGrp> thr = r12info()->thr();
  const int num_te_types = 3;
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

  int me = msg->me();
  int nproc = msg->n();
  const size_t mem_alloc = r12info()->memory();

  int nbasis = bs->nbasis();
  int nfuncmax = bs->max_nfunction_in_shell();
  int nshell = bs->nshell();
  int nocc = r12info()->nocc();
  int nocc_act = r12info()->nocc_act();
  int nfzc = r12info()->nfzc();
  int nfzv = r12info()->nfzv();
  int noso = r12info()->noso();
  int nvir  = noso - nocc;

  double pfac_xy_1, pfac_xy_2;
  if (two_basis_form &&
      ( abs_method == LinearR12::ABS_ABS ||
        abs_method == LinearR12::ABS_ABSPlus ) ) {
    pfac_xy_1 = 0.5;
    pfac_xy_2 = -0.5;
  }
  else {
    pfac_xy_1 = 0.5;
    pfac_xy_2 = 0.5;
  }

  ExEnv::out0() << endl << indent
	       << "Entered OBS A intermediates evaluator" << endl;
  ExEnv::out0() << indent << scprintf("nproc = %i", nproc) << endl;


  // Do a few preliminary tests to make sure the desired calculation
  // can be done (and appears to be meaningful!)

  if (nocc_act <= 0)
    throw std::runtime_error("There are no active occupied orbitals; program exiting");

  if (restart_orbital_) {
    ExEnv::out0() << indent
		  << scprintf("Restarting at orbital %d",
			      restart_orbital_)
		  << endl;
  }

  ////////////////////////////////////////////////////////
  // Compute batch size ni for mp2 loops;
  //
  // The following arrays are kept throughout (all of type double):
  //   scf_vector
  // The following objects are kept throughout:
  //   integrals evaluators
  // memory allocated for these arrays and objects is
  // called mem_static
  //
  ////////////////////////////////////////////////////////
  size_t mem_static = 0;
  int ni = 0;
  if (me == 0) {
    mem_static = nbasis*noso; // scf vector
    mem_static *= sizeof(double);
    int nthreads = thr->nthread();
    mem_static += nthreads * integral->storage_required_grt(bs); // integral evaluators
    ni = compute_transform_batchsize_(mem_alloc, mem_static, nocc_act-restart_orbital_, num_te_types); 
  }

  int max_norb = nocc_act - restart_orbital_;
  if (ni > max_norb)
    ni = max_norb;

  // Send value of ni and mem_static to other nodes
  msg->bcast(ni);
  double mem_static_double = (double) mem_static;
  msg->bcast(mem_static_double);
  mem_static = (size_t) mem_static_double;
  
  // Compute the storage remaining for the integral routines
  size_t dyn_mem = distsize_to_size(compute_transform_dynamic_memory_(ni,nocc_act,num_te_types));
  size_t mem_remaining = 0;
  if (mem_alloc <= (dyn_mem + mem_static)) mem_remaining += 0;
  else mem_remaining += mem_alloc - dyn_mem - mem_static;
  mem_remaining += thr->nthread() * integral->storage_required_grt(bs);
  
  ExEnv::out0() << indent
		<< "Memory available per node:      " << mem_alloc << " Bytes"
		<< endl;
  ExEnv::out0() << indent
		<< "Static memory used per node:    " << mem_static << " Bytes"
		<< endl;
  ExEnv::out0() << indent
		<< "Total memory used per node:     " << dyn_mem+mem_static << " Bytes"
		<< endl;
  ExEnv::out0() << indent
		<< "Memory required for one pass:   "
		<< compute_transform_dynamic_memory_(nocc_act,nocc_act,num_te_types)+mem_static
		<< " Bytes"
		<< endl;
  ExEnv::out0() << indent
		<< "Minimum memory required:        "
		<< compute_transform_dynamic_memory_(1,nocc_act,num_te_types)+mem_static
		<< " Bytes"
		<< endl;
  ExEnv::out0() << indent
		<< "Batch size:                     " << ni
		<< endl;
  
  if (ni == 0)
    throw std::runtime_error("R12IntEval_sbs_A: batch size is 0: more memory or processors are needed");
  
  if (r12info()->dynamic()) {
    ExEnv::out0() << indent << "Using dynamic load balancing." << endl;
  }

  int npass = 0;
  int rest = 0;
  if (ni == nocc_act-restart_orbital_) {
    npass = 1;
    rest = 0;
  }
  else {
    rest = (nocc_act-restart_orbital_)%ni;
    npass = (nocc_act-restart_orbital_ - rest)/ni + 1;
    if (rest == 0) npass--;
  }

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

  int nijmax = 0;
  {
    int index = 0;
    for (int i=0; i<ni; i++) {
      for (int j=0; j<nocc_act; j++) {
	if (index++ % nproc == me) nijmax++;
      }
    }
  }
  

  ////////////////////////////////////////////////
  // The scf vector is distributed between nodes;
  // put a copy of the scf vector on each node;
  ////////////////////////////////////////////////

  RefSCMatrix Scf_Vec = r12info()->scf_vec();
  RefDiagSCMatrix evalmat = r12info()->evals();
  int *mo_irrep = r12info()->orbsym();
  if (debug_ > 1) {
    evalmat.print("eigenvalues");
    Scf_Vec.print("eigenvectors");
  }

  double *scf_vector_dat = new double[nbasis*noso];
  Scf_Vec.t()->convert(scf_vector_dat);

  double* evals = new double[noso];
  double** scf_vector = new double*[nbasis];
  for (int i=0; i<nbasis; i++) {
    scf_vector[i] = &scf_vector_dat[i*noso];
  }
  for (int i=0; i<noso; i++) {
    evals[i] = evalmat(i);
  }
  Scf_Vec = 0;
  evalmat = 0;





  /////////////////////////////////////
  //  Begin MP2 loops
  /////////////////////////////////////

  // debug print
  if (debug_) {
    ExEnv::outn() << indent
		  << scprintf("node %i, begin loop over i-batches",me) << endl;
  }
  // end of debug print

  // Initialize the integrals
  integral->set_storage(mem_remaining);
  Ref<TwoBodyInt>* tbints_ = new Ref<TwoBodyInt>[thr->nthread()];
  for (int i=0; i<thr->nthread(); i++) {
    tbints_[i] = integral->grt();
  }
  ExEnv::out0() << indent
		<< scprintf("Memory used for integral storage:       %i Bytes",
			    integral->storage_used())
		<< endl;


  if (mem.null())
    throw std::runtime_error("R12IntEval_sbs_A: memory group not initialized");
  mem->set_localsize(num_te_types*nijmax*nbasis*nbasis*sizeof(double));
  ExEnv::out0() << indent
	       << "Size of global distributed array:       "
	       << mem->totalsize()
	       << " Bytes"
	       << endl;


  Ref<ThreadLock> lock = thr->new_lock();
  R12A_GRT_12Qtr** e12thread = new R12A_GRT_12Qtr*[thr->nthread()];
  for (int i=0; i<thr->nthread(); i++) {
    e12thread[i] = new R12A_GRT_12Qtr(i, thr->nthread(), me, nproc,
				      mem, msg, lock, bs, bs, tbints_[i],
				      nocc, nocc_act, scf_vector, tol, debug_,
				      r12info()->dynamic(),
                                      r12info()->print_percent());
  }

  
  ///////////////////////////////////////////////////////////
  // Figure out which integrals accumulator should be used
  ///////////////////////////////////////////////////////////

  Ref<R12IntsAcc> r12intsacc;
  R12IntEvalInfo::StoreMethod ints_method = r12info()->ints_method();
  char *r12ints_file = r12info()->ints_file();
  bool restart = (restart_orbital_ > 0);

  switch (ints_method) {

  case R12IntEvalInfo::mem_only:
    if (restart)
      throw std::runtime_error("R12IntEval_sbs_A::compute -- cannot use MemoryGrp-based accumulator when restarting");
    ExEnv::out0() << indent << "Will hold transformed integrals in memory" << endl;
    r12intsacc = new R12IntsAcc_MemoryGrp(mem,num_te_types,nbasis,nbasis,nocc,nfzc);
    break;

  case R12IntEvalInfo::mem_posix:
    if (npass == 1) {
      ExEnv::out0() << indent << "Will hold transformed integrals in memory" << endl;
      r12intsacc = new R12IntsAcc_MemoryGrp(mem,num_te_types,nbasis,nbasis,nocc,nfzc);
      break;
    }
    // else use the next case
      
  case R12IntEvalInfo::posix:
    ExEnv::out0() << indent << "Will use POSIX I/O on node 0 to handle transformed integrals" << endl;
    r12intsacc = new R12IntsAcc_Node0File(mem,r12ints_file,num_te_types,nbasis,nbasis,nocc,nfzc,restart);
    break;

#if HAVE_MPIIO
  case R12IntEvalInfo::mem_mpi:
    if (npass == 1) {
      ExEnv::out0() << indent << "Will hold transformed integrals in memory" << endl;
      r12intsacc = new R12IntsAcc_MemoryGrp(mem,num_te_types,nbasis,nbasis,nocc,nfzc);
      break;
    }
    // else use the next case

  case R12IntEvalInfo::mpi:
    ExEnv::out0() << indent << "Will use MPI-IO (individual I/O) to handle transformed integrals" << endl;
    r12intsacc = new R12IntsAcc_MPIIOFile_Ind(mem,r12ints_file,num_te_types,nbasis,nbasis,nocc,nfzc,restart);
    break;
#endif
  
  default:
    throw std::runtime_error("R12IntEval_sbs_A::compute -- invalid integrals store method");
  }
  free(r12ints_file);


  /////////////////////////////////////////////////
  // zero out data arrays prior to the FIRST pass
  /////////////////////////////////////////////////

  if (!restart) {
    emp2_aa.assign(0.0);
    emp2_ab.assign(0.0);
  }
  
  /*-----------------------------------

    Start the integrals transformation

   -----------------------------------*/
  tim_enter("mp2-r12/a passes");
  if (me == 0 && wfn->if_to_checkpoint() && r12intsacc->can_restart()) {
    StateOutBin stateout(wfn->checkpoint_file());
    SavableState::save_state(wfn,stateout);
    ExEnv::out0() << indent << "Checkpointed the wave function" << endl;
  }

  for (int pass=0; pass<npass; pass++) {

    ExEnv::out0() << indent << "Beginning pass " << pass+1 << endl;

    int i_offset = restart_orbital_ + pass*ni + nfzc;
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
	      if (i != j) {
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
	  msg->sum(ecorr_ij);
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
    r12intsacc->store_memorygrp(mem,ni);
    tim_exit("MO ints store");
    mem->sync();

    if (me == 0 && wfn->if_to_checkpoint() && r12intsacc->can_restart()) {
      current_orbital_ += ni;
      StateOutBin stateout(wfn->checkpoint_file());
      SavableState::save_state(wfn,stateout);
      ExEnv::out0() << indent << "Checkpointed the wave function" << endl;
    }

  } // end of loop over passes
  tim_exit("mp2-r12/a passes");
  if (debug_)
    ExEnv::out0() << indent << "End of mp2-r12/a transformation" << endl;
  // Done storing integrals - commit the content
  // WARNING: it is not safe to use mem until deactivate has been called on the accumulator
  //          After that deactivate the size of mem will be 0 [mem->set_localsize(0)]
  r12intsacc->commit();


  /*-----------------------------------------------
    Compute dipole and quadrupole moment integrals
   -----------------------------------------------*/
  RefSymmSCMatrix MX, MY, MZ, MXX, MYY, MZZ;
  r12info()->compute_multipole_ints(MX,MY,MZ,MXX,MYY,MZZ);
  if (debug_)
    ExEnv::out0() << indent << "Computed multipole moment integrals" << endl;
  
  /*--------------------------------
    Compute MP2-R12/A intermediates
    and collect on node0
   --------------------------------*/
  ExEnv::out0() << indent << "Begin computation of intermediates" << endl;
  tim_enter("mp2-r12a intermeds");
  int naa = (nocc_act*(nocc_act-1))/2;          // Number of alpha-alpha pairs (i > j)
  int nab = nocc_act*nocc_act;                  // Number of alpha-beta pairs
  if (debug_) {
    ExEnv::out0() << indent << "naa = " << naa << endl;
    ExEnv::out0() << indent << "nab = " << nab << endl;
  }
  double *Vaa_ijkl = new double[naa*naa];
  double *Taa_ijkl = new double[naa*naa];
  double *Xaa_ijkl = new double[naa*naa];
  double *Vab_ijkl = new double[nab*nab];
  double *Tab_ijkl = new double[nab*nab];
  double *Xab_ijkl = new double[nab*nab];
  if (debug_)
    ExEnv::out0() << indent << "Allocated intermediates V, X, and T" << endl;
  bzerofast(Vaa_ijkl,naa*naa);
  bzerofast(Taa_ijkl,naa*naa);
  bzerofast(Xaa_ijkl,naa*naa);
  bzerofast(Vab_ijkl,nab*nab);
  bzerofast(Tab_ijkl,nab*nab);
  bzerofast(Xab_ijkl,nab*nab);

  // Compute intermediates
  if (debug_)
    ExEnv::out0() << indent << "Ready to compute intermediates V, X, and T" << endl;
  const int pair_block_size = num_te_types*nbasis*nbasis;
  const double oosqrt2 = 1.0/sqrt(2.0);
  // Compute the number of tasks that have full access to the integrals
  // and split the work among them
  int nproc_with_ints = 0;
  for(int proc=0;proc<nproc;proc++)
    if (r12intsacc->has_access(proc)) nproc_with_ints++;
  int *proc_with_ints = new int[nproc];
  int count = 0;
  for(int proc=0;proc<nproc;proc++)
    if (r12intsacc->has_access(proc)) {
      proc_with_ints[proc] = count;
      count++;
    }
    else
      proc_with_ints[proc] = -1;
  ExEnv::out0() << indent << "Computing intermediates on " << nproc_with_ints
		<< " processors" << endl;

  
  //////////////////////////////////////////////////////////////
  //
  // Evaluation of the intermediates proceeds as follows:
  //
  //    loop over batches of kl, k >= l,  0<=k,l<nocc_act
  //      load (kl|xy), (kl| [T1,r12] |xy), and (lk| [T1,r12] |xy)
  //           (aka kl-sets) into memory
  //
  //      loop over batches of ij, i>=j, 0<=i,j<nocc_act
  //        load (ij|r12|xy) into memory
  //           (aka ij-sets) into memory
  //        compute V[ij][kl] and T[ij][kl] for all ij and kl in
  //                the "direct product" batch
  //      end ij loop
  //    end kl loop
  //
  /////////////////////////////////////////////////////////////////////////////////

  if (r12intsacc->has_access(me)) {
    int nij = nocc_act*(nocc_act+1)/2;
    int kl = 0;
    for(int k=0;k<nocc_act;k++)
      for(int l=0;l<=k;l++,kl++) {
	// Figure out if this task will handle this kl
        int kl_proc = kl%nproc_with_ints;
        if (kl_proc != proc_with_ints[me])
          continue;
	int kl_aa = k*(k-1)/2 + l;
	int kl_ab = k*nocc_act + l;
	int lk_ab = l*nocc_act + k;

        if (debug_)
          ExEnv::outn() << indent << "task " << me << ": working on (k,l) = " << k << "," << l << " " << endl;

	// Get (|1/r12|), (|r12|), and (|[r12,T1]|) integrals
        tim_enter("MO ints retrieve");
        double *klyx_buf_eri = r12intsacc->retrieve_pair_block(k,l,R12IntsAcc::eri);
        double *klyx_buf_r12 = r12intsacc->retrieve_pair_block(k,l,R12IntsAcc::r12);
        double *klyx_buf_r12t1 = r12intsacc->retrieve_pair_block(k,l,R12IntsAcc::r12t1);
	double *lkyx_buf_r12t1 = r12intsacc->retrieve_pair_block(l,k,R12IntsAcc::r12t1);
        tim_exit("MO ints retrieve");

	if (debug_)
          ExEnv::outn() << indent << "task " << me << ": obtained kl blocks" << endl;

	// to avoid every task hitting same ij at the same time, stagger ij-accesses, i.e. each kl task will start with ij=kl+1
	int i = k;
	int j = l;
	increment_ij(i,j,nocc_act);
	for(int ij_counter=0; ij_counter<nij; ij_counter++, increment_ij(i,j,nocc_act)) {

	    int ij_aa = i*(i-1)/2 + j;
	    int ij_ab = i*nocc_act + j;
	    int ji_ab = j*nocc_act + i;

            if (debug_)
              ExEnv::outn() << indent << "task " << me << ": (k,l) = " << k << "," << l << ": (i,j) = " << i << "," << j << endl;
            
            tim_enter("MO ints retrieve");
            double *ijyx_buf_r12 = r12intsacc->retrieve_pair_block(i,j,R12IntsAcc::r12);
            tim_exit("MO ints retrieve");

            if (debug_)
              ExEnv::outn() << indent << "task " << me << ": obtained ij blocks" << endl;

	    double *Vaa_ij = Vaa_ijkl + ij_aa*naa;
            double *Vab_ij = Vab_ijkl + ij_ab*nab;
            double *Vab_ji = Vab_ijkl + ji_ab*nab;
            double *Taa_ij = Taa_ijkl + ij_aa*naa;
            double *Tab_ij = Tab_ijkl + ij_ab*nab;
            double *Tab_ji = Tab_ijkl + ji_ab*nab;
            double *Xaa_ij = Xaa_ijkl + ij_aa*naa;
            double *Xab_ij = Xab_ijkl + ij_ab*nab;
            double *Xab_ji = Xab_ijkl + ji_ab*nab;
            

            tim_enter("MO ints contraction");
            double Vaa_ijkl, Vab_ijkl, Vab_jikl, Vab_ijlk, Vab_jilk;
            double Xaa_ijkl, Xab_ijkl, Xab_jikl, Xab_ijlk, Xab_jilk;
	    double Taa_ijkl, Tab_ijkl, Tab_jikl, Tab_ijlk, Tab_jilk;
	    Vaa_ijkl = Vab_ijkl = Vab_jikl = Vab_ijlk = Vab_jilk = 0.0;
	    Xaa_ijkl = Xab_ijkl = Xab_jikl = Xab_ijlk = Xab_jilk = 0.0;
	    Taa_ijkl = Tab_ijkl = Tab_jikl = Tab_ijlk = Tab_jilk = 0.0;

	    /*----------------------------------
	      Compute (r12)^2 contribution to X
	     ----------------------------------*/
	    double r1r1_ik = -1.0*(MXX->get_element(i,k) + MYY->get_element(i,k) + MZZ->get_element(i,k));
	    double r1r1_il = -1.0*(MXX->get_element(i,l) + MYY->get_element(i,l) + MZZ->get_element(i,l));
	    double r1r1_jk = -1.0*(MXX->get_element(j,k) + MYY->get_element(j,k) + MZZ->get_element(j,k));
	    double r1r1_jl = -1.0*(MXX->get_element(j,l) + MYY->get_element(j,l) + MZZ->get_element(j,l));
	    double r1r2_ijkl = MX->get_element(i,k)*MX->get_element(j,l) +
	      MY->get_element(i,k)*MY->get_element(j,l) +
	      MZ->get_element(i,k)*MZ->get_element(j,l);
	    double r1r2_ijlk = MX->get_element(i,l)*MX->get_element(j,k) +
	      MY->get_element(i,l)*MY->get_element(j,k) +
	      MZ->get_element(i,l)*MZ->get_element(j,k);
	    double delta_ik = (i==k ? 1.0 : 0.0);
	    double delta_il = (i==l ? 1.0 : 0.0);
	    double delta_jk = (j==k ? 1.0 : 0.0);
	    double delta_jl = (j==l ? 1.0 : 0.0);
	    Xab_ijkl += r1r1_ik * delta_jl + r1r1_jl * delta_ik - 2.0*r1r2_ijkl;
	    if (i != j)
	      Xab_jikl += r1r1_jk * delta_il + r1r1_il * delta_jk - 2.0*r1r2_ijlk;
	    if (k != l)
	      Xab_ijlk += r1r1_il * delta_jk + r1r1_jk * delta_il - 2.0*r1r2_ijlk;
	    if (i != j && k != l) {
	      Xaa_ijkl += r1r1_ik * delta_jl + r1r1_jl * delta_ik - 2.0*r1r2_ijkl -
		r1r1_jk * delta_il - r1r1_il * delta_jk + 2.0*r1r2_ijlk;
	      Xab_jilk += r1r1_ik * delta_jl + r1r1_jl * delta_ik - 2.0*r1r2_ijkl;
	    }

            for(int y=0;y<noso;y++) {
	      double pfac_xy;
	      if (y >= nocc)
		pfac_xy = pfac_xy_1;
	      else
		pfac_xy = pfac_xy_2;
              for(int x=0;x<noso;x++) {
                int yx_offset = y*nbasis+x;
                int xy_offset = x*nbasis+y;
                double ij_r12_xy = ijyx_buf_r12[yx_offset];
                double ij_r12_yx = ijyx_buf_r12[xy_offset];
                double kl_eri_xy = klyx_buf_eri[yx_offset];
                double kl_eri_yx = klyx_buf_eri[xy_offset];
                Vab_ijkl -= pfac_xy * (ij_r12_xy * kl_eri_xy + ij_r12_yx * kl_eri_yx);
		if (i != j)
		  Vab_jikl -= pfac_xy * (ij_r12_yx * kl_eri_xy + ij_r12_xy * kl_eri_yx);
		if (k != l)
		  Vab_ijlk -= pfac_xy * (ij_r12_xy * kl_eri_yx + ij_r12_yx * kl_eri_xy);
		if (i != j && k != l) {
		  Vaa_ijkl -= pfac_xy * (ij_r12_xy - ij_r12_yx)*(kl_eri_xy - kl_eri_yx);
		  Vab_jilk -= pfac_xy * (ij_r12_yx * kl_eri_yx + ij_r12_xy * kl_eri_xy);
		}
                double kl_r12_xy = klyx_buf_r12[yx_offset];
                double kl_r12_yx = klyx_buf_r12[xy_offset];
                Xab_ijkl -= pfac_xy * (ij_r12_xy * kl_r12_xy + ij_r12_yx * kl_r12_yx);
		if (i != j)
		  Xab_jikl -= pfac_xy * (ij_r12_yx * kl_r12_xy + ij_r12_xy * kl_r12_yx);
		if (k != l)
		  Xab_ijlk -= pfac_xy * (ij_r12_xy * kl_r12_yx + ij_r12_yx * kl_r12_xy);
		if (i != j && k != l) {
		  Xaa_ijkl -= pfac_xy * (ij_r12_xy - ij_r12_yx)*(kl_r12_xy - kl_r12_yx);
		  Xab_jilk -= pfac_xy * (ij_r12_yx * kl_r12_yx + ij_r12_xy * kl_r12_xy);
		}
                double kl_r12t1_xy = klyx_buf_r12t1[yx_offset];
                double kl_r12t1_yx = klyx_buf_r12t1[xy_offset];
                double lk_r12t1_xy = lkyx_buf_r12t1[yx_offset];
                double lk_r12t1_yx = lkyx_buf_r12t1[xy_offset];
                double kl_Tr12_xy = -kl_r12t1_xy-lk_r12t1_yx;
                double kl_Tr12_yx = -kl_r12t1_yx-lk_r12t1_xy;
                Tab_ijkl += pfac_xy * (ij_r12_xy * kl_Tr12_xy + ij_r12_yx * kl_Tr12_yx);
		if (i != j)
		  Tab_jikl += pfac_xy * (ij_r12_yx * kl_Tr12_xy + ij_r12_xy * kl_Tr12_yx);
		if (k != l)
		  Tab_ijlk += pfac_xy * (ij_r12_xy * kl_Tr12_yx + ij_r12_yx * kl_Tr12_xy);
		if (i != j && k != l) {
		  Taa_ijkl += pfac_xy * (ij_r12_xy - ij_r12_yx)*(kl_Tr12_xy - kl_Tr12_yx);
		  Tab_jilk += pfac_xy * (ij_r12_yx * kl_Tr12_yx + ij_r12_xy * kl_Tr12_xy);
		}
              }
	    }
            Vab_ij[kl_ab] += Vab_ijkl;
	    if (i != j)
	      Vab_ji[kl_ab] += Vab_jikl;
	    if (k != l)
	      Vab_ij[lk_ab] += Vab_ijlk;
	    if (i != j && k != l) {
	      Vaa_ij[kl_aa] += Vaa_ijkl;
	      Vab_ji[lk_ab] += Vab_jilk;
	    }
            Xab_ij[kl_ab] += Xab_ijkl;
	    if (i != j)
	      Xab_ji[kl_ab] += Xab_jikl;
	    if (k != l)
	      Xab_ij[lk_ab] += Xab_ijlk;
	    if (i != j && k != l) {
	      Xaa_ij[kl_aa] += Xaa_ijkl;
	      Xab_ji[lk_ab] += Xab_jilk;
	    }
            Tab_ij[kl_ab] += Tab_ijkl;
	    if (i != j)
	      Tab_ji[kl_ab] += Tab_jikl;
	    if (k != l)
	      Tab_ij[lk_ab] += Tab_ijlk;
	    if (i != j && k != l) {
	      Taa_ij[kl_aa] += Taa_ijkl;
	      Tab_ji[lk_ab] += Tab_jilk;
	    }
            tim_exit("MO ints contraction");

#if PRINT_R12_INTERMED
	    if (i != j && k != l)
	      printf("Vaa[%d][%d] = %lf\n",ij_aa,kl_aa,Vaa_ij[kl_aa]);
            printf("Vab[%d][%d] = %lf\n",ij_ab,kl_ab,Vab_ij[kl_ab]);
	    if (i != j)
	      printf("Vab[%d][%d] = %lf\n",ji_ab,kl_ab,Vab_ji[kl_ab]);
	    if (k != l)
	      printf("Vab[%d][%d] = %lf\n",ij_ab,lk_ab,Vab_ij[lk_ab]);
	    if (i != j && k != l)
	      printf("Vab[%d][%d] = %lf\n",ji_ab,lk_ab,Vab_ji[lk_ab]);
	    if (i != j && k != l)
	      printf("Xaa[%d][%d] = %lf\n",ij_aa,kl_aa,Xaa_ij[kl_aa]);
            printf("Xab[%d][%d] = %lf\n",ij_ab,kl_ab,Xab_ij[kl_ab]);
	    if (i != j)
	      printf("Xab[%d][%d] = %lf\n",ji_ab,kl_ab,Xab_ji[kl_ab]);
	    if (k != l)
	      printf("Xab[%d][%d] = %lf\n",ij_ab,lk_ab,Xab_ij[lk_ab]);
	    if (i != j && k != l)
	      printf("Xab[%d][%d] = %lf\n",ji_ab,lk_ab,Xab_ji[lk_ab]);
	    if (i != j && k != l)
	      printf("Taa[%d][%d] = %lf\n",ij_aa,kl_aa,Taa_ij[kl_aa]);
            printf("Tab[%d][%d] = %lf\n",ij_ab,kl_ab,Tab_ij[kl_ab]);
	    if (i != j)
	      printf("Tab[%d][%d] = %lf\n",ji_ab,kl_ab,Tab_ji[kl_ab]);
	    if (k != l)
	      printf("Tab[%d][%d] = %lf\n",ij_ab,lk_ab,Tab_ij[lk_ab]);
	    if (i != j && k != l)
	      printf("Tab[%d][%d] = %lf\n",ji_ab,lk_ab,Tab_ji[lk_ab]);
#endif
            r12intsacc->release_pair_block(i,j,R12IntsAcc::r12);
          }
        r12intsacc->release_pair_block(k,l,R12IntsAcc::eri);
        r12intsacc->release_pair_block(k,l,R12IntsAcc::r12);
        r12intsacc->release_pair_block(k,l,R12IntsAcc::r12t1);
	r12intsacc->release_pair_block(l,k,R12IntsAcc::r12t1);
      }
  }
  // Tasks that don't do any work here still need to create these timers
  tim_enter("MO ints retrieve");
  tim_exit("MO ints retrieve");
  tim_enter("MO ints contraction");
  tim_exit("MO ints contraction");

  delete[] proc_with_ints;
  tim_exit("mp2-r12a intermeds");
  ExEnv::out0() << indent << "End of computation of intermediates" << endl;
  r12intsacc->deactivate();

  
  // If running in distributed environment use Message group to collect intermediates and pair energies on node 0
  if (nproc > 1) {
    // collect contributions from the nodes that computed the intermediates and broadcast to all nodes
    msg->sum(Vaa_ijkl,naa*naa,0,-1);
    msg->sum(Vab_ijkl,nab*nab,0,-1);
    msg->sum(Xaa_ijkl,naa*naa,0,-1);
    msg->sum(Xab_ijkl,nab*nab,0,-1);
    msg->sum(Taa_ijkl,naa*naa,0,-1);
    msg->sum(Tab_ijkl,nab*nab,0,-1);

    int naa = emp2_aa.dim().n();
    int nab = emp2_ab.dim().n();
    double* epair_aa = new double[naa];
    double* epair_ab = new double[nab];
    bzerofast(epair_aa,naa);
    bzerofast(epair_ab,nab);
    emp2_aa.convert(epair_aa);
    emp2_ab.convert(epair_ab);
    msg->sum(epair_aa,naa,0,-1);
    msg->sum(epair_ab,nab,0,-1);
    msg->sync();
    emp2_aa.assign(epair_aa);
    emp2_ab.assign(epair_ab);
    delete[] epair_aa;
    delete[] epair_ab;
  }
  
  if (debug_)
    ExEnv::out0() << indent << "Gathered intermediates V, X, and T and MP2 pair energies" << endl;

  // Initialize global intermediates
  Vaa->unit();
  Vab->unit();
  Baa->unit();
  Bab->unit();
  Xaa.assign(0.0);
  Xab.assign(0.0);
  
  // Add intermediates contribution to their global values
  for(int ij=0;ij<naa;ij++)
    for(int kl=0;kl<=ij;kl++) {
      int ijkl = ij*naa+kl;
      int klij = kl*naa+ij;
      double velem = Vaa->get_element(ij,kl) + Vaa_ijkl[ijkl];
      Vaa->set_element(ij,kl,velem);
      if (ij != kl) {
        velem = Vaa->get_element(kl,ij) + Vaa_ijkl[klij];
        Vaa->set_element(kl,ij,velem);
      }
      double xelem = Xaa->get_element(ij,kl) + Xaa_ijkl[ijkl];
      Xaa->set_element(ij,kl,xelem);
      if (ij != kl) {
        xelem = Xaa->get_element(kl,ij) + Xaa_ijkl[klij];
        Xaa->set_element(kl,ij,xelem);
      }
      double belem = Baa->get_element(ij,kl) + 0.5*(Taa_ijkl[ijkl] + Taa_ijkl[klij]);
      Baa->set_element(ij,kl,belem);
      Baa->set_element(kl,ij,belem);
    }

  for(int ij=0;ij<nab;ij++)
    for(int kl=0;kl<=ij;kl++) {
      int ijkl = ij*nab+kl;
      int klij = kl*nab+ij;
      double velem = Vab->get_element(ij,kl) + Vab_ijkl[ijkl];
      Vab->set_element(ij,kl,velem);
      if (ij != kl) {
        velem = Vab->get_element(kl,ij) + Vab_ijkl[klij];
        Vab->set_element(kl,ij,velem);
      }
      double xelem = Xab->get_element(ij,kl) + Xab_ijkl[ijkl];
      Xab->set_element(ij,kl,xelem);
      if (ij != kl) {
        xelem = Xab->get_element(kl,ij) + Xab_ijkl[klij];
        Xab->set_element(kl,ij,xelem);
      }
      double belem = Bab->get_element(ij,kl) + 0.5*(Tab_ijkl[ijkl] + Tab_ijkl[klij]);
      Bab->set_element(ij,kl,belem);
      Bab->set_element(kl,ij,belem);
    }
  msg->sum(aoint_computed);

#if PRINT_BIGGEST_INTS
  biggest_ints_1.combine(msg);
  biggest_ints_2.combine(msg);
  biggest_ints_2s.combine(msg);
  biggest_ints_3a.combine(msg);
  biggest_ints_3.combine(msg);
#endif

#if PRINT_BIGGEST_INTS
  ExEnv::out0() << "biggest 1/4 transformed ints" << endl;
  for (int i=0; i<biggest_ints_1.ncontrib(); i++) {
    ExEnv::out0() << scprintf("%3d %3d %3d %3d %16.12f",
                              biggest_ints_1.indices(i)[0],
                              biggest_ints_1.indices(i)[1],
                              biggest_ints_1.indices(i)[2],
                              biggest_ints_1.indices(i)[3],
                              biggest_ints_1.val(i)
                              )
    << endl;
  }
  ExEnv::out0() << "biggest 2/4 transformed ints" << endl;
  for (int i=0; i<biggest_ints_2.ncontrib(); i++) {
    ExEnv::out0() << scprintf("%3d %3d %3d %3d %16.12f",
                              biggest_ints_2.indices(i)[0],
                              biggest_ints_2.indices(i)[1],
                              biggest_ints_2.indices(i)[2],
                              biggest_ints_2.indices(i)[3],
                              biggest_ints_2.val(i)
                              )
    << endl;
  }
  ExEnv::out0() << "restricted 2/4 transformed ints" << endl;
  for (int i=0; i<biggest_ints_2s.ncontrib(); i++) {
    ExEnv::out0() << scprintf("%3d %3d %3d %3d %16.12f",
                              biggest_ints_2s.indices(i)[0],
                              biggest_ints_2s.indices(i)[1],
                              biggest_ints_2s.indices(i)[2],
                              biggest_ints_2s.indices(i)[3],
                              biggest_ints_2s.val(i)
                              )
    << endl;
  }
  ExEnv::out0() << "biggest 3/4 transformed ints (in 3.)" << endl;
  for (int i=0; i<biggest_ints_3a.ncontrib(); i++) {
    ExEnv::out0() << scprintf("%3d %3d %3d %3d %16.12f",
                              biggest_ints_3a.indices(i)[0],
                              biggest_ints_3a.indices(i)[1],
                              biggest_ints_3a.indices(i)[2],
                              biggest_ints_3a.indices(i)[3],
                              biggest_ints_3a.val(i)
                              )
    << endl;
  }
  ExEnv::out0() << "biggest 3/4 transformed ints (in 4.)" << endl;
  for (int i=0; i<biggest_ints_3.ncontrib(); i++) {
  ExEnv::out0() << scprintf("%3d %3d %3d %3d %16.12f",
                            biggest_ints_3.indices(i)[0],
                            biggest_ints_3.indices(i)[1],
                            biggest_ints_3.indices(i)[2],
                            biggest_ints_3.indices(i)[3],
                            biggest_ints_3.val(i)
                            )
    << endl;
  }
#endif

  // Print out various energies etc.

  if (debug_) {
    ExEnv::out0() << indent << "Number of shell quartets for which AO integrals\n"
    << indent << "would have been computed without bounds checking: "
    << npass*nshell*nshell*(nshell+1)*(nshell+1)/2
    << endl;

    ExEnv::out0() << indent << "Number of shell quartets for which AO integrals\n"
    << indent << "were computed: " << aoint_computed
    << endl;
  }
  
  /*--------
    Cleanup
   --------*/
  delete[] Vaa_ijkl;
  delete[] Taa_ijkl;
  delete[] Xaa_ijkl;
  delete[] Vab_ijkl;
  delete[] Tab_ijkl;
  delete[] Xab_ijkl;
  
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

  if (me == 0 && wfn->if_to_checkpoint()) {
    StateOutBin stateout(wfn->checkpoint_file());
    SavableState::save_state(wfn,stateout);
    ExEnv::out0() << indent << "Checkpointed the wave function" << endl;
  }
  
  return;
}



///////////////////////////////////////////////////////
// Compute the batchsize for the transformation
//
// Only arrays allocated before exiting the loop over
// i-batches are included here  - only these arrays
// affect the batch size.
///////////////////////////////////////////////////////
int
R12IntEval_sbs_A::compute_transform_batchsize_(size_t mem_alloc, size_t mem_static, int nocc_act, const int num_te_types)
{
  // Check is have enough for even static objects
  size_t mem_dyn = 0;
  if (mem_alloc <= mem_static)
    return 0;
  else
    mem_dyn = mem_alloc - mem_static;

  // Determine if calculation is possible at all (i.e., if ni=1 possible)
  int ni = 1;
  distsize_t maxdyn = compute_transform_dynamic_memory_(ni, nocc_act, num_te_types);
  if (maxdyn > mem_dyn) {
    return 0;
  }

  ni = 2;
  while (ni<=nocc_act) {
    maxdyn = compute_transform_dynamic_memory_(ni, nocc_act, num_te_types);
    if (maxdyn >= mem_dyn) {
      ni--;
      break;
    }
    ni++;
  }
  if (ni > nocc_act) ni = nocc_act;

  return ni;
}

//////////////////////////////////////////////////////
// Compute required (dynamic) memory
// for a given batch size of the transformation
//
// Only arrays allocated before exiting the loop over
// i-batches are included here  - only these arrays
// affect the batch size.
//////////////////////////////////////////////////////
distsize_t
R12IntEval_sbs_A::compute_transform_dynamic_memory_(int ni, int nocc_act, const int num_te_types)
{
  int nproc = r12info()->msg()->n();

  ///////////////////////////////////////
  // the largest memory requirement will
  // occur just before
  // the end of the i-batch loop (mem)
  ///////////////////////////////////////

  // compute nij as nij on node 0, since nij on node 0 is >= nij on other nodes
  int index = 0;
  int nij = 0;
  for (int i=0; i<ni; i++) {
    for (int j=0; j<nocc_act; j++) {
      if (index++ % nproc == 0) nij++;
    }
  }

  int nbasis = r12info()->basis()->nbasis();
  int nfuncmax = r12info()->basis()->max_nfunction_in_shell();
  int nthread = r12info()->thr()->nthread();

  distsize_t memsize = sizeof(double)*(num_te_types*((distsize_t)nthread * ni * nbasis * nfuncmax * nfuncmax // iqrs
						     + (distsize_t)nij * 2 * nbasis * nfuncmax  // iqjs and iqjr buffers
						     + (distsize_t)nij * nbasis * nbasis // iqjs_contrib - buffer of half and higher
						     // transformed integrals
						     )
				       );
  return memsize;
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
