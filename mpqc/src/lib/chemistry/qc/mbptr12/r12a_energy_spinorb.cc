//
// r12a_energy_spinorb.cc
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

#include <stdlib.h>
#include <math.h>
#include <limits.h>

#include <scconfig.h>
#include <util/misc/formio.h>
#include <util/misc/timer.h>
#include <util/group/memory.h>
#include <util/group/message.h>
#include <util/class/class.h>
#include <util/state/state.h>
#include <util/state/state_text.h>
#include <util/state/state_bin.h>
#include <math/scmat/matrix.h>
#include <math/scmat/blocked.h>
#include <math/scmat/repl.h>
#include <math/scmat/local.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/basis/integral.h>
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

using namespace std;
using namespace sc;

#define SINGLE_THREAD_E12   0
#define SINGLE_THREAD_QBT34 0

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

/*-------------------------------------
  Based on MBPT2::compute_mp2_energy()
 -------------------------------------*/
void
MBPT2_R12::compute_r12a_spinorb_memgrp_()
{
  double escf = 0.0;
  double emp2 = 0.0;
  if (mp2_done_) {
    escf = hf_energy_;
    emp2 = mp2_corr_energy_;
  }
  else {
    escf = ref_energy();
    emp2 = MBPT2::corr_energy();
  }

  // Basis Sets
  Ref<GaussianBasisSet> bs = basis();
  Ref<GaussianBasisSet> bs_aux = aux_basis();

  Ref<SCMatrixKit> kit = bs->matrixkit();

  // log2 of the erep tolerance
  // (erep < 2^tol => discard)
  const int tol = (int) (-10.0/log10(2.0));  // discard ints smaller than 10^-20

  int nij;        // number of i,j pairs on a node (for e.g., mo_int)
  const int num_te_types = 3;
  enum te_types {eri=0, r12=1, r12t1=2};
  double *mo_int; // MO integrals of type (og|og)
  double *integral_ijsq; // half-transformed integrals

  int nbasis = bs->nbasis();
  int nbasis_aux = bs_aux->nbasis();
  int te_type;
  int nocc_act;
  int i, j;
  int ii, bb;
  int x, y;
  int a, b;
  int nshell;
  int offset;
  int ik_offset;
  int i_offset; 
  int npass, pass;
  int tmpint;
  int np, nq, nr, ns; 
  int P, Q, R, S;
  int p, q, r, s;
  int bf1, bf2, bf3, bf4;
  int index;
  int me;
  int nproc;
  int rest;
  int p_offset, q_offset, r_offset, s_offset;

  int aoint_computed = 0; 
  int xyz;
  int int_index;
  size_t mem_static;    // static memory in bytes
  int ij_proc;          // the processor which has ij pair
  int ij_index;         // of the ij pairs on a proc, this ij pair is number ij_index
                        // (i.e., ij_index < nij)
  int ik_proc;          // the processor which has ik pair
  int ik_index;
  int jloop, kloop;

  int ni;

  double *evals;              // scf eigenvalues
  double *iajb_ptr, *ibja_ptr, *iakb_ptr, *ibka_ptr;
  double *iajc_ptr, *ibjc_ptr, *icjb_ptr, *icja_ptr;
  double *ijkb_ptr, *ibkj_ptr;
  double pqrs;
  double *c_sa, c_rj;
  double *c_pi, *c_qi, *c_sj;
  double *c_qx, *c_qa, *c_sb, *c_pa, *c_pq, *c_sy;
  double delta_ijab, delta_ijbc, delta_ijac;
  double er12a=0.0;
  double emp2r12a=0.0;
  double *mo_intbuf;          // buffer used for sending mo integrals
  double tmpval, tmpval1;

  double *ijsx_tmp;      // three-quarter transformed two-el integrals
  double *ijyx_tmp;      // three-quarter transformed two-el integrals
  double *integral_ijsx;  // all three-quarter transformed two-el integrals
  double *integral_ijyx; // mo integrals (y = any MO)
  double *integral_ijsq_ptr;
  double *ijyx_ptr;
  double *ijsx_ptr;

  BiggestContribs biggest_coefs(5,10);
  CharacterTable ct = molecule()->point_group()->char_table();

#if PRINT_BIGGEST_INTS
  BiggestContribs biggest_ints_2(4,40);
  BiggestContribs biggest_ints_2s(4,40);
  BiggestContribs biggest_ints_3a(4,40);
  BiggestContribs biggest_ints_3(4,40);
#endif

  tim_enter("r12a-mem");

  nfuncmax = basis()->max_nfunction_in_shell();
  nshell = basis()->nshell();

  me = msg_->me();
  nproc = msg_->n();

  ExEnv::out0() << endl << indent
	       << "Entered memgrp based MP2-R12/A routine" << endl;
  ExEnv::out0() << indent << scprintf("nproc = %i", nproc) << endl;

  nocc = 0;
  for (i=0; i<oso_dimension()->n(); i++) {
    if (reference_->occupation(i) == 2.0) nocc++;
  }

  nocc_act = nocc - nfzc;
  nvir  = noso - nocc;

  // Do a few preliminary tests to make sure the desired calculation
  // can be done (and appears to be meaningful!)

  if (nocc_act <= 0) {
    ExEnv::err0() << "There are no active occupied orbitals; program exiting" << endl;
    abort();
  }

  if (nvir <= 0) {
    ExEnv::err0() << "There are no active virtual orbitals; program exiting" << endl;
    abort();
  }
    
  if (restart_orbital_r12_memgrp_) {
    ExEnv::out0() << indent
		 << scprintf("Restarting at orbital %d",
			     restart_orbital_r12_memgrp_)
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
  if (me == 0) {
    mem_static = nbasis*noso; // scf vector
    mem_static *= sizeof(double);
    int nthreads = thr_->nthread();
    mem_static += nthreads * integral()->storage_required(&Integral::grt,basis()); // integral evaluators
    ni = compute_r12atransform_batchsize_(mem_static, nocc_act-restart_orbital_r12_memgrp_, num_te_types); 
  }

  if (max_norb_ > 0 && ni > max_norb_) {
    ExEnv::out0() << indent
		 << "\"max_norb\" set: could have done "
		 << ni << " orbitals per pass otherwise."
		 << endl;
    ni = max_norb_;
  }

  // Send value of ni and mem_static to other nodes
  msg_->bcast(ni);
  double mem_static_double = (double) mem_static;
  msg_->bcast(mem_static_double);
  mem_static = (size_t) mem_static_double;
  
  // Compute the storage remaining for the integral routines
  size_t dyn_mem = distsize_to_size(compute_r12atransform_dynamic_memory_(ni,nocc_act,num_te_types));
  size_t mem_remaining = 0;
  if (mem_alloc <= (dyn_mem + mem_static)) mem_remaining += 0;
  else mem_remaining += mem_alloc - dyn_mem - mem_static;
  mem_remaining += thr_->nthread() * integral()->storage_required(&Integral::grt,basis());
  
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
	       << compute_r12atransform_dynamic_memory_(nocc_act,nocc_act,num_te_types)+mem_static
	       << " Bytes"
	       << endl;
  ExEnv::out0() << indent
	       << "Minimum memory required:        "
	       << compute_r12atransform_dynamic_memory_(1,nocc_act,num_te_types)+mem_static
	       << " Bytes"
	       << endl;
  ExEnv::out0() << indent
	       << "Batch size:                     " << ni
	       << endl;
  
  if (ni == 0) {
    ExEnv::err0() << "Batch size is 0: more memory or processors are needed"
		 << endl;
    abort();
  }
  
  if (dynamic_) {
    ExEnv::out0() << indent << "Using dynamic load balancing." << endl;
  }
  
  if (ni == nocc_act-restart_orbital_r12_memgrp_) {
    npass = 1;
    rest = 0;
  }
  else {
    rest = (nocc_act-restart_orbital_r12_memgrp_)%ni;
    npass = (nocc_act-restart_orbital_r12_memgrp_ - rest)/ni + 1;
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
  index = 0;
  for (i=0; i<ni; i++) {
    for (j=0; j<nocc_act; j++) {
      if (index++ % nproc == me) nijmax++;
    }
  }
  

  ////////////////////////////////////////////////
  // The scf vector is distributed between nodes;
  // put a copy of the scf vector on each node;
  ////////////////////////////////////////////////

  RefDiagSCMatrix occ;
  RefSCMatrix Scf_Vec;
  RefDiagSCMatrix evalmat;
  eigen(evalmat, Scf_Vec, occ);

  if (debug_ > 1) {
    evalmat.print("eigenvalues");
    Scf_Vec.print("eigenvectors");
  }

  double *scf_vector_dat = new double[nbasis*noso];
  Scf_Vec.t()->convert(scf_vector_dat);

  evals = new double[noso];
  double** scf_vector = new double*[nbasis];
  for (i=0; i<nbasis; i++) {
    scf_vector[i] = &scf_vector_dat[i*noso];
    }
  for (i=0; i<noso; i++) {
      evals[i] = evalmat(i);
    }

  Scf_Vec = 0;
  evalmat = 0;

  if (debug_ > 2 && me == 0) {
    for (j=0; j<noso; j++) {
      ExEnv::out0() << indent
           << scprintf("eigenvalue[%3d] = %15.10lf", j, evals[j]);
      if (j < nfzc) ExEnv::out0() << " (frozen docc)";
      else if (j < nocc_act + nfzc) ExEnv::out0() << " (active docc)";
      else if (j < nvir + nocc_act + nfzc) ExEnv::out0() << " (active uocc)";
      ExEnv::out0() << endl;
    }
  }

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
  integral()->set_storage(mem_remaining);
  tbints_ = new Ref<TwoBodyInt>[thr_->nthread()];
  for (i=0; i<thr_->nthread(); i++) {
    tbints_[i] = integral()->grt();
  }
  ExEnv::out0() << indent
	       << scprintf("Memory used for integral storage:       %i Bytes",
			   integral()->storage_used())
	       << endl;


  if (mem.null()) {
    ExEnv::errn() << "MBPT2: memory group not initialized" << endl;
    abort();
  }
  mem->set_localsize(num_te_types*nijmax*nbasis*nbasis*sizeof(double));
  ExEnv::out0() << indent
	       << "Size of global distributed array:       "
	       << mem->totalsize()
	       << " Bytes"
	       << endl;
  MemoryGrpBuf<double> membuf_remote(mem);


  Ref<ThreadLock> lock = thr_->new_lock();
  R12A_GRT_12Qtr** e12thread = new R12A_GRT_12Qtr*[thr_->nthread()];
  for (i=0; i<thr_->nthread(); i++) {
    e12thread[i] = new R12A_GRT_12Qtr(i, thr_->nthread(), me, nproc,
				      mem, msg_, lock, basis(), aux_basis(), tbints_[i],
				      nocc, nocc_act, scf_vector, tol, debug_,
				      dynamic_);
  }

  if (npass > 1 || restart_orbital_r12_memgrp_) {
    bool restart = restart_orbital_r12_memgrp_;
    // using File integrals accumulator when npass > 1 or restarting
#if HAVE_MPIIO
    ExEnv::out0() << indent << "Will use MPI-IO (individual I/O) to handle transformed integrals" << endl;
    r12intsacc_ = new R12IntsAcc_MPIIOFile_Ind(mem,r12ints_file_,num_te_types,nbasis,nbasis,nocc,nfzc,restart);
#else
    ExEnv::out0() << indent << "Will use POSIX I/O on node 0 to handle transformed integrals" << endl;
    r12intsacc_ = new R12IntsAcc_Node0File(mem,r12ints_file_,num_te_types,nbasis,nbasis,nocc,nfzc,restart);
#endif
  }
  else {
    // using MemoryGrp integrals accumulator when npass = 1 and not restarting
//    ExEnv::out0() << indent << "Will use MPI-IO (individual I/O) to handle transformed integrals" << endl;
//    bool restart = restart_orbital_r12_memgrp_;
//    r12intsacc_ = new R12IntsAcc_MPIIOFile_Ind(mem,r12ints_file_,num_te_types,nbasis,nbasis,nocc,nfzc,restart);
    ExEnv::out0() << indent << "Will hold transformed integrals in memory" << endl;
    r12intsacc_ = new R12IntsAcc_MemoryGrp(mem,num_te_types,nbasis,nbasis,nocc,nfzc);
  }


  /*-----------------------------------

    Start the integrals transformation

   -----------------------------------*/
  tim_enter("mp2-r12/a passes");
  if (me == 0 && if_to_checkpoint()) {
    StateOutBin stateout(checkpoint_file());
    SavableState::save_state(this,stateout);
  }
  for (pass=0; pass<npass; pass++) {

    ExEnv::out0() << indent << "Beginning pass " << pass+1 << endl;

    i_offset = restart_orbital_r12_memgrp_ + pass*ni + nfzc;
    if ((pass == npass - 1) && (rest != 0)) ni = rest;

    // Compute number of of i,j pairs on each node during current pass for
    // two-el integrals
    index = 0;
    nij = 0;
    for (i=0; i<ni; i++) {
      for (j=0; j<nocc_act; j++) {
        if (index++ % nproc == me) nij++;
      }
    }

    // debug print
    if (debug_)
      ExEnv::outn() << indent << "node " << me << ", nij = " << nij << endl;
    // end of debug print

    //mem->sync(); // This must be here or gamma non-sep will be wrong when running
                 // on multiple processors with more than one pass

    r_offset = 0;

    // Allocate and initialize some arrays
    // (done here to avoid having these arrays
    // overlap with arrays allocated later)

    // Allocate (and initialize) some arrays

    integral_ijsq = (double*) mem->localdata();
    bzerofast(integral_ijsq, (num_te_types*nij*nbasis*nbasis));
    integral_ijsq = 0;
    mem->sync();
    ExEnv::out0() << indent
		   << scprintf("Begin loop over shells (grt, 1.+2. q.t.)") << endl;

    // Do the two electron integrals and the first two quarter transformations
    tim_enter("grt+1.qt+2.qt");
    for (i=0; i<thr_->nthread(); i++) {
      e12thread[i]->set_i_offset(i_offset);
      e12thread[i]->set_ni(ni);
      thr_->add_thread(i,e12thread[i]);
#     if SINGLE_THREAD_E12
      e12thread[i]->run();
#     endif
      }
#   if !SINGLE_THREAD_E12
    thr_->start_threads();
    thr_->wait_threads();
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
	    for(te_type=0; te_type<PRINT_NUM_TE_TYPES; te_type++,integral_ij_offset+=nbasis*nbasis) {
	      for (int s = 0; s<nbasis; s++) {
		double *integral_ijsq_ptr = integral_ij_offset + s*nbasis;
		for (int q = 0; q<nbasis; q++) {
		  printf("2Q: (%d %d|%d %d) = %12.8f\n",
		       i,q,j+j_offset,s,*integral_ijsq_ptr);
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
    ijsx_tmp = new double[nbasis];

    ExEnv::out0() << indent << "Begin third q.t." << endl;
    tim_enter("3. q.t.");
    // Begin third quarter transformation;
    // generate (ix|js) for i act, j act, and x any MO
    index = 0;
    ij_index = 0;
    for (i=0; i<ni; i++) {
      for (j=0; j<nocc_act; j++) {
        if (index++ % nproc == me) {
	  double *integral_ij_offset = integral_ijsq + num_te_types*nbasis*nbasis*ij_index;
	  for(te_type=0; te_type<num_te_types; te_type++,integral_ij_offset+=nbasis*nbasis) {
	    for (s=0; s<nbasis; s++) {
	      integral_ijsq_ptr = integral_ij_offset + s*nbasis;
	      bzerofast(ijsx_tmp, nbasis);
	      for (q=0; q<nbasis; q++,integral_ijsq_ptr++) {
		ijsx_ptr = ijsx_tmp;
		c_qx = scf_vector[q];
		tmpval = *integral_ijsq_ptr;
#if PRINT_BIGGEST_INTS
		biggest_ints_2.insert(tmpval,i+i_offset,j,s,q);
		if ((i+i_offset==104 && j == 1)
		    ||(i+i_offset==104 && j == 2)) {
		  biggest_ints_2s.insert(tmpval,i+i_offset,j,s,q);
		}
#endif
		for (x=0; x<noso; x++) {
		  *ijsx_ptr++ += *c_qx++ * tmpval;
		}
              }   // exit q loop

	      // Put ixjs into integral_ijsq, while overwriting what was there;
	      // i.e., integral_ijsq will now contain three-quarter transformed
	      // integrals ijsx
	      integral_ijsq_ptr = integral_ij_offset + s*nbasis;
	      ijsx_ptr = ijsx_tmp;
	      for (x=0; x<noso; x++) {
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
    integral_ijsx = integral_ijsq;

#if PRINT3Q
    if (me == 0) {
      int index = 0;
      int ij_index = 0;
      int j_offset = nfzc;
      for (int i = 0; i<ni; i++) {
	for (int j = 0; j<nocc_act; j++) {
	  if (index++ % nproc == me) {
	    double *integral_ij_offset = integral_ijsx + num_te_types*nbasis*nbasis*ij_index;
	    for(te_type=0; te_type<PRINT_NUM_TE_TYPES; te_type++,integral_ij_offset+=nbasis*nbasis) {
	      for (int s = 0; s<nbasis; s++) {
		double *integral_ijsx_ptr = integral_ij_offset + s*nbasis;
		for (int x = 0; x<noso; x++) {
		  printf("3Q: (%d %d|%d %d) = %12.8f\n",
		       i,x,j+j_offset,s,*integral_ijsx_ptr);
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

    ijyx_tmp = new double[noso];
    // in ijyx: i act; x any MO, j act; y any MO.

    // Begin fourth quarter transformation
    ExEnv::out0() << indent << "Begin fourth q.t." << endl;
    tim_enter("4. q.t.");
    index = 0;
    ij_index = 0;
    for (i=0; i<ni; i++) {
      for (j=0; j<nocc_act; j++) {
        if (index++ % nproc == me) {
	  double *integral_ij_offset = integral_ijsx + num_te_types*nbasis*nbasis*ij_index;
	  for(te_type=0; te_type<num_te_types; te_type++,integral_ij_offset+=nbasis*nbasis) {
	    for (x=0; x<noso; x++) {
	      bzerofast(ijyx_tmp, noso);
	      ijsx_ptr = integral_ij_offset + x;
	      for (s=0; s<nbasis; s++) {
		c_sy = scf_vector[s];
		ijyx_ptr = ijyx_tmp;
		tmpval = *ijsx_ptr;
#if PRINT_BIGGEST_INTS
		biggest_ints_3.insert(tmpval,i+i_offset,j,s,x);
		if ((i+i_offset==105 && j == 2 && s == 170 && x == 3)
		    ||(i+i_offset==102 && j == 2 && s == 170 && x == 2)) {
		  ExEnv::outn() << scprintf("3/4: %3d %3d %3d %3d: %16.10f",
					   i+i_offset, j, s, x-nocc)
			       << endl;
                }
#endif
		for (y=0; y<noso; y++) {
		  *ijyx_ptr++ += *c_sy++ * tmpval;
		} // exit y loop
		ijsx_ptr += nbasis;
	      }   // exit s loop
	      // Put ijyx_tmp into ijsx for one i,x,j while
	      // overwriting elements of ijsx
	      ijsx_ptr = integral_ij_offset + x;
	      ijyx_ptr = ijyx_tmp;
	      for (y=0; y<noso; y++) {
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
    mo_int = integral_ijsx;
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
	  for(te_type=0; te_type<PRINT_NUM_TE_TYPES; te_type++,integral_ij_offset+=nbasis*nbasis) {
	    for (int y = 0; y < noso; y++) {
	      double *integral_ijyx_ptr = integral_ij_offset + y*nbasis;
	      for (int x = 0; x<noso; x++) {
		if (( symorb_irrep_[i+i_offset] ^
		       symorb_irrep_[j+j_offset] ^
		       symorb_irrep_[x] ^
		       symorb_irrep_[y]) ) {
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
	    for(te_type=0; te_type<PRINT_NUM_TE_TYPES; te_type++,integral_ij_offset+=nbasis*nbasis) {
	      for (int y = 0; y < noso; y++) {
		double *integral_ijyx_ptr = integral_ij_offset + y*nbasis;
		for (int x = 0; x<noso; x++) {
		  printf("4Q: type = %d (%d %d|%d %d) = %12.8f\n",
		       te_type,i,x,j+j_offset,y,*integral_ijyx_ptr);
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
    tim_enter("compute ecorr");

    index = 0;
    ij_index = 0;
    for (i=0; i<ni; i++) {
      double ecorr_i = 0.0;
      for (j=0; j<nocc_act; j++) {
        double ecorr_ij = 0.0;
        if (index++ % nproc == me) {

	  for (b=0; b<nvir; b++) {
	    iajb_ptr = &mo_int[nocc + nbasis*(b+nocc + nbasis*num_te_types*ij_index)];
	    ibja_ptr = &mo_int[b+nocc + nbasis*(nocc + nbasis*num_te_types*ij_index)];
	    for (a=0; a<nvir; a++) {
#if PRINT4Q_MP2
	      printf("4Q: (%d %d|%d %d) = %12.8f\n",
		     i,a+nocc,j+nfzc,b+nocc,*iajb_ptr);
#endif
	      delta_ijab = evals[i_offset+i]+evals[j+nfzc]-evals[nocc+a]-evals[nocc+b];
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
	      tmpval = *iajb_ptr*(2**iajb_ptr - *ibja_ptr)/delta_ijab;
	      er12a += tmpval;
	      if (debug_) ecorr_ij += tmpval;
	      iajb_ptr++;
	      ibja_ptr += nbasis;
	    } // exit a loop
	  }   // exit b loop
	  ij_index++;
	}     // endif
	if (debug_) {
	  msg_->sum(ecorr_ij);
	  ecorr_i += ecorr_ij;
	  ExEnv::out0() << indent
			<< scprintf("correlation energy for pair %3d %3d = %16.12f",
				    i+i_offset, j, ecorr_ij)
			<< endl;
	}
      }         // exit j loop
      if (debug_) {
	ExEnv::out0() << indent
		     << scprintf("correlation energy for orbital %3d = %16.12f",
				 i+i_offset, ecorr_i)
		     << endl;
      }
    }           // exit i loop
    tim_exit("compute ecorr");
  
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
    r12intsacc_->store_memorygrp(mem,ni);
    tim_exit("MO ints store");
    //restart_orbital_r12_memgrp_ += ni;
    current_orbital_ += ni;
    mem->sync();

    if (me == 0 && if_to_checkpoint()) {
      StateOutBin stateout(checkpoint_file());
      SavableState::save_state(this,stateout);
    }

  } // end of loop over passes
  tim_exit("mp2-r12/a passes");
  if (debug_)
    ExEnv::out0() << indent << "End of mp2-r12/a transformation" << endl;
  // Done storing integrals - commit the content
  // WARNING: it is not safe to use mem until deactivate has been called on the accumulator
  //          After that deactivate the size of mem will be 0 [mem->set_localsize(0)]
  r12intsacc_->commit();

  /*--------------------------------
    Compute MP2-R12/A intermediates
    and collect on node0
   --------------------------------*/
  tim_enter("mp2-r12a intermeds");
  int naa = (nocc_act*(nocc_act-1))/2;          // Number of alpha-alpha pairs
  int nab = nocc_act*nocc_act;                  // Number of alpha-beta pairs
  if (debug_) {
    ExEnv::out0() << indent << "naa = " << naa << endl;
    ExEnv::out0() << indent << "nab = " << nab << endl;
  }
  double *Vaa_ijkl = new double[naa*naa];
  double *Taa_ijkl = new double[naa*naa];
  double *Vab_ijkl = new double[nab*nab];
  double *Tab_ijkl = new double[nab*nab];
  if (debug_)
    ExEnv::out0() << indent << "Allocated intermediates V and T" << endl;
  bzerofast(Vaa_ijkl,naa*naa);
  bzerofast(Taa_ijkl,naa*naa);
  bzerofast(Vab_ijkl,nab*nab);
  bzerofast(Tab_ijkl,nab*nab);

  // Compute intermediates
  if (debug_)
    ExEnv::out0() << indent << "Ready to compute intermediates V and T" << endl;
  const int pair_block_size = num_te_types*nbasis*nbasis;
  const double oosqrt2 = 1.0/sqrt(2.0);
  // Compute the number of tasks that have full access to the integrals
  // and split the work among them
  int nproc_with_ints = 0;
  for(int proc=0;proc<nproc;proc++)
    if (r12intsacc_->has_access(proc)) nproc_with_ints++;
  int *proc_with_ints = new int[nproc];
  int count = 0;
  for(int proc=0;proc<nproc;proc++)
    if (r12intsacc_->has_access(proc)) {
      proc_with_ints[proc] = count;
      count++;
    }
    else
      proc_with_ints[proc] = -1;
  if (debug_)
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

  if (r12intsacc_->has_access(me)) {
    int kl = 0;
    int kl_aa=-1;
    for(int k=0;k<nocc_act;k++)
      for(int l=0;l<=k;l++,kl++) {
        double pfac_kl = (k==l) ? oosqrt2 : 1.0;
        int kl_proc = kl%nproc_with_ints;
        if (kl_proc != proc_with_ints[me])
          continue;
	if (k != l)
	  ++kl_aa;
	int kl_ab = k*nocc_act + l;
	int lk_ab = l*nocc_act + k;
        
        // Get (|r12|) and (|[r12,T1]|) integrals only
        tim_enter("MO ints retrieve");
        double *klyx_buf_eri, *klyx_buf_r12t1;
        klyx_buf_eri = r12intsacc_->retrieve_pair_block(k,l,R12IntsAcc::eri);
        klyx_buf_r12t1 = r12intsacc_->retrieve_pair_block(k,l,R12IntsAcc::r12t1);
        
        double *lkyx_buf_r12t1;
        if (k!=l) {
          // Get (lx|[r12,T1]|ky)=(ky|[r12,T2]|lx) integrals only
          lkyx_buf_r12t1 = r12intsacc_->retrieve_pair_block(l,k,R12IntsAcc::r12t1);
	}
        else
          lkyx_buf_r12t1 = klyx_buf_r12t1;
        tim_exit("MO ints retrieve");

	int ij = 0;
        int ij_aa=-1;
        for(i=0;i<nocc_act;i++)
          for(j=0;j<=i;j++,ij++) {

            double pfac_ij = (i==j) ? oosqrt2 : 1.0;
	    if (i != j)
	      ++ij_aa;
	    int ij_ab = i*nocc_act + j;
	    int ji_ab = j*nocc_act + i;
            
            tim_enter("MO ints retrieve");
            double *ijyx_buf_r12 = r12intsacc_->retrieve_pair_block(i,j,R12IntsAcc::r12);
            tim_exit("MO ints retrieve");

	    double *Vaa_ij = Vaa_ijkl + ij_aa*naa;
            double *Vab_ij = Vab_ijkl + ij_ab*nab;
            double *Vab_ji = Vab_ijkl + ji_ab*nab;
            double *Taa_ij = Taa_ijkl + ij_aa*naa;
            double *Tab_ij = Tab_ijkl + ij_ab*nab;
            double *Tab_ji = Tab_ijkl + ji_ab*nab;
            
            tim_enter("MO ints contraction");
            if (ij == kl) {
              Vab_ij[ij_ab] = 1.0;
	      if (i != j) {
		Vaa_ij[ij_aa] = 1.0;
		Vab_ji[ji_ab] = 1.0;
	      }
              Tab_ij[ij_ab] = 1.0;
	      if (i != j) {
		Taa_ij[ij_aa] = 1.0;
		Tab_ji[ji_ab] = 1.0;
	      }
            }
            double Vaa_ijkl, Vab_ijkl, Vab_jikl, Vab_ijlk, Vab_jilk;
	    double Taa_ijkl, Tab_ijkl, Tab_jikl, Tab_ijlk, Tab_jilk;
	    Vaa_ijkl = Vab_ijkl = Vab_jikl = Vab_ijlk = Vab_jilk = 0.0;
	    Taa_ijkl = Tab_ijkl = Tab_jikl = Tab_ijlk = Tab_jilk = 0.0;
            for(y=0;y<noso;y++)
              for(x=0;x<=y;x++) {
                double pfac_xy = (x==y) ? 0.5 : 1.0;
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
            Vab_ij[kl_ab] += Vab_ijkl;
	    if (i != j)
	      Vab_ji[kl_ab] += Vab_jikl;
	    if (k != l)
	      Vab_ij[lk_ab] += Vab_ijlk;
	    if (i != j && k != l) {
	      Vaa_ij[kl_aa] += Vaa_ijkl;
	      Vab_ji[lk_ab] += Vab_jilk;
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
	      printf("Taa[%d][%d] = %lf\n",ij_aa,kl_aa,Taa_ij[kl_aa]);
            printf("Tab[%d][%d] = %lf\n",ij_ab,kl_ab,Tab_ij[kl_ab]);
	    if (i != j)
	      printf("Tab[%d][%d] = %lf\n",ji_ab,kl_ab,Tab_ji[kl_ab]);
	    if (k != l)
	      printf("Tab[%d][%d] = %lf\n",ij_ab,lk_ab,Tab_ij[lk_ab]);
	    if (i != j && k != l)
	      printf("Tab[%d][%d] = %lf\n",ji_ab,lk_ab,Tab_ji[lk_ab]);
#endif
            r12intsacc_->release_pair_block(i,j,R12IntsAcc::r12);
          }
        r12intsacc_->release_pair_block(k,l,R12IntsAcc::eri);
        r12intsacc_->release_pair_block(k,l,R12IntsAcc::r12t1);
        if (k != l) {
          r12intsacc_->release_pair_block(l,k,R12IntsAcc::r12t1);
        }
      }
  }
  delete[] proc_with_ints;
  tim_exit("mp2-r12a intermeds");
  r12intsacc_->deactivate();
  if (debug_)
    ExEnv::out0() << indent << "Computed V and T matrices" << endl;

  if (nproc > 1) {
    // Use MemoryGrp to send all contributions to intermediates V and T to node 0
    msg_->sum(Vaa_ijkl,naa*naa,0,0);
    msg_->sum(Vab_ijkl,nab*nab,0,0);
    msg_->sum(Taa_ijkl,naa*naa,0,0);
    msg_->sum(Tab_ijkl,nab*nab,0,0);
  }

  if (debug_)
    ExEnv::out0() << indent << "Gathered V and T matrices on node 0" << endl;

  //
  // Final stretch: compute the MP2-R12/A energy
  //
  tim_enter("mp2-r12a energy");

  // Convert V's into local SCMatrices on node 0
  // Then compute B (local SCMatrices as well)
  Ref<LocalSCMatrixKit> local_matrix_kit = new LocalSCMatrixKit();
  Ref<SCMatrix> Vaa, Vab, Baa, Bab;
  Ref<SCVector> epair_aa, epair_ab;
  RefSCDimension npair_ab = new SCDimension(nocc_act*nocc_act);
  RefSCDimension npair_aa = new SCDimension((nocc_act*(nocc_act-1))/2);
  if (me == 0) {
    Vaa = local_matrix_kit->matrix(npair_aa,npair_aa);
    Baa = local_matrix_kit->matrix(npair_aa,npair_aa);
    Vab = local_matrix_kit->matrix(npair_ab,npair_ab);
    Bab = local_matrix_kit->matrix(npair_ab,npair_ab);

    if (debug_)
      ExEnv::out0() << indent << "Allocated V and B" << endl;
    
    for(int ij=0;ij<naa;ij++)
      for(int kl=0;kl<=ij;kl++) {
	int ijkl = ij*naa+kl;
	int klij = kl*naa+ij;
	Vaa->set_element(ij,kl,Vaa_ijkl[ijkl]);
	Vaa->set_element(kl,ij,Vaa_ijkl[klij]);
	double belem = 0.5*(Taa_ijkl[ijkl] + Taa_ijkl[klij]);
	Baa->set_element(ij,kl,belem);
	Baa->set_element(kl,ij,belem);
      }
    
    for(int ij=0;ij<nab;ij++)
      for(int kl=0;kl<=ij;kl++) {
	int ijkl = ij*nab+kl;
	int klij = kl*nab+ij;
	Vab->set_element(ij,kl,Vab_ijkl[ijkl]);
	Vab->set_element(kl,ij,Vab_ijkl[klij]);
	double belem = 0.5*(Tab_ijkl[ijkl] + Tab_ijkl[klij]);
	Bab->set_element(ij,kl,belem);
	Bab->set_element(kl,ij,belem);
      }
    
    if (debug_) {
      Vaa->print("Alpha-alpha V matrix");
      Baa->print("Alpha-alpha B matrix");
      Vab->print("Alpha-beta V matrix");
      Bab->print("Alpha-beta B matrix");
    }

    //
    // Compute basis set completeness
    //
    double traceV_aa = Vaa->trace();
    double traceB_aa = Baa->trace();
    double traceV_ab = Vab->trace();
    double traceB_ab = Bab->trace();

    ExEnv::out0() << endl;
    ExEnv::out0() << indent << "Basis Set completeness diagnostics:" << endl;
    ExEnv::out0() << indent
		 << "-Tr(V)/Tr(B) for alpha-alpha pairs:" << indent <<
      scprintf("%10.6lf",(-1.0)*traceV_aa/traceB_aa) << endl;
    ExEnv::out0() << indent
		 << "-Tr(V)/Tr(B) for alpha-beta pairs:" << indent <<
      scprintf("%10.6lf",(-1.0)*traceV_ab/traceB_ab) << endl;

    //
    // Evaluate pair energies
    //
    epair_aa = local_matrix_kit->vector(npair_aa);
    epair_ab = local_matrix_kit->vector(npair_ab);

    // For some reason invert_this doesn't work here
    Baa->gen_invert_this();
    Bab->gen_invert_this();
    if (debug_ > 1) {
      Baa->print("Inverse alpha-alpha B matrix");
      Bab->print("Inverse alpha-beta B matrix");
    }
  
    double er12a_aa = 0.0;
    double er12a_ab = 0.0;
    
    ExEnv::out0() << endl << indent << "Alpha-alpha MBPT2-R12/A pair energies:" << endl;
    ExEnv::out0() << indent << scprintf("    i       j         e(ij)") << endl;
    ExEnv::out0() << indent << scprintf("  -----   -----   ------------") << endl;
    int ij, kl;
    for(i=0,ij=0;i<nocc_act;i++)
      for(j=0;j<i;j++,ij++) {
	double eij = 0.0;
	for(kl=0;kl<naa;kl++)
	  for(int mn=0;mn<naa;mn++) {
	    eij -= 2.0*Vaa->get_element(kl,ij)*Vaa->get_element(mn,ij)*Baa->get_element(kl,mn);
	  }
	epair_aa->set_element(ij,eij);
	er12a_aa += eij;
	ExEnv::out0() << indent << scprintf("  %3d     %3d     %12.9lf",i+1,j+1,eij) << endl;
      }
    
    ExEnv::out0() << endl << indent << "Alpha-beta MBPT2-R12/A pair energies:" << endl;
    ExEnv::out0() << indent << scprintf("    i       j         e(ij)") << endl;
    ExEnv::out0() << indent << scprintf("  -----   -----   ------------") << endl;
    for(i=0,ij=0;i<nocc_act;i++)
      for(j=0;j<nocc_act;j++,ij++) {
	double eij = 0.0;
	for(kl=0;kl<nab;kl++)
	  for(int mn=0;mn<nab;mn++) {
	    eij -= 1.0*Vab->get_element(kl,ij)*Vab->get_element(mn,ij)*Bab->get_element(kl,mn);
	  }
	epair_ab->set_element(ij,eij);
	er12a_ab += eij;
	ExEnv::out0() << indent << scprintf("  %3d     %3d     %12.9lf",i+1,j+1,eij) << endl;
      }

    er12a = er12a_aa + er12a_ab;
    r12_corr_energy_ = er12a;
  }
  tim_exit("mp2-r12a energy");

  ///////////////////////////////////////////////////////////////
  // The computation of the MP2 energy is now complete on each
  // node;
  ///////////////////////////////////////////////////////////////
  msg_->sum(aoint_computed);

#if PRINT_BIGGEST_INTS
  biggest_ints_1.combine(msg_);
  biggest_ints_2.combine(msg_);
  biggest_ints_2s.combine(msg_);
  biggest_ints_3a.combine(msg_);
  biggest_ints_3.combine(msg_);
#endif

  if (me == 0) {
    emp2r12a = escf + emp2 + er12a;
    
#if PRINT_BIGGEST_INTS
    ExEnv::out0() << "biggest 1/4 transformed ints" << endl;
    for (i=0; i<biggest_ints_1.ncontrib(); i++) {
      ExEnv::outn() << scprintf("%3d %3d %3d %3d %16.12f",
			       biggest_ints_1.indices(i)[0],
			       biggest_ints_1.indices(i)[1],
			       biggest_ints_1.indices(i)[2],
			       biggest_ints_1.indices(i)[3],
			       biggest_ints_1.val(i)
			       )
		   << endl;
    }
    ExEnv::out0() << "biggest 2/4 transformed ints" << endl;
    for (i=0; i<biggest_ints_2.ncontrib(); i++) {
      ExEnv::outn() << scprintf("%3d %3d %3d %3d %16.12f",
			       biggest_ints_2.indices(i)[0],
			       biggest_ints_2.indices(i)[1],
			       biggest_ints_2.indices(i)[2],
			       biggest_ints_2.indices(i)[3],
			       biggest_ints_2.val(i)
			       )
		   << endl;
    }
    ExEnv::out0() << "restricted 2/4 transformed ints" << endl;
    for (i=0; i<biggest_ints_2s.ncontrib(); i++) {
      ExEnv::outn() << scprintf("%3d %3d %3d %3d %16.12f",
			       biggest_ints_2s.indices(i)[0],
			       biggest_ints_2s.indices(i)[1],
			       biggest_ints_2s.indices(i)[2],
			       biggest_ints_2s.indices(i)[3],
			       biggest_ints_2s.val(i)
			       )
		   << endl;
    }
    ExEnv::out0() << "biggest 3/4 transformed ints (in 3.)" << endl;
    for (i=0; i<biggest_ints_3a.ncontrib(); i++) {
      ExEnv::outn() << scprintf("%3d %3d %3d %3d %16.12f",
			       biggest_ints_3a.indices(i)[0],
			       biggest_ints_3a.indices(i)[1],
			       biggest_ints_3a.indices(i)[2],
			       biggest_ints_3a.indices(i)[3],
			       biggest_ints_3a.val(i)
			       )
		   << endl;
    }
    ExEnv::out0() << "biggest 3/4 transformed ints (in 4.)" << endl;
    for (i=0; i<biggest_ints_3.ncontrib(); i++) {
      ExEnv::outn() << scprintf("%3d %3d %3d %3d %16.12f",
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
    
    ExEnv::out0()<<endl<<indent
		 <<scprintf("RHF energy [au]:                           %17.12lf\n", escf);
    ExEnv::out0()<<indent
		 <<scprintf("MP2 correlation energy [au]:               %17.12lf\n", emp2);
    ExEnv::out0()<<indent
		 <<scprintf("(MBPT2)-R12/A correlation energy [au]:     %17.12lf\n", er12a);
    ExEnv::out0()<<indent
		<<scprintf("MBPT2-R12/A energy [au]:                   %17.12lf\n", emp2r12a);
    ExEnv::out0().flush();
  }
  set_energy(emp2r12a);
  set_actual_value_accuracy(reference_->actual_value_accuracy()
                            *ref_to_mp2_acc);

  /*--------------------------
    Cleanup
   --------------------------*/
  delete[] Vaa_ijkl;
  delete[] Taa_ijkl;
  delete[] Vab_ijkl;
  delete[] Tab_ijkl;
  
  for (i=0; i<thr_->nthread(); i++) {
    delete e12thread[i];
  }
  delete[] e12thread;
  
  
  delete[] tbints_; tbints_ = 0;
  delete[] scf_vector;
  delete[] scf_vector_dat;
  delete[] evals;
  tim_exit("r12a-mem");
  return;
}


////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
