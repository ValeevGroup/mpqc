//
// r12a_intermed_spinorb_abs.cc
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
#include <chemistry/qc/basis/obint.h>
#include <chemistry/qc/basis/symmint.h>
#include <chemistry/qc/mbpt/util.h>
#include <chemistry/qc/mbpt/bzerofast.h>
#include <chemistry/qc/mbptr12/mbptr12.h>
#include <chemistry/qc/mbptr12/trans123_r12a_abs.h>
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

#define PRINT3Q 0
#define PRINT4Q 0
#define PRINT4Q_MP2 0
#define PRINT_NUM_TE_TYPES 4
#define PRINT_R12_INTERMED 0

#define LINDEP_TOL 1.e-6

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
MBPT2_R12::compute_r12a_intermed_spinorb_ABS_()
{
  // Basis Sets
  Ref<GaussianBasisSet> bs = basis();
  Ref<GaussianBasisSet> bs_aux = aux_basis();
  integral()->set_basis(bs,bs,bs,bs_aux);

  // log2 of the erep tolerance
  // (erep < 2^tol => discard)
  const int tol = (int) (-10.0/log10(2.0));  // discard ints smaller than 10^-20

  int nij;        // number of i,j pairs on a node (for e.g., mo_int)
  const int num_te_types = 4;
  enum te_types {eri=0, r12=1, r12t1=2, r12t2=3};
  double *mo_int; // MO integrals of type (oo|og')
  double *integral_ijrs; // half-transformed integrals

  int nbasis = bs->nbasis();
  int nbasis_aux = bs_aux->nbasis();
  int nfuncmax = bs->max_nfunction_in_shell();
  int nfuncmax_aux = bs_aux->max_nfunction_in_shell();
  int nshell = bs->nshell();
  int nshell_aux = bs_aux->nshell();

  int nocc_act;
  int npass, me, nproc, rest;
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

  double *integral_ijsk;  // all three-quarter transformed two-el integrals
  double *integral_ijkx;  // mo integrals (y = any MO)

  BiggestContribs biggest_coefs(5,10);
  CharacterTable ct = molecule()->point_group()->char_table();

#if PRINT_BIGGEST_INTS
  BiggestContribs biggest_ints_2(4,40);
  BiggestContribs biggest_ints_2s(4,40);
  BiggestContribs biggest_ints_3a(4,40);
  BiggestContribs biggest_ints_3(4,40);
#endif

  tim_enter("r12a-abs-mem");

  me = msg_->me();
  nproc = msg_->n();

  ExEnv::out0() << endl << indent
	       << "Entered memgrp based MP2-R12/A routine" << endl;
  ExEnv::out0() << indent << scprintf("nproc = %i", nproc) << endl;

  nocc = 0;
  for (int i=0; i<oso_dimension()->n(); i++) {
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
    mem_static += nthreads * integral()->storage_required(&Integral::grt,bs,bs,bs,bs_aux); // integral evaluators
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
  mem_remaining += thr_->nthread() * integral()->storage_required(&Integral::grt,bs,bs,bs,bs_aux);
  
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
  int index = 0;
  for (int i=0; i<ni; i++) {
    for (int j=0; j<nocc_act; j++) {
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
  for (int i=0; i<nbasis; i++) {
    scf_vector[i] = &scf_vector_dat[i*noso];
    }
  for (int i=0; i<noso; i++) {
      evals[i] = evalmat(i);
    }

  Scf_Vec = 0;
  evalmat = 0;

  if (debug_ > 2 && me == 0) {
    for (int j=0; j<noso; j++) {
      ExEnv::out0() << indent
           << scprintf("eigenvalue[%3d] = %15.10lf", j, evals[j]);
      if (j < nfzc) ExEnv::out0() << " (frozen docc)";
      else if (j < nocc_act + nfzc) ExEnv::out0() << " (active docc)";
      else if (j < nvir + nocc_act + nfzc) ExEnv::out0() << " (active uocc)";
      ExEnv::out0() << endl;
    }
  }


  //////////////////////////////////////////////////
  // For the auxiliary basis form an orthogonalizer
  // and distribute to nodes
  // Use canonical orthogonalization
  // ( does Wavefunction call that symmetric? )
  ////////////////////////////////////////////////

  RefSCMatrix orthog_aux;
  RefSCDimension osodim_aux;
  {
  //
  // Compute overlap matrix
  //

  // Make an Integral and initialize with bs_aux
  Ref<Integral> integral_aux = integral()->clone();
  integral_aux->set_basis(bs_aux);
  Ref<PetiteList> pl_aux = integral_aux->petite_list();
  Ref<OneBodyInt> ov_aux_engine = integral_aux->overlap();

  // form skeleton s matrix
  RefSymmSCMatrix s(bs_aux->basisdim(), bs_aux->matrixkit());
  Ref<SCElementOp> ov =
    new OneBodyIntOp(new SymmOneBodyIntIter(ov_aux_engine, pl_aux));
  
  s.assign(0.0);
  s.element_op(ov);
  ov=0;
  if (debug_ > 1) s.print("AO skeleton overlap (auxiliary basis)");
  
  // then symmetrize it
  RefSCDimension sodim_aux = pl_aux->SO_basisdim();
  RefSymmSCMatrix overlap_aux(sodim_aux, bs_aux->so_matrixkit());
  pl_aux->symmetrize(s,overlap_aux);

  // and clean up a bit
  ov_aux_engine = 0;
  s = 0;
  integral_aux = 0;

  //
  // Compute canonical orthogonalizer
  //

  RefSCMatrix overlap_aux_eigvec;
  RefDiagSCMatrix overlap_isqrt_eigval;
  RefDiagSCMatrix overlap_sqrt_eigval;

  // diagonalize a copy of overlap_
  RefSymmSCMatrix M = overlap_aux.copy();
  RefSCMatrix U(M.dim(), M.dim(), M.kit());
  RefDiagSCMatrix m(M.dim(), M.kit());
  M.diagonalize(m,U);
  M = 0;
  Ref<SCElementMaxAbs> maxabsop = new SCElementMaxAbs;
  m.element_op(maxabsop.pointer());
  double maxabs = maxabsop->result();
  double s_tol = LINDEP_TOL * maxabs;

  double minabs = maxabs;
  BlockedDiagSCMatrix *bm = dynamic_cast<BlockedDiagSCMatrix*>(m.pointer());
  if (bm == 0) {
    ExEnv::out0() << "R12A_intermed_spinorb_abs: orthog_aux: expected blocked overlap" << endl;
  }
  
  double *pm_sqrt = new double[bm->dim()->n()];
  double *pm_isqrt = new double[bm->dim()->n()];
  int *pm_index = new int[bm->dim()->n()];
  int *nfunc = new int[bm->nblocks()];
  int nfunctot = 0;
  int nlindep = 0;
  for (int i=0; i<bm->nblocks(); i++) {
      nfunc[i] = 0;
      if (bm->block(i).null()) continue;
      int n = bm->block(i)->dim()->n();
      double *pm = new double[n];
      bm->block(i)->convert(pm);
      for (int j=0; j<n; j++) {
          if (pm[j] > s_tol) {
	    if (pm[j] < minabs) { minabs = pm[j]; }
	    pm_sqrt[nfunctot] = sqrt(pm[j]);
	    pm_isqrt[nfunctot] = 1.0/pm_sqrt[nfunctot];
	    pm_index[nfunctot] = j;
	    nfunc[i]++;
	    nfunctot++;
	  }
          else {
	    nlindep++;
	  }
      }
      delete[] pm;
  }

  if (debug_)
    ExEnv::out0() << indent << "Removed " << nlindep << " linearly dependent functions from the auxiliary basis" << endl;

  // make sure all nodes end up with exactly the same data
  MessageGrp::get_default_messagegrp()->bcast(nfunctot);
  MessageGrp::get_default_messagegrp()->bcast(nfunc, bm->nblocks());
  MessageGrp::get_default_messagegrp()->bcast(pm_sqrt,nfunctot);
  MessageGrp::get_default_messagegrp()->bcast(pm_isqrt,nfunctot);
  MessageGrp::get_default_messagegrp()->bcast(pm_index,nfunctot);

  osodim_aux = new SCDimension(nfunctot, bm->nblocks(),
			       nfunc, "ortho aux SO (canonical)");
  for (int i=0; i<bm->nblocks(); i++) {
    osodim_aux->blocks()->set_subdim(i, new SCDimension(nfunc[i]));
  }

  overlap_aux_eigvec = bs_aux->so_matrixkit()->matrix(sodim_aux, osodim_aux);
  BlockedSCMatrix *bev
    = dynamic_cast<BlockedSCMatrix*>(overlap_aux_eigvec.pointer());
  BlockedSCMatrix *bU
    = dynamic_cast<BlockedSCMatrix*>(U.pointer());
  int ifunc = 0;
  for (int i=0; i<bev->nblocks(); i++) {
    if (bev->block(i).null()) continue;
    for (int j=0; j<nfunc[i]; j++) {
      bev->block(i)->assign_column(
				   bU->block(i)->get_column(pm_index[ifunc]),j
				   );
      ifunc++;
    }
  }

  overlap_sqrt_eigval = bs_aux->so_matrixkit()->diagmatrix(osodim_aux);
  overlap_sqrt_eigval->assign(pm_sqrt);
  overlap_isqrt_eigval = bs_aux->so_matrixkit()->diagmatrix(osodim_aux);
  overlap_isqrt_eigval->assign(pm_isqrt);

  delete[] nfunc;
  delete[] pm_sqrt;
  delete[] pm_isqrt;
  delete[] pm_index;
  
  if (debug_ > 1) {
    overlap_aux.print("S aux");
    overlap_aux_eigvec.print("S aux eigvec");
    overlap_isqrt_eigval.print("s^(-1/2) aux eigval");
    overlap_sqrt_eigval.print("s^(1/2) aux eigval");
  }

  RefSCMatrix orthog_aux_so = overlap_isqrt_eigval * overlap_aux_eigvec.t();
  orthog_aux_so = orthog_aux_so.t();
  orthog_aux = pl_aux->evecs_to_AO_basis(orthog_aux_so);
  orthog_aux_so = 0;

  if (debug_ > 1)
    orthog_aux.print("Aux orthogonalizer");
  RefSCMatrix tmp = orthog_aux.t() * pl_aux->to_AO_basis(overlap_aux) * orthog_aux;
  if (debug_ > 1)
    tmp.print("Xt * S * X (Aux)");
  }

  int noso_aux = orthog_aux.coldim().n();
  double *orthog_aux_vector = new double[nbasis_aux*noso_aux];
  orthog_aux.convert(orthog_aux_vector);
  orthog_aux = 0;
  int *symorb_irrep_aux = new int[noso_aux];
  int orbnum=0;
  for (int i=0; i<osodim_aux->blocks()->nblock(); i++) {
    for (int j=0; j<osodim_aux->blocks()->size(i); j++, orbnum++) {
      symorb_irrep_aux[orbnum] = i;
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
  for (int i=0; i<thr_->nthread(); i++) {
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
  mem->set_localsize(num_te_types*nijmax*nocc*nbasis_aux*sizeof(double));
  ExEnv::out0() << indent
	       << "Size of global distributed array:       "
	       << mem->totalsize()
	       << " Bytes"
	       << endl;
  MemoryGrpBuf<double> membuf_remote(mem);


  Ref<ThreadLock> lock = thr_->new_lock();
  R12A_ABS_123Qtr** e123thread = new R12A_ABS_123Qtr*[thr_->nthread()];
  for (int i=0; i<thr_->nthread(); i++) {
    e123thread[i] = new R12A_ABS_123Qtr(i, thr_->nthread(), me, nproc,
					mem, msg_, lock, bs, bs_aux, tbints_[i],
					nocc, nocc_act, scf_vector, tol, debug_,
					dynamic_);
  }

  if (npass > 1 || restart_orbital_r12_memgrp_) {
    bool restart = restart_orbital_r12_memgrp_;
    // using File integrals accumulator when npass > 1 or restarting
#if HAVE_MPIIO
    ExEnv::out0() << indent << "Will use MPI-IO (individual I/O) to handle transformed integrals" << endl;
    r12intsacc_ = new R12IntsAcc_MPIIOFile_Ind(mem,r12ints_file_,num_te_types,nocc,nbasis_aux,nocc,nfzc,restart);
#else
    ExEnv::out0() << indent << "Will use POSIX I/O on node 0 to handle transformed integrals" << endl;
    r12intsacc_ = new R12IntsAcc_Node0File(mem,r12ints_file_,num_te_types,nocc,nbasis_aux,nocc,nfzc,restart);
#endif
  }
  else {
    // using MemoryGrp integrals accumulator when npass = 1 and not restarting
//    ExEnv::out0() << indent << "Will use MPI-IO (individual I/O) to handle transformed integrals" << endl;
//    bool restart = restart_orbital_r12_memgrp_;
//    r12intsacc_ = new R12IntsAcc_MPIIOFile_Ind(mem,r12ints_file_,num_te_types,nocc,nbasis_aux,nocc,nfzc,restart);
    ExEnv::out0() << indent << "Will hold transformed integrals in memory" << endl;
    r12intsacc_ = new R12IntsAcc_MemoryGrp(mem,num_te_types,nocc,nbasis_aux,nocc,nfzc);
  }


  /*-----------------------------------

    Start the integrals transformation

   -----------------------------------*/
  tim_enter("mp2-r12/a passes");
  if (me == 0 && if_to_checkpoint()) {
    StateOutBin stateout(checkpoint_file());
    SavableState::save_state(this,stateout);
  }
  for (int pass=0; pass<npass; pass++) {

    ExEnv::out0() << indent << "Beginning pass " << pass+1 << endl;

    int i_offset = restart_orbital_r12_memgrp_ + pass*ni + nfzc;
    if ((pass == npass - 1) && (rest != 0)) ni = rest;

    // Compute number of of i,j pairs on each node during current pass for
    // two-el integrals
    index = 0;
    nij = 0;
    for (int i=0; i<ni; i++) {
      for (int j=0; j<nocc_act; j++) {
        if (index++ % nproc == me) nij++;
      }
    }

    if (debug_)
      ExEnv::outn() << indent << "node " << me << ", nij = " << nij << endl;

    // Allocate and initialize some arrays
    // (done here to avoid having these arrays
    // overlap with arrays allocated later)

    // Allocate (and initialize) some arrays
    integral_ijsk = (double*) mem->localdata();
    bzerofast(integral_ijsk, (num_te_types*nij*nocc*nbasis_aux));
    integral_ijsk = 0;
    mem->sync();
    ExEnv::out0() << indent
		   << scprintf("Begin loop over shells (grt, 1.+2. q.t.)") << endl;

    // Do the two electron integrals and the first two quarter transformations
    tim_enter("grt+1.qt+2.qt");
    for (int i=0; i<thr_->nthread(); i++) {
      e123thread[i]->set_i_offset(i_offset);
      e123thread[i]->set_ni(ni);
      thr_->add_thread(i,e123thread[i]);
#     if SINGLE_THREAD_E12
      e123thread[i]->run();
#     endif
      }
#   if !SINGLE_THREAD_E12
    thr_->start_threads();
    thr_->wait_threads();
#   endif
    tim_exit("grt+1.qt+2.qt");
    ExEnv::out0() << indent << "End of loop over shells" << endl;

    mem->sync();  // Make sure ijsk is complete on each node before continuing
    integral_ijsk = (double*) mem->localdata();

#if PRINT3Q
    if (me == 0) {
      int index = 0;
      int ij_index = 0;
      int j_offset = nfzc;
      for (int i = 0; i<ni; i++) {
	for (int j = 0; j<nocc_act; j++) {
	  if (index++ % nproc == me) {
	    double *ij_offset = integral_ijsk + num_te_types*nocc*nbasis_aux*ij_index;
	    for(int te_type=0; te_type<PRINT_NUM_TE_TYPES; ++te_type, ij_offset+=nocc*nbasis_aux) {
	      double *ijsk_ptr = ij_offset;
	      for (int s = 0; s<nbasis_aux; s++) {
		for (int k = 0; k<nocc; k++) {
		  printf("3Q: type = %d (%d %d|%d %d) = %12.8f\n",
		       te_type,i+i_offset,k,j+j_offset,s,*ijsk_ptr);
		ijsk_ptr++;
		}
	      }
	    }
	    ij_index++;
	  }
	}
      }
    }
#endif

    double *kx_tmp = new double[nocc*noso_aux];
    // in ikjx (stored as ijkx): i act; j act; k occ; x any MO.

    // Begin fourth quarter transformation
    ExEnv::out0() << indent << "Begin fourth q.t." << endl;
    tim_enter("4. q.t.");
    index = 0;
    ij_index = 0;
    for (int i=0; i<ni; i++) {
      for (int j=0; j<nocc_act; j++) {
        if (index++ % nproc == me) {
	  double *ij_offset = integral_ijsk + num_te_types*nocc*nbasis_aux*ij_index;

	  for(int te_type=0; te_type<num_te_types; te_type++,ij_offset+=nocc*nbasis_aux) {
	    bzerofast(kx_tmp, nocc*noso_aux);
	    double *ijs_offset = ij_offset;

	    int sx = 0;
	    for (int s=0; s<nbasis_aux; s++, ijs_offset+=nocc) {

	      for (int x=0; x<noso_aux; x++, ++sx) {
		double c_sx = orthog_aux_vector[sx];
		double *kx_ptr = kx_tmp + x;
		double *sk_ptr = ijs_offset;

		for (int k=0; k<nocc; ++k) {
		  *kx_ptr += c_sx * *sk_ptr;
		  ++sk_ptr;
		  kx_ptr += noso_aux;
		} // end of k loop
	      } // end of x loop
	    } // end of s loop

	    // Put kx_tmp into ijsk for one i,j while
	    // overwriting elements of ijsk
	    double *ijsk_ptr = ij_offset;
	    double *ijkx_ptr = kx_tmp;
	    int nkx = nocc*noso_aux;
	    for (int kx=0; kx<nkx; ++kx) {
	      *ijsk_ptr = *ijkx_ptr;
	      ++ijsk_ptr;
	      ++ijkx_ptr;
	    }
	  } // end of te_type loop
          ij_index++;
	}   // endif
      }     // exit j loop
    }       // exit i loop
    // end of fourth quarter transformation
    tim_exit("4. q.t.");
    ExEnv::out0() << indent << "End of fourth q.t." << endl;

    // The array integral_ijsk has now been overwritten by MO integrals ijkx
    // rename the array mo_int
    mo_int = integral_ijsk;
    delete[] kx_tmp;

    // Zero out nonsymmetric integrals
    {
    int index = 0;
    int ij_index = 0;
    int j_offset = nfzc;
    for (int i = 0; i<ni; i++) {
      for (int j = 0; j<nocc_act; j++) {
	if (index++ % nproc == me) {
	  double *ij_offset = mo_int + num_te_types*nocc*nbasis_aux*ij_index;
	  for(int te_type=0; te_type<num_te_types; ++te_type,ij_offset+=nocc*nbasis_aux) {
	    double *ijkx_ptr = ij_offset;
	    for (int k=0; k<nocc; ++k) {
	      for (int x=0; x<noso_aux; ++x) {
		if (( symorb_irrep_[i+i_offset] ^
		       symorb_irrep_[j+j_offset] ^
		       symorb_irrep_[k] ^
		       symorb_irrep_aux[x]) ) {
		  *ijkx_ptr = 0.0;
		}
		ijkx_ptr++;
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
	    double *ij_offset = mo_int + num_te_types*nocc*nbasis_aux*ij_index;
	    for(int te_type=0; te_type<PRINT_NUM_TE_TYPES; te_type++,ij_offset+=nocc*nbasis_aux) {
	      double *ijkx_ptr = ij_offset;
	      for (int k=0; k<nocc; ++k) {
		for (int x=0; x<noso_aux; ++x) {
		  printf("4Q: type = %d (%d %d|%d %d) = %12.8f\n",
		       te_type,i+i_offset,k,j+j_offset,x,*ijkx_ptr);
		  ijkx_ptr++;
		}
	      }
	    }
	    ij_index++;
	  }
	}
      }
    }
#endif

    integral_ijsk = 0;
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
  double *Xaa_ijkl = new double[naa*naa];
  double *Taa_ijkl = new double[naa*naa];
  double *Vab_ijkl = new double[nab*nab];
  double *Xab_ijkl = new double[nab*nab];
  double *Tab_ijkl = new double[nab*nab];
  if (debug_)
    ExEnv::out0() << indent << "Allocated intermediates V, X, and T" << endl;
  bzerofast(Vaa_ijkl,naa*naa);
  bzerofast(Xaa_ijkl,naa*naa);
  bzerofast(Taa_ijkl,naa*naa);
  bzerofast(Vab_ijkl,nab*nab);
  bzerofast(Xab_ijkl,nab*nab);
  bzerofast(Tab_ijkl,nab*nab);

  // Compute intermediates
  if (debug_)
    ExEnv::out0() << indent << "Ready to compute intermediates V, X, and T" << endl;
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
        double *klox_buf_eri = r12intsacc_->retrieve_pair_block(k,l,R12IntsAcc::eri);
        double *klox_buf_r12 = r12intsacc_->retrieve_pair_block(k,l,R12IntsAcc::r12);
        double *klox_buf_r12t1 = r12intsacc_->retrieve_pair_block(k,l,R12IntsAcc::r12t1);
        double *klox_buf_r12t2 = r12intsacc_->retrieve_pair_block(k,l,R12IntsAcc::r12t2);
        
	double *lkox_buf_eri = r12intsacc_->retrieve_pair_block(l,k,R12IntsAcc::eri);
	double *lkox_buf_r12 = r12intsacc_->retrieve_pair_block(l,k,R12IntsAcc::r12);
	double *lkox_buf_r12t1 = r12intsacc_->retrieve_pair_block(l,k,R12IntsAcc::r12t1);
	double *lkox_buf_r12t2 = r12intsacc_->retrieve_pair_block(l,k,R12IntsAcc::r12t2);
        tim_exit("MO ints retrieve");

	int ij = 0;
        int ij_aa=-1;
        for(int i=0;i<nocc_act;i++)
          for(int j=0;j<=i;j++,ij++) {

            double pfac_ij = (i==j) ? oosqrt2 : 1.0;
	    if (i != j)
	      ++ij_aa;
	    int ij_ab = i*nocc_act + j;
	    int ji_ab = j*nocc_act + i;
            
            tim_enter("MO ints retrieve");
            double *ijox_buf_r12 = r12intsacc_->retrieve_pair_block(i,j,R12IntsAcc::r12);
	    double *jiox_buf_r12 = r12intsacc_->retrieve_pair_block(j,i,R12IntsAcc::r12);
            tim_exit("MO ints retrieve");

	    double *Vaa_ij = Vaa_ijkl + ij_aa*naa;
            double *Vab_ij = Vab_ijkl + ij_ab*nab;
            double *Vab_ji = Vab_ijkl + ji_ab*nab;
	    double *Xaa_ij = Xaa_ijkl + ij_aa*naa;
            double *Xab_ij = Xab_ijkl + ij_ab*nab;
            double *Xab_ji = Xab_ijkl + ji_ab*nab;
            double *Taa_ij = Taa_ijkl + ij_aa*naa;
            double *Tab_ij = Tab_ijkl + ij_ab*nab;
            double *Tab_ji = Tab_ijkl + ji_ab*nab;
            
            tim_enter("MO ints contraction");
            double Vaa_ijkl, Vab_ijkl, Vab_jikl, Vab_ijlk, Vab_jilk;
            double Xaa_ijkl, Xab_ijkl, Xab_jikl, Xab_ijlk, Xab_jilk;
	    double Taa_ijkl, Tab_ijkl, Tab_jikl, Tab_ijlk, Tab_jilk;
	    Vaa_ijkl = Vab_ijkl = Vab_jikl = Vab_ijlk = Vab_jilk = 0.0;
	    Xaa_ijkl = Xab_ijkl = Xab_jikl = Xab_ijlk = Xab_jilk = 0.0;
	    Taa_ijkl = Tab_ijkl = Tab_jikl = Tab_ijlk = Tab_jilk = 0.0;
            for(int o=0;o<nocc;o++) {
	      double pfac_xy = 1.0;
              for(int x=0;x<noso_aux;x++) {
                int ox_offset = o*nbasis_aux+x;
                double ij_r12_ox = ijox_buf_r12[ox_offset];
                double ji_r12_ox = jiox_buf_r12[ox_offset];
                double kl_eri_ox = klox_buf_eri[ox_offset];
                double lk_eri_ox = lkox_buf_eri[ox_offset];
                Vab_ijkl -= pfac_xy * (ij_r12_ox * kl_eri_ox + ji_r12_ox * lk_eri_ox);
		if (i != j)
		  Vab_jikl -= pfac_xy * (ji_r12_ox * kl_eri_ox + ij_r12_ox * lk_eri_ox);
		if (k != l)
		  Vab_ijlk -= pfac_xy * (ij_r12_ox * lk_eri_ox + ji_r12_ox * kl_eri_ox);
		if (i != j && k != l) {
		  Vaa_ijkl -= pfac_xy * (ij_r12_ox - ji_r12_ox)*(kl_eri_ox - lk_eri_ox);
		  Vab_jilk -= pfac_xy * (ij_r12_ox * kl_eri_ox + ji_r12_ox * lk_eri_ox);
		}
                double kl_r12_ox = klox_buf_r12[ox_offset];
                double lk_r12_ox = lkox_buf_r12[ox_offset];
                Xab_ijkl -= pfac_xy * (ij_r12_ox * kl_r12_ox + ji_r12_ox * lk_r12_ox);
		if (i != j)
		  Xab_jikl -= pfac_xy * (ji_r12_ox * kl_r12_ox + ij_r12_ox * lk_r12_ox);
		if (k != l)
		  Xab_ijlk -= pfac_xy * (ij_r12_ox * lk_r12_ox + ji_r12_ox * kl_r12_ox);
		if (i != j && k != l) {
		  Xaa_ijkl -= pfac_xy * (ij_r12_ox - ji_r12_ox)*(kl_r12_ox - lk_r12_ox);
		  Xab_jilk -= pfac_xy * (ij_r12_ox * kl_r12_ox + ji_r12_ox * lk_r12_ox);
		}
                double kl_r12t1_ox = klox_buf_r12t1[ox_offset];
                double kl_r12t2_ox = klox_buf_r12t2[ox_offset];
                double lk_r12t1_ox = lkox_buf_r12t1[ox_offset];
                double lk_r12t2_ox = lkox_buf_r12t2[ox_offset];
                double kl_Tr12_ox = -kl_r12t1_ox-kl_r12t2_ox;
                double lk_Tr12_ox = -lk_r12t1_ox-lk_r12t2_ox;
                Tab_ijkl += pfac_xy * (ij_r12_ox * kl_Tr12_ox + ji_r12_ox * lk_Tr12_ox);
		if (i != j)
		  Tab_jikl += pfac_xy * (ji_r12_ox * kl_Tr12_ox + ij_r12_ox * lk_Tr12_ox);
		if (k != l)
		  Tab_ijlk += pfac_xy * (ij_r12_ox * lk_Tr12_ox + ji_r12_ox * kl_Tr12_ox);
		if (i != j && k != l) {
		  Taa_ijkl += pfac_xy * (ij_r12_ox - ji_r12_ox)*(kl_Tr12_ox - lk_Tr12_ox);
		  Tab_jilk += pfac_xy * (ij_r12_ox * kl_Tr12_ox + ji_r12_ox * lk_Tr12_ox);
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
	    r12intsacc_->release_pair_block(i,j,R12IntsAcc::r12);
	    r12intsacc_->release_pair_block(j,i,R12IntsAcc::r12);
          }
        r12intsacc_->release_pair_block(k,l,R12IntsAcc::eri);
        r12intsacc_->release_pair_block(k,l,R12IntsAcc::r12);
        r12intsacc_->release_pair_block(k,l,R12IntsAcc::r12t1);
        r12intsacc_->release_pair_block(k,l,R12IntsAcc::r12t2);
	r12intsacc_->release_pair_block(l,k,R12IntsAcc::eri);
	r12intsacc_->release_pair_block(l,k,R12IntsAcc::r12);
	r12intsacc_->release_pair_block(l,k,R12IntsAcc::r12t1);
	r12intsacc_->release_pair_block(l,k,R12IntsAcc::r12t2);
      }
  }
  delete[] proc_with_ints;
  tim_exit("mp2-r12a intermeds");
  r12intsacc_->deactivate();
  if (debug_)
    ExEnv::out0() << indent << "Computed intermediates V, X, and T" << endl;

  if (nproc > 1) {
    // Use MemoryGrp to send all contributions to intermediates V, X, and T to node 0
    msg_->sum(Vaa_ijkl,naa*naa,0,0);
    msg_->sum(Vab_ijkl,nab*nab,0,0);
    msg_->sum(Xaa_ijkl,naa*naa,0,0);
    msg_->sum(Xab_ijkl,nab*nab,0,0);
    msg_->sum(Taa_ijkl,naa*naa,0,0);
    msg_->sum(Tab_ijkl,nab*nab,0,0);
  }

  if (debug_)
    ExEnv::out0() << indent << "Gathered intermediates V, X, and T on node 0" << endl;

  // Add intermediates contribution to their global values
  if (me == 0) {
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
  }
  msg_->sum(aoint_computed);

#if PRINT_BIGGEST_INTS
  biggest_ints_1.combine(msg_);
  biggest_ints_2.combine(msg_);
  biggest_ints_2s.combine(msg_);
  biggest_ints_3a.combine(msg_);
  biggest_ints_3.combine(msg_);
#endif

  if (me == 0) {
#if PRINT_BIGGEST_INTS
    ExEnv::out0() << "biggest 1/4 transformed ints" << endl;
    for (int i=0; i<biggest_ints_1.ncontrib(); i++) {
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
    for (int i=0; i<biggest_ints_2.ncontrib(); i++) {
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
    for (int i=0; i<biggest_ints_2s.ncontrib(); i++) {
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
    for (int i=0; i<biggest_ints_3a.ncontrib(); i++) {
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
    for (int i=0; i<biggest_ints_3.ncontrib(); i++) {
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
  }

  /*--------------------------
    Cleanup
   --------------------------*/
  delete[] Vaa_ijkl;
  delete[] Vab_ijkl;
  delete[] Xaa_ijkl;
  delete[] Xab_ijkl;
  delete[] Taa_ijkl;
  delete[] Tab_ijkl;
  
  for (int i=0; i<thr_->nthread(); i++) {
    delete e123thread[i];
  }
  delete[] e123thread;
  
  
  delete[] tbints_; tbints_ = 0;
  delete[] scf_vector;
  delete[] scf_vector_dat;
  delete[] evals;
  tim_exit("r12a-abs-mem");
  return;
}


////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
