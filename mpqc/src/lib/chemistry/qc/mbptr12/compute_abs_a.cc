//
// compute_abs_a.cc
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
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/obint.h>
#include <chemistry/qc/basis/symmint.h>
#include <chemistry/qc/basis/orthog.h>
#include <chemistry/qc/mbpt/util.h>
#include <chemistry/qc/mbpt/bzerofast.h>
#include <chemistry/qc/mbptr12/trans123_r12a_abs.h>
#include <chemistry/qc/mbptr12/r12ia.h>
#include <chemistry/qc/mbptr12/r12ia_memgrp.h>
#include <chemistry/qc/mbptr12/r12ia_node0file.h>
#ifdef HAVE_MPIIO
  #include <chemistry/qc/mbptr12/r12ia_mpiiofile.h>
#endif
#include <chemistry/qc/mbptr12/vxb_eval_abs_a.h>

using namespace std;
using namespace sc;

#define SINGLE_THREAD_E123   0
#define PRINT3Q 0
#define PRINT4Q 0
#define PRINT4Q_MP2 0
#define PRINT_NUM_TE_TYPES 4
#define PRINT_R12_INTERMED 0
#define LINDEP_TOL 1.e-6

#define USE_GLOBAL_ORTHOG 1

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
R12IntEval_abs_A::compute(RefSCMatrix& Vaa, RefSCMatrix& Xaa, RefSCMatrix& Baa,
			  RefSCMatrix& Vab, RefSCMatrix& Xab, RefSCMatrix& Bab)
{
  int debug_ = r12info()->debug_level();

  MolecularEnergy* mole = r12info()->mole();
  Ref<Integral> integral = r12info()->integral();
  Ref<GaussianBasisSet> bs = r12info()->basis();
  Ref<GaussianBasisSet> bs_aux = r12info()->basis_aux();
  bool two_basis_form = (bs != bs_aux);
  if (!two_basis_form)
    throw std::runtime_error("R12IntEval_abs_A::compute called when basis sets are identical");
  integral->set_basis(bs,bs,bs,bs_aux);

  Ref<SCMatrixKit> matrixkit_aux = bs_aux->matrixkit();
  Ref<SCMatrixKit> so_matrixkit_aux = bs_aux->so_matrixkit();
  Ref<MessageGrp> msg = r12info()->msg();
  Ref<MemoryGrp> mem = r12info()->mem();
  Ref<ThreadGrp> thr = r12info()->thr();
  const int num_te_types = 4;
  enum te_types {eri=0, r12=1, r12t1=2, r12t2=3};

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

  tim_enter("r12a-abs-mem");

  int me = msg->me();
  int nproc = msg->n();
  const size_t mem_alloc = r12info()->memory();

  int nbasis = bs->nbasis();
  int nbasis_aux = bs_aux->nbasis();
  int nfuncmax = bs->max_nfunction_in_shell();
  int nfuncmax_aux = bs_aux->max_nfunction_in_shell();
  int nshell = bs->nshell();
  int nshell_aux = bs_aux->nshell();
  int nocc = r12info()->nocc();
  int nocc_act = r12info()->nocc_act();
  int nfzc = r12info()->nfzc();
  int nfzv = r12info()->nfzv();
  int noso = r12info()->noso();
  int nvir  = noso - nocc;

  ExEnv::out0() << endl << indent
	       << "Entered ABS A intermediates evaluator" << endl;
  ExEnv::out0() << indent << scprintf("nproc = %i", nproc) << endl;


  // Do a few preliminary tests to make sure the desired calculation
  // can be done (and appears to be meaningful!)

  if (nocc_act <= 0)
    throw std::runtime_error("There are no active occupied orbitals; program exiting");

  if (nvir-nfzv <= 0)
    throw std::runtime_error("There are no active virtual orbitals; program exiting");
    
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
    mem_static = nbasis*noso + nbasis_aux*nbasis_aux; // scf vector + aux_basis orthogonalizer
    mem_static *= sizeof(double);
    int nthreads = thr->nthread();
    mem_static += nthreads * integral->storage_required_grt(bs,bs,bs,bs_aux); // integral evaluators
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
  mem_remaining += thr->nthread() * integral->storage_required_grt(bs,bs,bs,bs_aux);
  
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
    throw std::runtime_error("R12IntEval_abs_A: batch size is 0: more memory or processors are needed");
  
  if (r12info()->dynamic())
    ExEnv::out0() << indent << "Using dynamic load balancing." << endl;

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
		<< scprintf(" npass  rest  nbasis  nshell  nfuncmax  nbasis(ABS) nshell(ABS) nfuncmax(ABS)") << endl;
  ExEnv::out0() << indent
		<< scprintf("  %-4i   %-3i   %-5i    %-4i     %-3i   %-5i        %-4i         %-3i",
			    npass,rest,nbasis,nshell,nfuncmax,nbasis_aux,nshell_aux,nfuncmax_aux)
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
    Ref<Integral> integral_aux = integral->clone();
    integral_aux->set_basis(bs_aux);
    Ref<PetiteList> pl_aux = integral_aux->petite_list();
    Ref<OneBodyInt> ov_aux_engine = integral_aux->overlap();

    // form skeleton s matrix
    RefSymmSCMatrix s(bs_aux->basisdim(), matrixkit_aux);
    Ref<SCElementOp> ov =
      new OneBodyIntOp(new SymmOneBodyIntIter(ov_aux_engine, pl_aux));
    
    s.assign(0.0);
    s.element_op(ov);
    ov=0;
    if (debug_ > 1) s.print("AO skeleton overlap (auxiliary basis)");
    
    // then symmetrize it
    RefSCDimension sodim_aux = pl_aux->SO_basisdim();
    RefSymmSCMatrix overlap_aux(sodim_aux, so_matrixkit_aux);
    pl_aux->symmetrize(s,overlap_aux);
    
    // and clean up a bit
    ov_aux_engine = 0;
    s = 0;
    integral_aux = 0;

    //
    // Compute orthogonalizer for the auxiliary basis
    //
#if USE_GLOBAL_ORTHOG
    OverlapOrthog orthog(OverlapOrthog::Canonical,
                         overlap_aux,
                         so_matrixkit_aux,
                         LINDEP_TOL,
                         0);

    RefSCMatrix orthog_aux_so = orthog.basis_to_orthog_basis();
    orthog_aux_so = orthog_aux_so.t();
    osodim_aux = orthog_aux_so.coldim();
#else    
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
    msg->bcast(nfunctot);
    msg->bcast(nfunc, bm->nblocks());
    msg->bcast(pm_sqrt,nfunctot);
    msg->bcast(pm_isqrt,nfunctot);
    msg->bcast(pm_index,nfunctot);
    osodim_aux = new SCDimension(nfunctot, bm->nblocks(),
				 nfunc, "ortho aux SO (canonical)");
    for (int i=0; i<bm->nblocks(); i++) {
      osodim_aux->blocks()->set_subdim(i, new SCDimension(nfunc[i]));
    }

    overlap_aux_eigvec = so_matrixkit_aux->matrix(sodim_aux, osodim_aux);
    BlockedSCMatrix *bev
      = dynamic_cast<BlockedSCMatrix*>(overlap_aux_eigvec.pointer());
    BlockedSCMatrix *bU
      = dynamic_cast<BlockedSCMatrix*>(U.pointer());

    int ifunc = 0;
    for (int i=0; i<bev->nblocks(); i++) {
      if (bev->block(i).null()) continue;
      RefSCMatrix bev_block = bev->block(i);
      RefSCMatrix bU_block = bU->block(i);
      for (int j=0; j<nfunc[i]; j++) {
// For some reason this produces a bunch of temp objects which are never destroyed
// and default_matrixkit is not detroyed at the end
//	bev->block(i)->assign_column(
//				     bU->block(i)->get_column(pm_index[ifunc]),j
//				     );
        int nk = bev_block->rowdim().n();
        for (int k=0; k<nk; k++)
          bev_block->set_element(k,j,bU_block->get_element(k,pm_index[ifunc]));
	ifunc++;
      }
    }

    overlap_sqrt_eigval = so_matrixkit_aux->diagmatrix(osodim_aux);
    overlap_sqrt_eigval->assign(pm_sqrt);
    overlap_isqrt_eigval = so_matrixkit_aux->diagmatrix(osodim_aux);
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
#endif

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
  integral->set_storage(mem_remaining);
  Ref<TwoBodyInt>* tbints = new Ref<TwoBodyInt>[thr->nthread()];
  for (int i=0; i<thr->nthread(); i++) {
    tbints[i] = integral->grt();
  }
  ExEnv::out0() << indent
		<< scprintf("Memory used for integral storage:       %i Bytes",
			    integral->storage_used())
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
  
  Ref<ThreadLock> lock = thr->new_lock();
  R12A_ABS_123Qtr** e123thread = new R12A_ABS_123Qtr*[thr->nthread()];
  for (int i=0; i<thr->nthread(); i++) {
    e123thread[i] = new R12A_ABS_123Qtr(i, thr->nthread(), me, nproc,
					mem, msg, lock, bs, bs_aux, tbints[i],
					nocc, nocc_act, scf_vector, tol, debug_,
					r12info()->dynamic());
  }

  Ref<R12IntsAcc> r12intsacc;
  if (npass > 1 || restart_orbital_) {
    bool restart = (restart_orbital_ > 0);
    const char *r12ints_file = r12info()->ints_file();
    // using File integrals accumulator when npass > 1 or restarting
#if HAVE_MPIIO
    ExEnv::out0() << indent << "Will use MPI-IO (individual I/O) to handle transformed integrals" << endl;
    r12intsacc = new R12IntsAcc_MPIIOFile_Ind(mem,r12ints_file,num_te_types,nocc,nbasis_aux,nocc,nfzc,restart);
#else
    ExEnv::out0() << indent << "Will use POSIX I/O on node 0 to handle transformed integrals" << endl;
    r12intsacc = new R12IntsAcc_Node0File(mem,r12ints_file,num_te_types,nocc,nbasis_aux,nocc,nfzc,restart);
#endif
  }
  else {
    // using MemoryGrp integrals accumulator when npass = 1 and not restarting
//    ExEnv::out0() << indent << "Will use MPI-IO (individual I/O) to handle transformed integrals" << endl;
//    bool restart = restart_orbital_;
//    const char *r12ints_file = r12info()->ints_file();
//    r12intsacc = new R12IntsAcc_MPIIOFile_Ind(mem,r12ints_file,num_te_types,nocc,nbasis_aux,nocc,nfzc,restart);
    ExEnv::out0() << indent << "Will hold transformed integrals in memory" << endl;
    r12intsacc = new R12IntsAcc_MemoryGrp(mem,num_te_types,nocc,nbasis_aux,nocc,nfzc);
  }


  /*-----------------------------------

    Start the integrals transformation

   -----------------------------------*/
  tim_enter("mp2-r12/a passes");
  if (me == 0 && mole->if_to_checkpoint() && r12intsacc->can_restart()) {
    StateOutBin stateout(mole->checkpoint_file());
    SavableState::save_state(mole,stateout);
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

    if (debug_)
      ExEnv::outn() << indent << "node " << me << ", nij = " << nij << endl;

    // Allocate and initialize some arrays
    // (done here to avoid having these arrays
    // overlap with arrays allocated later)
    
    // Allocate (and initialize) some arrays
    double* integral_ijsk = (double*) mem->localdata();
    bzerofast(integral_ijsk, (num_te_types*nij*nocc*nbasis_aux));
    integral_ijsk = 0;
    mem->sync();
    ExEnv::out0() << indent
		  << scprintf("Begin loop over shells (grt, 1.+2. q.t.)") << endl;

    // Do the two electron integrals and the first two quarter transformations
    tim_enter("grt+1.qt+2.qt");
    for (int i=0; i<thr->nthread(); i++) {
      e123thread[i]->set_i_offset(i_offset);
      e123thread[i]->set_ni(ni);
      thr->add_thread(i,e123thread[i]);
#     if SINGLE_THREAD_E123
      e123thread[i]->run();
#     endif
    }
#   if !SINGLE_THREAD_E123
    thr->start_threads();
    thr->wait_threads();
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
    int index = 0;
    int ij_index = 0;
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
    double* mo_int = integral_ijsk;
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
		if (( mo_irrep[i+i_offset] ^
		       mo_irrep[j+j_offset] ^
		       mo_irrep[k] ^
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
    r12intsacc->store_memorygrp(mem,ni);
    tim_exit("MO ints store");
    mem->sync();

    if (me == 0 && mole->if_to_checkpoint() && r12intsacc->can_restart()) {
      current_orbital_ += ni;
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
  r12intsacc->commit();

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

  if (r12intsacc->has_access(me)) {
    int kl = 0;
    for(int k=0;k<nocc_act;k++)
      for(int l=0;l<=k;l++,kl++) {
        double pfac_kl = (k==l) ? oosqrt2 : 1.0;
        int kl_proc = kl%nproc_with_ints;
        if (kl_proc != proc_with_ints[me])
          continue;
	int kl_aa = k*(k-1)/2 + l;
	int kl_ab = k*nocc_act + l;
	int lk_ab = l*nocc_act + k;
        
        // Get (|r12|) and (|[r12,T1]|) integrals only
        tim_enter("MO ints retrieve");
        double *klox_buf_eri = r12intsacc->retrieve_pair_block(k,l,R12IntsAcc::eri);
        double *klox_buf_r12 = r12intsacc->retrieve_pair_block(k,l,R12IntsAcc::r12);
        double *klox_buf_r12t1 = r12intsacc->retrieve_pair_block(k,l,R12IntsAcc::r12t1);
        double *klox_buf_r12t2 = r12intsacc->retrieve_pair_block(k,l,R12IntsAcc::r12t2);
        
	double *lkox_buf_eri = r12intsacc->retrieve_pair_block(l,k,R12IntsAcc::eri);
	double *lkox_buf_r12 = r12intsacc->retrieve_pair_block(l,k,R12IntsAcc::r12);
	double *lkox_buf_r12t1 = r12intsacc->retrieve_pair_block(l,k,R12IntsAcc::r12t1);
	double *lkox_buf_r12t2 = r12intsacc->retrieve_pair_block(l,k,R12IntsAcc::r12t2);
        tim_exit("MO ints retrieve");

	int ij = 0;
        for(int i=0;i<nocc_act;i++)
          for(int j=0;j<=i;j++,ij++) {

            double pfac_ij = (i==j) ? oosqrt2 : 1.0;
	    int ij_aa = i*(i-1)/2 + j;
	    int ij_ab = i*nocc_act + j;
	    int ji_ab = j*nocc_act + i;
            
            tim_enter("MO ints retrieve");
            double *ijox_buf_r12 = r12intsacc->retrieve_pair_block(i,j,R12IntsAcc::r12);
	    double *jiox_buf_r12 = r12intsacc->retrieve_pair_block(j,i,R12IntsAcc::r12);
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
                int ox_offset = o*noso_aux + x;
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
	    r12intsacc->release_pair_block(i,j,R12IntsAcc::r12);
	    r12intsacc->release_pair_block(j,i,R12IntsAcc::r12);
          }
        r12intsacc->release_pair_block(k,l,R12IntsAcc::eri);
        r12intsacc->release_pair_block(k,l,R12IntsAcc::r12);
        r12intsacc->release_pair_block(k,l,R12IntsAcc::r12t1);
        r12intsacc->release_pair_block(k,l,R12IntsAcc::r12t2);
	r12intsacc->release_pair_block(l,k,R12IntsAcc::eri);
	r12intsacc->release_pair_block(l,k,R12IntsAcc::r12);
	r12intsacc->release_pair_block(l,k,R12IntsAcc::r12t1);
	r12intsacc->release_pair_block(l,k,R12IntsAcc::r12t2);
      }
  }
  delete[] proc_with_ints;
  tim_exit("mp2-r12a intermeds");
  r12intsacc->deactivate();
  if (debug_)
    ExEnv::out0() << indent << "Computed intermediates V, X, and T" << endl;

  if (nproc > 1) {
    // Use MemoryGrp to send all contributions to intermediates V, X, and T to node 0
    msg->sum(Vaa_ijkl,naa*naa,0,0);
    msg->sum(Vab_ijkl,nab*nab,0,0);
    msg->sum(Xaa_ijkl,naa*naa,0,0);
    msg->sum(Xab_ijkl,nab*nab,0,0);
    msg->sum(Taa_ijkl,naa*naa,0,0);
    msg->sum(Tab_ijkl,nab*nab,0,0);
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
  msg->sum(aoint_computed);

  if (me == 0 && mole->if_to_checkpoint() && r12intsacc->can_restart()) {
    StateOutBin stateout(mole->checkpoint_file());
    SavableState::save_state(mole,stateout);
    ExEnv::out0() << indent << "Checkpointed the wave function" << endl;
  }

#if PRINT_BIGGEST_INTS
  biggest_ints_1.combine(msg);
  biggest_ints_2.combine(msg);
  biggest_ints_2s.combine(msg);
  biggest_ints_3a.combine(msg);
  biggest_ints_3.combine(msg);
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
		    << npass*nshell*(nshell+1)*nshell*nshell_aux/2
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
  
  for (int i=0; i<thr->nthread(); i++) {
    delete e123thread[i];
  }
  delete[] e123thread;
  
  
  delete[] tbints; tbints = 0;
  delete[] scf_vector;
  delete[] scf_vector_dat;
  delete[] evals;
  tim_exit("r12a-abs-mem");
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
R12IntEval_abs_A::compute_transform_batchsize_(size_t mem_alloc, size_t mem_static, int nocc_act, const int num_te_types)
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
R12IntEval_abs_A::compute_transform_dynamic_memory_(int ni, int nocc_act, const int num_te_types)
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
  int nbasis_aux = r12info()->basis_aux()->nbasis();
  int nfuncmax_aux = r12info()->basis_aux()->max_nfunction_in_shell();
  int nthread = r12info()->thr()->nthread();
  int nocc = r12info()->nocc();

  distsize_t memsize = sizeof(double)*(num_te_types*((distsize_t)nthread * ni * nbasis * nfuncmax * nfuncmax_aux // iqrs
						     + (distsize_t)ni * nocc * nfuncmax_aux * nfuncmax_aux  // ikrs
						     + (distsize_t)nij * nocc * nbasis_aux // ikjs - buffer of 3 q.t. and higher
						     // transformed integrals
						     )
				       + (distsize_t)nocc * nfuncmax_aux // ks
				       + (distsize_t)nocc * nbasis_aux // kx
				       );
  return memsize;
}



////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
