//
// csgrad.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Ida Nielsen <ida@kemi.aau.dk>
// Maintainer: LPS
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

#include <util/misc/formio.h>
#include <util/misc/regtime.h>
#include <util/group/memory.h>
#include <util/group/message.h>
#include <util/class/class.h>
#include <util/state/state.h>
#include <math/scmat/matrix.h>
#include <math/scmat/blocked.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/mbpt/bzerofast.h>
#include <chemistry/qc/mbpt/mbpt.h>
#include <chemistry/qc/mbpt/util.h>
#include <chemistry/qc/mbpt/csgrade12.h>
#include <chemistry/qc/mbpt/csgrad34qb.h>
#include <chemistry/qc/mbpt/csgrads2pdm.h>

using namespace std;
using namespace sc;

#define SINGLE_THREAD_E12   0
#define SINGLE_THREAD_QBT34 0
#define SINGLE_THREAD_S2PDM 0

#define PRINT2Q 0
#define PRINT3Q 0
#define PRINT4Q 0
#if PRINT_BIGGEST_INTS
BiggestContribs biggest_ints_1(4,40);
#endif

#define WRITE_DOUBLES 0

static void sum_gradients(const Ref<MessageGrp>& msg, double **f, int n1, int n2);
static void zero_gradients(double **f, int n1, int n2);
static void accum_gradients(double **g, double **f, int n1, int n2);

#define PRINT1Q 0

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

void
MBPT2::compute_cs_grad()
{

  // New version of MP2 gradient program which uses the full
  // permutational symmetry of the two-electron integral derivatives

  Ref<SCMatrixKit> kit = basis()->matrixkit();

  int do_d2_ = 1;  // if true, compute d2 diagnostic

  int nij;        // number of i,j pairs on a node (for e.g., mo_int)
  double *mo_int; // MO integrals of type (ov|ov)
                  // (and these integrals divided by
                  // orbital energy denominators)
  double *integral_iqjs; // half-transformed integrals

  int nocc_act, nvir_act;
  int i, j, k;
  int ii, bb;
  int x, y;
  int a, b, c;
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
  int aointder_computed = 0;
  int xyz;
  int natom = molecule()->natom();     // the number of atoms
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
  double ecorr_mp2 = 0.0;
  double escf;
  double emp2=0.0;
  int tol;                    // log2 of the erep tolerance
                              // (erep < 2^tol => discard)
  double *Wkj=0,*Wab=0,*Waj=0;// occ-occ, vir-vir and vir-occ parts of
                              // second order correction to MP2
                              // energy weighted density matrix
  double *Pkj=0,*Pab=0;       // occ-occ and vir-vir parts of second order
                              // correction to MP2 density matrix
  double *d2occ_mat, *d2vir_mat; // matrices for computation of D2 diagnostic
  double *Laj=0;              // MP2 Lagrangian
  double *Lpi;                // contrib to MP2 Lagrangian partially in AO basis
  double *pkj_ptr=0, *pab_ptr;
  double *d2occ_mat_ptr;
  double *d2vir_mat_ptr;
  double *wkj_ptr, *wjk_ptr, *wab_ptr, *wba_ptr, *waj_ptr=0;
  double *laj_ptr, *lpi_ptr, *lqi_ptr;
  double *gamma_iajs, *gamma_iajs_tmp;
                              // partially back-transformed non-sep 2PDM's
  double *gamma_iqjs_tmp;
  double *gamma_iajs_ptr;
  double *gamma_iqjs_ptr;
  double *gammabuf;           // buffer used for sending elements of gamma_iqjs
  double *mo_intbuf;          // buffer used for sending mo integrals
  double tmpval, tmpval1;
  double *P2AO, *W2AO;
  double *p2ao_ptr, *w2ao_ptr;
  double *PHF, *WHF;
  double *phf_ptr, *whf_ptr;
  double *PMP2, *WMP2;
  double *pmp2_ptr, *wmp2_ptr;

  double *ixjs_tmp;      // three-quarter transformed two-el integrals
  double *integral_ixjs;  // all three-quarter transformed two-el integrals
  double *integral_iajy; // mo integrals (y = any MO)
  double *integral_ikja; // mo integrals
  double *integral_iqjs_ptr;
  double *iajy_ptr;
  double *ixjs_ptr;
  double *ikja_ptr;
  double *iajs_ptr, *ikjs_ptr;

  double **gradient=0, *gradient_dat=0;  // The MP2 gradient
  double **hf_gradient=0, *hf_gradient_dat=0;  // The HF gradient
  double **ginter=0;    // Intermediates for the MP2 gradient
  double **hf_ginter=0;    // Intermediates for the HF gradient
  double d2o, d2v, d2_diag;

  BiggestContribs biggest_coefs(5,10);
  CharacterTable ct = molecule()->point_group()->char_table();

#if PRINT_BIGGEST_INTS
  BiggestContribs biggest_ints_2(4,40);
  BiggestContribs biggest_ints_2s(4,40);
  BiggestContribs biggest_ints_3a(4,40);
  BiggestContribs biggest_ints_3(4,40);
#endif

  int dograd = gradient_needed();

  Timer tim("mp2-mem");

  nfuncmax = basis()->max_nfunction_in_shell();

  nshell = basis()->nshell();

  me = msg_->me();

  if (me == 0) {
    ExEnv::out0() << endl << indent
         << "Entered memgrp based MP2 routine" << endl;
    }

  nproc = msg_->n();
  if (me == 0)
    ExEnv::out0() << indent << scprintf("nproc = %i", nproc) << endl;

  tol = (int) (-10.0/log10(2.0));  // discard ereps smaller than 10^-10

  nocc = 0;
  for (i=0; i<oso_dimension()->n(); i++) {
    if (reference_->occupation(i) == 2.0) nocc++;
    }

  nocc_act = nocc - nfzc;
  nvir  = noso - nocc;
  nvir_act = nvir - nfzv;

  // Do a few preliminary tests to make sure the desired calculation
  // can be done (and appears to be meaningful!)

  if (nocc_act <= 0) {
    if (me == 0) {
      ExEnv::err0() << "There are no active occupied orbitals; program exiting" << endl;
      }
    abort();
    }

  if (nvir_act <= 0) {
    if (me == 0) {
      ExEnv::err0() << "There are no active virtual orbitals; program exiting" << endl;
      }
    abort();
    }

  if (restart_orbital_memgrp_) {
    if (!dograd && !do_d1_ && !do_d2_) {
      ExEnv::out0() << indent
           << scprintf("Restarting at orbital %d with partial energy %18.14f",
                       restart_orbital_memgrp_, restart_ecorr_)
           << endl;
      ecorr_mp2 = restart_ecorr_;
      }
    else {
      ExEnv::out0() << indent
           << "Restart requested but not possible with gradients, D1, or D2"
           << endl;
      restart_ecorr_ = 0.0;
      restart_orbital_memgrp_ = 0;
      }
    }
  else {
      restart_ecorr_ = 0.0;
    }

  ////////////////////////////////////////////////////////
  // Compute batch size ni for mp2 loops;
  //
  // The following arrays are kept throughout (all of type double):
  //   scf_vector, gradient, ginter, Pkj, Pab, Wkj, Wab, Waj, Laj
  // and memory allocated for these arrays  and integral evaluators
  // is called mem_static
  //
  ////////////////////////////////////////////////////////
  if (me == 0) {
    mem_static = nbasis*noso; // scf vector
    mem_static += 2*nbasis*nfuncmax; // iqjs & iqjr
    if (dograd) {
      mem_static += 9*natom; // gradient & ginter & hf_ginter
      mem_static += (nocc*(nocc+1))/2; // Pkj
      mem_static += (nvir*(nvir+1))/2; // Pab
      mem_static += nocc*nocc; // Wkj
      mem_static += nvir*nvir; // Wab
      mem_static += 2*nocc*nvir; // Waj & Laj
      if (do_d2_) {
        mem_static += (nocc_act*(nocc_act+1))/2; // d2occ_mat
        mem_static += (nvir_act*(nvir_act+1))/2; // d2vir_mat
        }
      }
    else if (do_d1_) {
      mem_static += nocc*nvir; // partial Laj
      }
    mem_static *= sizeof(double);
    int nthreads = thr_->nthread();
    mem_static += nthreads * integral()->storage_required_eri(basis()); // integral evaluators
    ni = compute_cs_batchsize(mem_static, nocc_act-restart_orbital_memgrp_);
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
  double dmem_static = mem_static;
  msg_->bcast(dmem_static);
  mem_static = size_t(dmem_static);

  // Compute the storage to be used by the integral routines (required plus optional)
  size_t dyn_mem = distsize_to_size(compute_cs_dynamic_memory(ni,nocc_act));
  int mem_remaining;
  if (mem_alloc <= (dyn_mem + mem_static)) mem_remaining = 0;
  else mem_remaining = mem_alloc - dyn_mem - mem_static;
  mem_remaining += thr_->nthread() * integral()->storage_required_eri(basis());

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
       << compute_cs_dynamic_memory(nocc_act,nocc_act)+mem_static
       << " Bytes"
       << endl;
  ExEnv::out0() << indent
       << "Minimum memory required:        "
       << compute_cs_dynamic_memory(1,nocc_act)+mem_static
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

  if (ni == nocc_act-restart_orbital_memgrp_) {
    npass = 1;
    rest = 0;
    }
  else {
    rest = (nocc_act-restart_orbital_memgrp_)%ni;
    npass = (nocc_act-restart_orbital_memgrp_ - rest)/ni + 1;
    if (rest == 0) npass--;
    }

  if (me == 0) {
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
    }

  int nijmax = 0;
  index = 0;
  for (i=0; i<ni; i++) {
      for (j=0; j<nocc; j++) {
          if (index++ % nproc == me) nijmax++;
        }
    }

  ////////////////////////////////////////////////
  // The scf vector is distributed between nodes;
  // put a copy of the scf vector on each node;
  ////////////////////////////////////////////////

  escf = reference_->energy();
  hf_energy_ = escf;

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

  //////////////////////////////////////////////////////////////
  // Allocate storage for various arrays needed for gradients
  // (Pkj and Pab are symmetric, so store only lower triangle)
  //////////////////////////////////////////////////////////////

  if (dograd) {
    Pkj            = (double*) malloc((nocc*(nocc+1)/2)*sizeof(double));
    Pab            = (double*) malloc((nvir*(nvir+1)/2)*sizeof(double));
    Wkj            = (double*) malloc(nocc*nocc*sizeof(double));
    Wab            = (double*) malloc(nvir*nvir*sizeof(double));
    Waj            = (double*) malloc(nvir*nocc*sizeof(double));
    if (do_d2_) {
      d2occ_mat =  new double[nocc_act*(nocc_act+1)/2];
      d2vir_mat =  new double[nvir_act*(nvir_act+1)/2];
      }

    gradient_dat = new double[natom*3];
    gradient = new double*[natom];
    for (i=0; i<natom; i++) {
      gradient[i] = &gradient_dat[i*3];
      }

    hf_gradient_dat = new double[natom*3];
    hf_gradient = new double*[natom];
    for (i=0; i<natom; i++) {
      hf_gradient[i] = &hf_gradient_dat[i*3];
      }

    ginter = new double*[natom];
    for (i=0; i<natom; i++) {
      ginter[i] = new double[3];
      for (xyz=0; xyz<3; xyz++) ginter[i][xyz] = 0;
      }

    hf_ginter = new double*[natom];
    for (i=0; i<natom; i++) {
      hf_ginter[i] = new double[3];
      for (xyz=0; xyz<3; xyz++) hf_ginter[i][xyz] = 0;
      }

    //////////////////////////////
    // Initialize various arrays
    //////////////////////////////

    bzerofast(Pkj,nocc*(nocc+1)/2);
    bzerofast(Wkj,nocc*nocc);
    bzerofast(Pab,nvir*(nvir+1)/2);
    bzerofast(Wab,nvir*nvir);
    bzerofast(Waj,nvir*nocc);
    if (do_d2_) {
      bzerofast(d2occ_mat,nocc_act*(nocc_act+1)/2);
      bzerofast(d2vir_mat,nvir_act*(nvir_act+1)/2);
      }

    if (me == 0) zero_gradients(gradient, natom, 3);
    if (me == 0) zero_gradients(hf_gradient, natom, 3);
    }

  if (dograd || do_d1_) {
    Laj = (double*) malloc(nvir*nocc*sizeof(double));
    bzerofast(Laj,nvir*nocc);
    }

  if (debug_ > 2 && me == 0) {
    for (j=0; j<noso; j++) {
      ExEnv::out0() << indent
           << scprintf("eigenvalue[%3d] = %15.10lf", j, evals[j]);
      if (j < nfzc) ExEnv::out0() << " (frozen docc)";
      else if (j < nocc_act + nfzc) ExEnv::out0() << " (active docc)";
      else if (j < nvir_act + nocc_act + nfzc) ExEnv::out0() << " (active uocc)";
      else ExEnv::out0() << " (frozen uocc)";
      ExEnv::out0() << endl;
      }
    }

  /////////////////////////////////////
  //  Begin MP2 loops
  /////////////////////////////////////

  // debug print
  if (debug_ && me == 0) {
    ExEnv::out0() << indent
         << scprintf("node %i, begin loop over i-batches",me) << endl;
    }
  // end of debug print

  // Initialize the integrals
  integral()->set_storage(mem_remaining);
  tbints_ = new Ref<TwoBodyInt>[thr_->nthread()];
  for (i=0; i<thr_->nthread(); i++) {
      tbints_[i] = integral()->electron_repulsion();
    }
  if (dograd || do_d1_) {
    tbintder_ = new Ref<TwoBodyDerivInt>[thr_->nthread()];
    for (i=0; i<thr_->nthread(); i++) {
      tbintder_[i] = integral()->electron_repulsion_deriv();
      }
    }

  int mem_integral_intermediates = integral()->storage_used();
  int mem_integral_storage = (mem_remaining - mem_integral_intermediates) / thr_->nthread();
  if (mem_integral_storage<0) mem_integral_storage = 0;
  for (i=0; i<thr_->nthread(); i++) {
      tbints_[i]->set_integral_storage(mem_integral_storage);
    }

  ExEnv::out0() << endl << indent
       << scprintf("Memory used for integral intermediates: %i Bytes",
                   mem_integral_intermediates)
       << endl;
  ExEnv::out0() << indent
       << scprintf("Memory used for integral storage:       %i Bytes",
                   mem_integral_storage)
       << endl;

  if (mem == 0) {
      ExEnv::errn() << "MBPT2: memory group not initialized" << endl;
      abort();
    }

  mem->set_localsize(size_t(nijmax)*nbasis*nbasis*sizeof(double));
  ExEnv::out0() << indent
       << "Size of global distributed array:       "
       << mem->totalsize()
       << " Bytes"
       << endl;

  MemoryGrpBuf<double> membuf_remote(mem);

  int usep4 = !dograd;

  Ref<ThreadLock> lock = thr_->new_lock();
  CSGradErep12Qtr** e12thread = new CSGradErep12Qtr*[thr_->nthread()];
  DistShellPair::SharedData sp_e_data, sp_g_data;
  for (i=0; i<thr_->nthread(); i++) {
    e12thread[i] = new CSGradErep12Qtr(i, thr_->nthread(), me, nproc,
                                       mem, msg_, lock, basis(), tbints_[i],
                                       nocc, scf_vector, tol, debug_,
                                       dynamic_, print_percent_,
                                       &sp_e_data, usep4);
    }

    CSGrad34Qbtr** qbt34thread;
    if (dograd || do_d1_) {
      qbt34thread = new CSGrad34Qbtr*[thr_->nthread()];
      for (i=0; i<thr_->nthread(); i++) {
        qbt34thread[i] = new CSGrad34Qbtr(i, thr_->nthread(), me, nproc,
                                          mem, msg_, lock, basis(), tbints_[i],
                                          tbintder_[i], nocc, nfzc, scf_vector,
                                          tol, debug_, dynamic_, print_percent_,
                                          &sp_g_data, dograd, natom);
        }
      }

  tim.enter("mp2 passes");
  for (pass=0; pass<npass; pass++) {

    if (me == 0) {
      ExEnv::out0() << indent << "Beginning pass " << pass+1 << endl;
      }

    i_offset = restart_orbital_memgrp_ + pass*ni + nfzc;
    if ((pass == npass - 1) && (rest != 0)) ni = rest;

    // Compute number of of i,j pairs on each node for
    // two-el integrals and non-sep 2PDM elements
    index = 0;
    nij = 0;
    for (i=0; i<ni; i++) {
      for (j=0; j<nocc; j++) {
        if (index++ % nproc == me) nij++;
        }
      }

    // debug print
    if (debug_)
      ExEnv::outn() << indent << "node " << me << ", nij = " << nij << endl;
    // end of debug print

    mem->sync(); // This must be here or gamma non-sep will be wrong when running
                 // on multiple processors with more than one pass

    r_offset = 0;

    // Allocate and initialize some arrays
    // (done here to avoid having these arrays
    // overlap with arrays allocated later)

    // Allocate (and initialize) some arrays

    integral_iqjs = (double*) mem->localdata();

    bzerofast(integral_iqjs, nij*nbasis*nbasis);

    integral_iqjs = 0;

    mem->sync();

    index = 0;

    if (me == 0) {
      ExEnv::out0() << indent
           << scprintf("Begin loop over shells (erep, 1.+2. q.t.)") << endl;
      }

    // Do the two eletron integrals and the first two quarter transformations
    tim.enter("erep+1.qt+2.qt");
    sp_e_data.init();
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
    tim.exit("erep+1.qt+2.qt");

    if (me == 0) {
      ExEnv::out0() << indent << "End of loop over shells" << endl;
      }

    mem->sync();  // Make sure iqjs is complete on each node before continuing

    integral_iqjs = (double*) mem->localdata();

#if PRINT2Q
    if (me == 0) {
      int index = 0;
      int ij_index = 0;
      for (int i = 0; i<ni; i++) {
	for (int j = 0; j<nocc; j++) {
	  if (index++ % nproc == me) {
	    if (j >= nfzc) {
	      double *integral_ij_offset = integral_iqjs + nbasis*nbasis*ij_index;
	      for (int s = 0; s<nbasis; s++) {
		double *integral_ijsq_ptr = integral_ij_offset + s*nbasis;
		for (int q = 0; q<nbasis; q++) {
		  printf("2Q: (%d %d|%d %d) = %12.8f\n",
			 i,q,j,s,*integral_ijsq_ptr);
		  *integral_ijsq_ptr++;
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
    ixjs_tmp = new double[nbasis];

    if (me == 0) {
      ExEnv::out0() << indent << "Begin third q.t." << endl;
      }

    tim.enter("3. q.t.");
    // Begin third quarter transformation;
    // generate (ix|js) for i act, j act or frz, and x any MO (act or frz)
    index = 0;
    ij_index = 0;
    for (i=0; i<ni; i++) {
      for (j=0; j<nocc; j++) {
        if (index++ % nproc == me) {

          for (s=0; s<nbasis; s++) {

            bzerofast(ixjs_tmp, nbasis);
            for (q=0; q<nbasis; q++) {
              integral_iqjs_ptr = &integral_iqjs[q + nbasis*(s + nbasis*ij_index)];
              ixjs_ptr = ixjs_tmp;
              c_qx = scf_vector[q];
              tmpval = *integral_iqjs_ptr;
#if PRINT_BIGGEST_INTS
              biggest_ints_2.insert(tmpval,i+i_offset,j,s,q);
              if ((i+i_offset==104 && j == 1)
                  ||(i+i_offset==104 && j == 2)) {
                biggest_ints_2s.insert(tmpval,i+i_offset,j,s,q);
                }
#endif
              for (x=0; x<noso; x++) {
                *ixjs_ptr++ += *c_qx++ * tmpval;
                }
              }   // exit q loop

            // Put ixjs into integral_iqjs, while overwriting what was there;
            // i.e., integral_iqjs will now contain three-quarter transformed
            // integrals ixjs
            integral_iqjs_ptr = &integral_iqjs[nbasis*(s + nbasis*ij_index)];
            ixjs_ptr = ixjs_tmp;
            for (x=0; x<noso; x++) {
#if PRINT_BIGGEST_INTS
              if (x>=nocc) {
                biggest_ints_3a.insert(*ixjs_ptr,i+i_offset,j,s,x-nocc);
                }
#endif
              *integral_iqjs_ptr++ = *ixjs_ptr++;
              }
            }   // exit s loop
          ij_index++;
          }     // endif
        }       // exit j loop
      }         // exit i loop
    // end of third quarter transformation
    tim.exit("3. q.t.");

    if (me == 0) {
      ExEnv::out0() << indent << "End of third q.t." << endl;
      }

    delete[] ixjs_tmp;

    // The array of half-transformed integrals integral_iqjs has now
    // been overwritten by three-quarter transformed integrals ixjs;
    // rename the array integral_ixjs, where x = any MO
    integral_ixjs = integral_iqjs;

#if PRINT3Q
    if (me == 0) {
      int index = 0;
      int ij_index = 0;
      for (int i = 0; i<ni; i++) {
	for (int j = 0; j<nocc; j++) {
	  if (index++ % nproc == me) {
	    if (j >= nfzc) {
	      double *integral_ij_offset = integral_ixjs + nbasis*nbasis*ij_index;
	      for (int s = 0; s<nbasis; s++) {
		double *integral_ijsx_ptr = integral_ij_offset + s*nbasis;
		for (int x = 0; x<noso; x++) {
		  printf("3Q: (%d %d|%d %d) = %12.8f\n",
			 i,x,j,s,*integral_ijsx_ptr);
		  *integral_ijsx_ptr++;
		}
	      }
	    }
	    ij_index++;
	  }
	}
      }
    }
#endif

    integral_iajy = new double[noso];
    // in iajy: i act; a,j act or frz; y act or frz occ or virt.
    integral_ikja = new double[nvir_act];
    // in ikja: i,j act; k act or frz; a act.

    if (me == 0) {
      ExEnv::out0() << indent << "Begin fourth q.t." << endl;
      }

    // Begin fourth quarter transformation
    // generating MO integrals (ov|ov), (ov|oo) and (oo|ov)
    tim.enter("4. q.t.");
    index = 0;
    ij_index = 0;
    for (i=0; i<ni; i++) {
      for (j=0; j<nocc; j++) {
        if (index++ % nproc == me) {

          for (a=0; a<nvir; a++) {
            bzerofast(integral_iajy, noso);
            iajs_ptr = &integral_ixjs[a+nocc + nbasis*nbasis*ij_index];
            for (s=0; s<nbasis; s++) {
              c_sy = scf_vector[s];
              iajy_ptr = integral_iajy;
              tmpval = *iajs_ptr;
#if PRINT_BIGGEST_INTS
              biggest_ints_3.insert(tmpval,i+i_offset,j,s,a);
              if ((i+i_offset==105 && j == 2 && s == 170 && a == 3)
                  ||(i+i_offset==102 && j == 2 && s == 170 && a == 2)) {
                ExEnv::outn() << scprintf("3/4: %3d %3d %3d %3d: %16.10f",
                                 i+i_offset, j, s, x-nocc)
                     << endl;
                }
#endif
              for (y=0; y<noso; y++) {
                *iajy_ptr++ += *c_sy++ * tmpval;
                } // exit y loop
              iajs_ptr += nbasis;
              }   // exit s loop
            // Put integral_iajy into ixjs for one i,a,j while
            // overwriting elements of ixjs
            iajs_ptr = &integral_ixjs[a+nocc + nbasis*nbasis*ij_index];
            iajy_ptr = integral_iajy;
            for (y=0; y<noso; y++) {
              *iajs_ptr = *iajy_ptr++;
              iajs_ptr += nbasis;
              } // exit y loop
            }   // exit a loop

          // this is only needed for gradients
          if ((dograd || do_d1_) && j >= nfzc) {
            for (k=0; k<nocc; k++) {
              bzerofast(integral_ikja, nvir_act);
              ikjs_ptr = &integral_ixjs[k + nbasis*nbasis*ij_index];
              for (s=0; s<nbasis; s++) {
                c_sa = &scf_vector[s][nocc];
                ikja_ptr = integral_ikja;
                tmpval = *ikjs_ptr;
                for (a=0; a<nvir_act; a++) {
                  *ikja_ptr++ += *c_sa++ * tmpval;
                  } // exit a loop
                ikjs_ptr += nbasis;
                }   // exit s loop
              // Put integral_ikja into ixjs for one i,k,j while
              // overwriting elements of ixjs
              ikjs_ptr = &integral_ixjs[k + nbasis*(nocc + nbasis*ij_index)];
              ikja_ptr = integral_ikja;
              for (a=0; a<nvir_act; a++) {
                *ikjs_ptr = *ikja_ptr++;
                ikjs_ptr += nbasis;
                } // exit a loop
              }   // exit k loop
            }     //endif

          ij_index++;
          }   // endif
        }     // exit j loop
      }       // exit i loop
    // end of fourth quarter transformation
    tim.exit("4. q.t.");

    if (me == 0) {
      ExEnv::out0() << indent << "End of fourth q.t." << endl;
      }

    // The array integral_ixjs has now been overwritten by MO integrals
    // iajy and ikja, so rename the array mo_int
    mo_int = integral_ixjs;

    delete[] integral_iajy;
    delete[] integral_ikja;

    // Divide the (ia|jb) MO integrals by the term
    // evals[i]+evals[j]-evals[a]-evals[b]
    // and keep these integrals in mo_int;
    // do this only for active i, j, a, and b
    tim.enter("divide (ia|jb)'s");

    index = 0;
    ij_index = 0;
    for (i=0; i<ni; i++) {
      ii = i+i_offset;
      for (j=0; j<nocc; j++) {
        if (index++ % nproc == me) {
          if (j>=nfzc) {
            for (b=0; b<nvir_act; b++) {
              iajb_ptr = &mo_int[nocc + nbasis*(b+nocc + nbasis*ij_index)];
              bb = b+nocc;
              for (a=0; a<nvir_act; a++) {
#if PRINT4Q
	      printf("4Q: (%d %d|%d %d) = %12.8f\n",
		     i,a+nocc,j,b+nocc,*iajb_ptr);
#endif
		// Zero out nonsymmetric integral, else divide by denominators
	        if (usep4 && ( symorb_irrep_[ii] ^
		      symorb_irrep_[j] ^
		      symorb_irrep_[a+nocc] ^
		      symorb_irrep_[bb]) ) {
		  *iajb_ptr++ = 0.0;
		}
		else
		  *iajb_ptr++ /= evals[ii]+evals[j]-evals[a+nocc]-evals[bb];
                } // exit a loop
              }   // exit b loop
            }     // endif
          ij_index++;
          }       // endif
        }         // exit j loop
      }           // exit i loop
    tim.exit("divide (ia|jb)'s");

    // We now have the fully transformed integrals (ia|jb)
    // divided by the proper orbital energy denominators
    // for one batch of i, all j<nocc, and all a<nvir and b<nvir,
    // where i, j, a, and b are all active;
    // compute contribution to the MP2 correlation energy
    // from these integrals

#if WRITE_DOUBLES
    if (nproc > 1 || npass > 1) {
      ExEnv::outn() << "csgrad.cc: WRITE_DOUBLES set but case not allowed" << endl;
      abort();
      }
    ExEnv::outn() << "csgrad.cc: WRITING DOUBLES: CHECK ORDER" << endl;
    char *doutname = SCFormIO::fileext_to_filename(".mp2");
    FILE *dout = fopen(doutname,"w");
    delete[] doutname;
    fwrite(&nocc_act, sizeof(int), 1, dout);
    fwrite(&nvir_act, sizeof(int), 1, dout);
    for (j=nfzc; j<nocc; j++) {
      for (b=0; b<nvir_act; b++) {
        for (i=0; i<ni; i++) {
          ij_index = nocc*i + j;
          iajb_ptr = &mo_int[nocc + nbasis*(b+nocc + nbasis*ij_index)];
          fwrite(iajb_ptr, sizeof(double), nvir_act, dout);
          for (a=0; a<nvir_act; a++) {
            if (fabs(iajb_ptr[a])>1.0e-8) {
              ExEnv::outn() << scprintf(" Djbia(%2d %2d %2d %2d) = %12.8f",
                               j+1-nfzc,b+1,i+1,a+1,iajb_ptr[a])
                   << endl;
              }
            }
          }
        }
      }
    fclose(dout);
#endif

    tim.enter("compute ecorr");

    index = 0;
    ij_index = 0;
    for (i=0; i<ni; i++) {
      double ecorr_i = 0.0;
      for (j=0; j<nocc; j++) {
        double ecorr_ij = 0.0;
        if (index++ % nproc == me) {

          if (j>=nfzc) {
            for (b=0; b<nvir_act; b++) {
              iajb_ptr = &mo_int[nocc + nbasis*(b+nocc + nbasis*ij_index)];
              ibja_ptr = &mo_int[b+nocc + nbasis*(nocc + nbasis*ij_index)];
              for (a=0; a<nvir_act; a++) {
                delta_ijab = evals[i_offset+i]+evals[j]-evals[nocc+a]-evals[nocc+b];
                // only include determinants with unique coefficients
                if (a>=b && i_offset+i>=j) {
                  if (a>b && i_offset+i>j) {
                    // aaaa or bbbb
                    biggest_coefs.insert(*iajb_ptr - *ibja_ptr,
                                         i_offset+i,j,a,b,1111);
                    // aabb or bbaa or abba or baab
                    biggest_coefs.insert(*ibja_ptr,i_offset+i,j,b,a,1212);
                    } // endif
                  // aabb or bbaa or abba or baab
                  biggest_coefs.insert(*iajb_ptr,i_offset+i,j,a,b,1212);
                  } // endif

                tmpval = *iajb_ptr*(2**iajb_ptr - *ibja_ptr)*delta_ijab;
                ecorr_mp2 += tmpval;
                if (debug_) ecorr_ij += tmpval;
                iajb_ptr++;
                ibja_ptr += nbasis;;
                } // exit a loop
              }   // exit b loop
            }     // endif
          ij_index++;
          }       // endif
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
    tim.exit("compute ecorr");

    // debug print
    if (debug_ && me == 0) {
      ExEnv::out0() << indent << "End of ecorr" << endl;
      }
    // end of debug print

    if (npass > 1 && pass < npass - 1) {
      double passe = ecorr_mp2;
      msg_->sum(passe);
      ExEnv::out0() << indent
           << "Partial correlation energy for pass " << pass << ":" << endl;
      ExEnv::out0() << indent
           << scprintf("  restart_ecorr          = %18.14f", passe)
           << endl;
      ExEnv::out0() << indent
           << scprintf("  restart_orbital_memgrp = %d", ((pass+1) * ni))
           << endl;
      }

    integral_iqjs = 0;
    mem->sync(); // Make sure MO integrals are complete on all nodes before continuing

    // don't go beyond this point if only the energy is needed
    if (!dograd && !do_d1_) continue;

    mo_int = (double*) mem->localdata();

    if (!dograd) goto compute_L;

    // Update the matrices Pkj and Wkj with
    // contributions from (occ vir|occ vir) integrals
    index = 0;
    ij_index = 0;
    tim.enter("Pkj and Wkj");
    for (i=0; i<ni; i++) {
      for (j=0; j<nocc; j++) {
        if (index++ % nproc == me) {

          if (j>=nfzc) {
            for (kloop=me; kloop<me+nocc; kloop++) {
              // stagger k's to minimize contention
              k = kloop%nocc;
              if (k<=j) pkj_ptr = &Pkj[j*(j+1)/2 + k];
              if (do_d2_) {
                if (k<=j && k>=nfzc) {
                  d2occ_mat_ptr = &d2occ_mat[(j-nfzc)*(j-nfzc+1)/2 + k-nfzc];
                  }
                }
              wjk_ptr = &Wkj[j*nocc + k];
              // Send for iakb, if necessary
              ik_index = (i*nocc + k)/nproc;
              ik_proc = (i*nocc + k)%nproc;
              ik_offset = nocc + nocc*nbasis + nbasis*nbasis*ik_index;
              mo_intbuf = (double*) membuf_remote.readonly_on_node(ik_offset,
                                                                   nbasis*nvir-nocc,
                                                                   ik_proc);
              for (a=0; a<nvir_act; a++) {
                ibja_ptr = &mo_int[nocc + nbasis*(a+nocc + nbasis*ij_index)];
                iajb_ptr = &mo_int[a+nocc + nbasis*(nocc + nbasis*ij_index)];
                iakb_ptr = &mo_intbuf[a];
                for (b=0; b<nvir_act; b++) {
                  tmpval1 = *iakb_ptr * *iajb_ptr;
                  tmpval = 2**iakb_ptr * *ibja_ptr++ - 4*tmpval1;
                  // tmpval = 2**iakb_ptr * (*ibja_ptr++ - 2 * *iajb_ptr);
                  iakb_ptr += nbasis;
                  iajb_ptr += nbasis;
                  if (k<nfzc && k<=j) *pkj_ptr += tmpval/(evals[k]-evals[j]);
                  else if (k<=j) {
                    *pkj_ptr += tmpval;
                    if (do_d2_) *d2occ_mat_ptr += tmpval1;
                    }
                  if (k>=nfzc) {
                    delta_ijab = evals[i_offset+i]+evals[j]-evals[nocc+a]-evals[nocc+b];
                    *wjk_ptr += tmpval*delta_ijab;
                    }
                  } // exit b loop
                }   // exit a loop
              mo_intbuf = 0;
              membuf_remote.release();
              }     // end kloop loop
            }       // endif

          ij_index++;
          }         // endif
        }           // exit j loop
      }             // exit i loop
    tim.exit("Pkj and Wkj");

    // debug print
    if (debug_ && me == 0) {
      ExEnv::out0() << indent << "End of Pkj and Wkj" << endl;
      }
    // end of debug print

    // Update the matrices Pab and Wab with
    // contributions from (occ vir|occ vir) integrals
    tim.enter("Pab and Wab");
    index = 0;
    ij_index = 0;
    for (i=0; i<ni; i++) {
      for (j=0; j<nocc; j++) {
        if (index++ % nproc == me) {
          if (j>=nfzc) {

            offset = nocc + nocc*nbasis + nbasis*nbasis*ij_index;
            for (a=0; a<nvir_act; a++) {
              pab_ptr = &Pab[a*(a+1)/2];
              if (do_d2_) d2vir_mat_ptr = &d2vir_mat[a*(a+1)/2];
              for (b=0; b<=a; b++) {  // active-active part of Pab and Wab
                wab_ptr = &Wab[a*nvir + b];
                wba_ptr = &Wab[b*nvir + a];
                ibjc_ptr = &mo_int[offset + b];
                icjb_ptr = &mo_int[offset + b*nbasis];
                iajc_ptr = &mo_int[offset + a];
                icja_ptr = &mo_int[offset + a*nbasis];
                for (c=0; c<nvir_act; c++) {
                  tmpval = *iajc_ptr**ibjc_ptr;
                  if (do_d2_) *d2vir_mat_ptr += tmpval;
                  *pab_ptr += 4*tmpval - 2**iajc_ptr * *icjb_ptr;
                  // *pab_ptr += 2**iajc_ptr * (2 * *ibjc_ptr - *icjb_ptr);

                  if (a == b) {
                    delta_ijac = evals[i_offset+i]+evals[j]-evals[nocc+a]-evals[nocc+c];
                    *wab_ptr += 2**ibjc_ptr*(*icja_ptr - 2 * *iajc_ptr)*delta_ijac;
                    }
                  else {
                    delta_ijbc = evals[i_offset+i]+evals[j]-evals[nocc+b]-evals[nocc+c];
                    delta_ijac = evals[i_offset+i]+evals[j]-evals[nocc+a]-evals[nocc+c];
                    *wab_ptr += 2**ibjc_ptr * (*icja_ptr - 2 * *iajc_ptr)*delta_ijac;
                    *wba_ptr += 2**iajc_ptr * (*icjb_ptr - 2 * *ibjc_ptr)*delta_ijbc;
                    }
                  iajc_ptr += nbasis;
                  ibjc_ptr += nbasis;
                  icja_ptr++;
                  icjb_ptr++;
                  } // exit c loop
                pab_ptr++;
                d2vir_mat_ptr++;
                }   // exit b loop
              }     // exit a loop

            for (a=0; a<nfzv; a++) {        // active-frozen part of Pab
              pab_ptr = &Pab[(a+nvir_act)*(a+nvir_act+1)/2];
              for (b=0; b<nvir_act; b++) {
                tmpval = evals[nocc+b] - evals[nocc+nvir_act+a];
                ibjc_ptr = &mo_int[offset+b];
                iajc_ptr = &mo_int[offset+a+nvir_act];
                icja_ptr = &mo_int[offset+(a+nvir_act)*nbasis];
                for (c=0; c<nvir_act; c++) {
                  *pab_ptr += 2**ibjc_ptr*(2**iajc_ptr - *icja_ptr++)/tmpval;
                  ibjc_ptr += nbasis;
                  iajc_ptr += nbasis;
                  }  // exit c loop
                pab_ptr++;
                }    // exit b loop
              }      // exit a loop

            }        // endif
          ij_index++;
          }     // endif
        }       // exit j loop
      }         // exit i loop
    tim.exit("Pab and Wab");

    // debug print
    if (debug_ && me == 0) {
      ExEnv::out0() << indent << "End of Pab and Wab" << endl;
      }
    // end of debug print

    compute_L:

    ///////////////////////////////////////
    // Update Waj and Laj with contrib.
    // from (oo|ov) and (ov|oo) integrals;
    // here a is active and j is active or
    // frozen
    ///////////////////////////////////////
    tim.enter("Waj and Laj");

    // (oo|ov) contribution
    index = 0;
    ik_index = 0;
    for (i=0; i<ni; i++) {
      for (k=0; k<nocc; k++) {
        if (index++ % nproc == me) {
          if (k>=nfzc) {
            offset = nbasis*nocc + nbasis*nbasis*ik_index;
            for (j=0; j<nocc; j++) {
              for (b=0; b<nvir_act; b++) {
                ibka_ptr = &mo_int[b+nocc + offset];
                ijkb_ptr = &mo_int[j + nbasis*b + offset];
                if (dograd) waj_ptr = &Waj[j*nvir]; // order as j*nvir+a to make loops more efficient
                laj_ptr = &Laj[j*nvir];
                for (a=0; a<nvir_act; a++) {
                  tmpval = 2**ibka_ptr * *ijkb_ptr;
                  ibka_ptr += nbasis;
                  if (dograd) *waj_ptr++ += tmpval;
                  *laj_ptr++ -= tmpval; // This term had the wrong sign in Frisch's paper
                  } // exit a loop
                }   // exit b loop
              }     // exit j loop
            }       // endif
          ik_index++;
          }         // endif
        }           // exit k loop
      }             // exit i loop

    // (ov|oo) contribution
    index = 0;
    ik_index = 0;
    for (i=0; i<ni; i++) {
      for (k=0; k<nocc; k++) {
        if (index++ % nproc == me) {
          if (k>=nfzc) {
            offset = nocc + nbasis*nbasis*ik_index;
            for (b=0; b<nvir_act; b++) {
              for (j=0; j<nocc; j++) {
                ibkj_ptr = &mo_int[offset + b + j*nbasis];
                ibka_ptr = &mo_int[offset + b + nocc*nbasis];
                if (dograd) waj_ptr = &Waj[j*nvir];
                laj_ptr = &Laj[j*nvir];
                for (a=0; a<nvir_act; a++) {
                  tmpval = 4 * *ibka_ptr * *ibkj_ptr;
                  ibka_ptr += nbasis;
                  if (dograd) *waj_ptr++ -= tmpval;
                  *laj_ptr++ += tmpval; // This term had the wrong sign in Frisch's paper
                  } // exit a loop
                }   // exit j loop
              }     // exit b loop
            }       // endif
          ik_index++;
          }       // endif
        }         // exit k loop
      }           // exit i loop

    tim.exit("Waj and Laj");
    /////////////////////////////
    // End of Waj and Laj update
    /////////////////////////////

    // debug print
    if (debug_ && me == 0) {
      ExEnv::out0() << indent << "End of Paj and Waj" << endl;
      }
    // end of debug print

    mo_int = 0;

    mem->sync(); // Need to synchronize before deleting mo_intbuf

    mo_int = (double*) mem->localdata();

    gamma_iajs_tmp = new double[nbasis*nvir_act];
    if (!gamma_iajs_tmp) {
      ExEnv::outn() << indent << "Could not allocate gamma_iajs_tmp" << endl;
      abort();
      }

    // debug print
    if (debug_ && me == 0) {
      ExEnv::out0() << indent << "Begin first and second q.b.t." << endl;
      }
    // end of debug print

    ///////////////////////////////////////////////////////////
    // Perform first and second quarter back-transformation.
    // Each node produces gamma_iajs, and gamma_iqjs
    // for a subset of i and j, all a and all s;
    // the back-transf. is done only for active i, j, a, and b
    ///////////////////////////////////////////////////////////

    // Begin first quarter back-transformation
    tim.enter("1. q.b.t.");
    index = 0;
    ij_index = 0;
    for (i=0; i<ni; i++) {
      for (j=0; j<nocc; j++) {
        if (index++ % nproc == me) {
          if (j>=nfzc) {
            bzerofast(gamma_iajs_tmp,nbasis*nvir_act);
            offset = nocc + nocc*nbasis + nbasis*nbasis*ij_index;

            for (a=0; a<nvir_act; a++) {
              for (s=0; s<nbasis; s++) {
                c_sb = &scf_vector[s][nocc];
                gamma_iajs_ptr = &gamma_iajs_tmp[s*nvir_act + a];
                ibja_ptr = &mo_int[a*nbasis + offset];
                iajb_ptr = &mo_int[a + offset];

                tmpval = 0.0;
                for (b=0; b<nvir_act; b++) {
                  tmpval += 2**c_sb++ * (2**iajb_ptr - *ibja_ptr++);
                  iajb_ptr += nbasis;
                  } // exit b loop
                *gamma_iajs_ptr += tmpval;
                }   // exit s loop
              }     // exit a loop
            // Put gamma_iajs_tmp into mo_int for one i,j
            // while overwriting mo_int
            gamma_iajs_ptr = gamma_iajs_tmp;
            for (y=0; y<nbasis; y++) {
              iajy_ptr = &mo_int[nocc + nbasis*(y + nbasis*ij_index)];
              for (a=0; a<nvir_act; a++) {
                *iajy_ptr++ = *gamma_iajs_ptr++;
                }
              }

            }     // endif
          ij_index++;
          }       // endif
        }         // exit j loop
      }           // exit i loop
    // end of first quarter back-transformation
    tim.exit("1. q.b.t.");

    // debug print
    if (debug_ && me == 0) {
      ExEnv::out0() << indent << "End of first q.b.t." << endl;
      }
    // end of debug print

    mo_int = 0;

    mem->sync(); // Make sure all nodes are done with gamma_iajs_tmp before renaming

    delete[] gamma_iajs_tmp;

    // The array mo_int has now been overwritten by the quarter
    // back-transformed non-sep 2PDM gamma_iajs, so rename
    gamma_iajs = (double*) mem->localdata();

    gamma_iqjs_tmp = new double[nbasis];
    if (!gamma_iqjs_tmp) {
      ExEnv::errn() << "Could not allocate gamma_iqjs_tmp" << endl;
      abort();
      }

    if (debug_ && me == 0) {
      ExEnv::out0() << indent << "Begin second q.b.t." << endl;
      }

    // Begin second quarter back-transformation
    // (gamma_iqjs elements ordered as i,j,s,q,
    // i.e., q varies fastest)
    tim.enter("2. q.b.t.");
    index = 0;
    ij_index = 0;
    for (i=0; i<ni; i++) {
      for (j=0; j<nocc; j++) {
        if (index++ % nproc == me) {
          if (j>=nfzc) {
            offset = nbasis*nbasis*ij_index;

            for (s=0; s<nbasis; s++) {
              bzerofast(gamma_iqjs_tmp,nbasis);
              for (q=0; q<nbasis; q++) {
                gamma_iqjs_ptr = &gamma_iqjs_tmp[q];
                gamma_iajs_ptr = &gamma_iajs[nocc + s*nbasis + offset];
                c_qa = &scf_vector[q][nocc];

                tmpval = 0.0;
                for (a=0; a<nvir_act; a++) {
                  tmpval += *c_qa++ * *gamma_iajs_ptr++;
                  } // exit a loop
                *gamma_iqjs_ptr += tmpval;
                }   // exit q loop
              // Put gamma_iqjs_tmp into gamma_iajs for one i,j,s
              // while overwriting gamma_iajs
              gamma_iajs_ptr = &gamma_iajs[s*nbasis + offset];
              gamma_iqjs_ptr = gamma_iqjs_tmp;
              for (q=0; q<nbasis; q++) {
                *gamma_iajs_ptr++ = *gamma_iqjs_ptr++;
                }
              }   // exit s loop

            }     // endif
          ij_index++;
          }       // endif
        }         // exit j loop
      }           // exit i loop
    tim.exit("2. q.b.t.");
    // end of second quarter back-transformation

    if (debug_ && me == 0) {
      ExEnv::out0() << indent << "End of second q.b.t." << endl;
      }

    gamma_iajs = 0;

    mem->sync(); // Keep this here to make sure all nodes have gamma_iqjs
                 // before it is needed below, and that gamma_iajs is not
                 // deleted prematurely

    // The quarter back-transformed elements gamma_iajs have now been
    // overwritten by the half back-transformed elements gamma_iqjs

    delete[] gamma_iqjs_tmp;

    /////////////////////////////////////////////////
    // End of 1. and 2. quarter back-transformation
    /////////////////////////////////////////////////

    Lpi = new double[nbasis*ni];
    bzerofast(Lpi,nbasis*ni);

    if (me == 0) {
      ExEnv::out0() << indent << "Begin third and fourth q.b.t." << endl;
      }

    //////////////////////////////////////////////////////////
    // Perform third and fourth quarter back-transformation
    // and compute contribution to gradient from non-sep 2PDM
    //////////////////////////////////////////////////////////

    tim.enter("3.qbt+4.qbt+non-sep contrib.");
    sp_g_data.init();
    for (i=0; i<thr_->nthread(); i++) {
      qbt34thread[i]->set_i_offset(i_offset);
      qbt34thread[i]->set_ni(ni);
      thr_->add_thread(i,qbt34thread[i]);
#     if SINGLE_THREAD_QBT34
      qbt34thread[i]->run();
#     endif
      }
#   if !SINGLE_THREAD_QBT34
    thr_->start_threads();
    thr_->wait_threads();
#   endif
    tim.exit("3.qbt+4.qbt+non-sep contrib.");
    // Add thread contributions to Lpi and ginter
    for (i=0; i<thr_->nthread(); i++) {
      double *Lpi_thread = qbt34thread[i]->get_Lpi();
      double **ginter_thread = qbt34thread[i]->get_ginter();
      for (j=0; j<nbasis*ni; j++) Lpi[j] += Lpi_thread[j];
      for (j=0; j<natom; j++) {
        for (k=0; k<3; k++) {
          ginter[j][k] += ginter_thread[j][k];
          }
        }
      aointder_computed += qbt34thread[i]->get_aointder_computed();
      }

    if (me == 0) {
      ExEnv::out0() << indent << "End of third and fourth q.b.t." << endl;
      }

    mem->sync(); // Make sure all nodes are done before deleting arrays

    if (debug_ > 1) {
      RefSCDimension ni_dim(new SCDimension(ni,1));
      ni_dim->blocks()->set_subdim(0, new SCDimension(ni));
      RefSCDimension nbasis_dim(new SCDimension(nbasis,1));
      nbasis_dim->blocks()->set_subdim(0, new SCDimension(nbasis));
      RefSCMatrix Lpi_mat(nbasis_dim, ni_dim, kit);
      Lpi_mat->assign(Lpi);
      Lpi_mat.print("Lpi");
      }

    if (debug_ && me == 0) {
      ExEnv::out0() << indent << "Back-transform Lpi" << endl;
      }

    // Back-transform Lpi to MO basis
    lpi_ptr = Lpi;
    for (p=0; p<nbasis; p++) {
      for (i=0; i<ni; i++) {
        c_pa = &scf_vector[p][nocc];
        laj_ptr = &Laj[nvir*(i_offset + i)];
        for (a=0; a<nvir; a++) {
          *laj_ptr++ += *c_pa++ * *lpi_ptr;
          } // exit a loop
        lpi_ptr++;
        }   // exit i loop
      }     // exit p loop

//  malloc_chain_check(1);

    delete[] Lpi;

    if (me == 0) {
      ExEnv::out0() << indent << "Done with pass " << pass+1 << endl;
      }
    }           // exit loop over i-batches (pass)
  tim.exit("mp2 passes");

  if (dograd || do_d1_) {
    for (i=0; i<thr_->nthread(); i++) {
      delete qbt34thread[i];
    }
    delete[] qbt34thread;
  }

  mem->set_localsize(0);

  // debug print
  if (debug_ && me == 0) {
    ExEnv::out0() << indent << "Exited loop over i-batches" << endl;
    }
  // end of debug print

  ///////////////////////////////////////////////////////////////
  // The computation of the MP2 energy is now complete on each
  // node; add the nodes' contributions and print out the energy
  ///////////////////////////////////////////////////////////////
  msg_->sum(ecorr_mp2);
  msg_->sum(aoint_computed);
  msg_->sum(aointder_computed);

  biggest_coefs.combine(msg_);
#if PRINT_BIGGEST_INTS
  biggest_ints_1.combine(msg_);
  biggest_ints_2.combine(msg_);
  biggest_ints_2s.combine(msg_);
  biggest_ints_3a.combine(msg_);
  biggest_ints_3.combine(msg_);
#endif

  if (me == 0) {
    emp2 = escf + ecorr_mp2;

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
    ExEnv::outn() << "biggest 2/4 transformed ints" << endl;
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
    ExEnv::outn() << "restricted 2/4 transformed ints" << endl;
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
    ExEnv::outn() << "biggest 3/4 transformed ints (in 3.)" << endl;
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
    ExEnv::outn() << "biggest 3/4 transformed ints (in 4.)" << endl;
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

    if (restart_orbital_memgrp_ == 0) {
      if (biggest_coefs.ncontrib()) {
        ExEnv::out0() << endl << indent
             << "Largest first order coefficients (unique):"
             << endl;
        }
      for (i=0; i<biggest_coefs.ncontrib(); i++) {
        int i0 = biggest_coefs.indices(i)[0];
        int i1 = biggest_coefs.indices(i)[1];
        int i2 = biggest_coefs.indices(i)[2] + nocc;
        int i3 = biggest_coefs.indices(i)[3] + nocc;
        int spincase = biggest_coefs.indices(i)[4];
        ExEnv::out0() << indent
             << scprintf("  %2d %12.8f %2d %3s %2d %3s -> %2d %3s %2d %3s (%s)",
                         i+1, biggest_coefs.val(i),
                         symorb_num_[i0]+1,
                         ct.gamma(symorb_irrep_[i0]).symbol(),
                         symorb_num_[i1]+1,
                         ct.gamma(symorb_irrep_[i1]).symbol(),
                         symorb_num_[i2]+1,
                         ct.gamma(symorb_irrep_[i2]).symbol(),
                         symorb_num_[i3]+1,
                         ct.gamma(symorb_irrep_[i3]).symbol(),
                         (spincase==1111?"++++":"+-+-")
               )
             << endl;
        }
      }

    // Print out various energies etc.

    if (debug_) {
      ExEnv::out0() << indent << "Number of shell quartets for which AO integrals\n"
           << indent << "(or integral derivatives) would have been computed\n"
           << indent << "without bounds checking: "
           << npass*nshell*nshell*(nshell+1)*(nshell+1)/2
           << endl;

      ExEnv::out0() << indent << "Number of shell quartets for which AO integrals\n"
           << indent << "were computed: " << aoint_computed
           << endl;

      if (dograd) {
        ExEnv::out0() << indent
             << "Number of shell quartets for which AO integral derivatives\n"
             << indent << "were computed: " << aointder_computed
             << endl;
        }
      }

    ExEnv::out0()<<endl<<indent
        <<scprintf("RHF energy [au]:                   %17.12lf\n", escf);
    ExEnv::out0()<<indent
        <<scprintf("MP2 correlation energy [au]:       %17.12lf\n", ecorr_mp2);
    ExEnv::out0()<<indent
        <<scprintf("MP2 energy [au]:                   %17.12lf\n", emp2);
    ExEnv::out0().flush();
    }
  if (not (method_ == "mp" || method_ == "default")) {
    ExEnv::out0() << indent
         << "MBPT2: bad method for closed shell case: " << method_
         << ", using mp" << endl;
    }
  set_energy(emp2);
  set_actual_value_accuracy(reference_->actual_value_accuracy()
                            *ref_to_mp2_acc());

  RefSCDimension nocc_act_dim(new SCDimension(nocc_act,1));
  nocc_act_dim->blocks()->set_subdim(0, new SCDimension(nocc_act));
  RefSCDimension nvir_act_dim(new SCDimension(nvir_act,1));
  nvir_act_dim->blocks()->set_subdim(0, new SCDimension(nvir_act));
  RefSCDimension nocc_dim(new SCDimension(nocc,1));
  nocc_dim->blocks()->set_subdim(0, new SCDimension(nocc));
  RefSCDimension nvir_dim(new SCDimension(nvir,1));
  nvir_dim->blocks()->set_subdim(0, new SCDimension(nvir));
  RefSCDimension nbasis_dim = ao_dimension()->blocks()->subdim(0);
  RefSCDimension noso_dim(new SCDimension(noso,1));

  if (dograd || do_d1_) {
    msg_->sum(Laj,nvir*nocc);

    RefSCMatrix T1_mat(nocc_act_dim, nvir_act_dim, kit);
    // the elements of T1_mat are the single-substitution amplitudes
    BiggestContribs biggest_t1(2,10);
    // compute the S2 norm of Lee et al. (s2_diag)
    double s2_diag = 0.0;
    for (j=nfzc; j<nocc; j++) {
      laj_ptr = &Laj[j*nvir];
      for (a=0; a<nvir_act; a++) {
        tmpval = 0.5**laj_ptr++/(evals[a+nocc]-evals[j]);
        T1_mat.set_element(j-nfzc,a,tmpval);
        biggest_t1.insert(tmpval,j,a);
        s2_diag += tmpval*tmpval;
        }
      }
    s2_diag = sqrt(s2_diag/(2*nocc_act));
    // compute the T1 matrix 1-norm
    double t1onenorm = 0.0;
    Ref<SCElementKNorm> sumabs = new SCElementKNorm(1);
    Ref<SCElementOp> genop = sumabs.pointer();
    for (a=0; a < nvir_act; a++) {
      sumabs->init();
      T1_mat.get_column(a).element_op(genop);
      if (t1onenorm < sumabs->result()) t1onenorm = sumabs->result();
      }
    // compute the T1 matrix inf-norm
    double t1infnorm = 0.0;
    for (j=0; j < nocc_act; j++) {
      sumabs->init();
      T1_mat.get_row(j).element_op(genop);
      if (t1infnorm < sumabs->result()) t1infnorm = sumabs->result();
      }
    // compute the T1 matrix 2-norm ( = D1(MP2) )
    RefSymmSCMatrix D1_mat(nocc_act_dim,kit);
    D1_mat.assign(0.0);
    D1_mat.accumulate_symmetric_product(T1_mat);
    T1_mat = 0;
    double d1_diag = sqrt(D1_mat.eigvals().get_element(nocc_act-1));
    D1_mat = 0;
    // print the norms
    ExEnv::out0()<<endl;
    ExEnv::out0()
         <<indent<<scprintf("D1(MP2)                = %12.8f", d1_diag)
         <<endl
         <<indent<<scprintf("S2 matrix 1-norm       = %12.8f", t1onenorm)
         <<endl
         <<indent<<scprintf("S2 matrix inf-norm     = %12.8f", t1infnorm)
         <<endl
         <<indent<<scprintf("S2 diagnostic          = %12.8f", s2_diag)
         <<endl;
    if (biggest_t1.ncontrib()) {
      ExEnv::out0() << endl
           << indent << "Largest S2 values (unique determinants):" << endl;
      }
    for (i=0; i<biggest_t1.ncontrib(); i++) {
      int i0 = biggest_t1.indices(i)[0];
      int i1 = biggest_t1.indices(i)[1] + nocc;
      ExEnv::out0()
           << indent << scprintf("  %2d %12.8f  %3d %3s -> %4d %3s",
                                 i+1, biggest_t1.val(i),
                                 symorb_num_[i0]+1,
                                 ct.gamma(symorb_irrep_[i0]).symbol(),
                                 symorb_num_[i1]+1,
                                 ct.gamma(symorb_irrep_[i1]).symbol())
           << endl;
      }

    } // if (dograd || do_d1_)

  for (i=0; i<thr_->nthread(); i++) {
      delete e12thread[i];
    }
  delete[] e12thread;

  // quit here if only the energy is needed
  if (!dograd) {
    delete[] tbints_; tbints_ = 0;
    if (do_d1_) free(Laj);
    delete[] scf_vector;
    delete[] scf_vector_dat;
    delete[] evals;
    tim.exit("mp2-mem");
    return;
    }

  // Accumulate intermediate gradients on node 0
  sum_gradients(msg_, ginter, natom, 3);

  // Add intermediate gradients to the gradient on node 0
  if (me==0)
    accum_gradients(gradient, ginter, natom, 3);

  // Print out contribution to the gradient from non-sep. 2PDM
  if (debug_) {
    print_natom_3(ginter,
              "Contribution to MP2 gradient from non-separable 2PDM [au]:");
    }

  ////////////////////////////////////////////////////////
  // Add contributions from all nodes to various matrices
  ////////////////////////////////////////////////////////
  tmpint = (nvir > nocc ? nvir:nocc);
  double *tmpmat = new double[tmpint*tmpint];
  msg_->sum(Pkj,nocc*(nocc+1)/2,tmpmat); // Pkj is now complete
  msg_->sum(Pab,nvir*(nvir+1)/2,tmpmat); // Pab is now complete
  msg_->sum(Wab,nvir*nvir,tmpmat);
  msg_->sum(Wkj,nocc*nocc,tmpmat);
  msg_->sum(Waj,nvir*nocc,tmpmat);
  if (do_d2_) {
    msg_->sum(d2occ_mat,nocc_act*(nocc_act+1)/2,tmpmat);
    msg_->sum(d2vir_mat,nvir_act*(nvir_act+1)/2,tmpmat);
    }
  delete[] tmpmat;

  // Compute D2 diagnostic (d2_diag) from matrices d2_occ_mat and d2_vir_mat
  if (do_d2_) {
    RefSymmSCMatrix D2occ_mat(nocc_act_dim, kit);
    RefSymmSCMatrix D2vir_mat(nvir_act_dim,kit);
    D2occ_mat->assign(d2occ_mat);
    D2vir_mat->assign(d2vir_mat);
    d2o = sqrt(D2occ_mat.eigvals().get_element(nocc_act-1));
    d2v = sqrt(D2vir_mat.eigvals().get_element(nvir_act-1));
    d2_diag = (d2o > d2v ? d2o:d2v);
    ExEnv::out0() << endl
         << indent <<  scprintf("D2(MP1) = %12.8f", d2_diag) << endl << endl;
    delete[] d2occ_mat;
    delete[] d2vir_mat;
  }

  // Finish computation of Wab
  tim.enter("Pab and Wab");
  pab_ptr = Pab;
  for (a=0; a<nvir_act; a++) {  // active-active part of Wab
    wba_ptr = &Wab[a];
    wab_ptr = &Wab[a*nvir];
    for (b=0; b<=a; b++) {
      if (a==b) {
        *wab_ptr -= evals[nocc+a]**pab_ptr;
        }
      else {
        *wab_ptr -= evals[nocc+a]**pab_ptr;
        *wba_ptr -= evals[nocc+b]**pab_ptr;
        }
      pab_ptr++;
      wab_ptr++;
      wba_ptr += nvir;
      } // exit b loop
    }   // exit a loop
  for (a=0; a<nfzv; a++) {  // active-frozen part of Wab
    wba_ptr = &Wab[nvir_act+a];
    wab_ptr = &Wab[(nvir_act+a)*nvir];
    pab_ptr = &Pab[(nvir_act+a)*(nvir_act+a+1)/2];
    for (b=0; b<nvir_act; b++) {
      *wab_ptr -= evals[nocc+b]**pab_ptr;
      *wba_ptr -= evals[nocc+b]**pab_ptr;
      pab_ptr++;
      wab_ptr++;
      wba_ptr += nvir;
      } // exit b loop
    }   // exit a loop
  // Wab is now complete
  tim.exit("Pab and Wab");
  RefSCMatrix Wab_matrix(nvir_dim, nvir_dim, kit);
  Wab_matrix->assign(Wab); // Put elements of Wab into Wab_matrix
  free(Wab);

  // Update Wkj with contribution from Pkj
  tim.enter("Pkj and Wkj");
  pkj_ptr = Pkj;
  for (k=0; k<nocc; k++) {
    wjk_ptr = &Wkj[k];
    wkj_ptr = &Wkj[k*nocc];
    for (j=0; j<=k; j++) {
      if (k<nfzc && j<nfzc) {   // don't want both j and k frozen
        wkj_ptr++;
        wjk_ptr += nocc;
        pkj_ptr++;
        continue;
        }
      if (j==k) {
        *wkj_ptr++ -= evals[k]**pkj_ptr++;
        }
      else if (j<nfzc) {
        *wkj_ptr++ -= evals[k]**pkj_ptr;
        *wjk_ptr   -= evals[k]**pkj_ptr;
        pkj_ptr++;
        }
      else {
        *wkj_ptr++ -= evals[k]**pkj_ptr;
        *wjk_ptr   -= evals[j]**pkj_ptr;
        pkj_ptr++;
        }
      wjk_ptr += nocc;
      }  // exit j loop
    }    // exit k loop
  tim.exit("Pkj and Wkj");

  /////////////////////////////////
  // Finish the computation of Laj
  /////////////////////////////////

  tim.enter("Laj");

  RefSCMatrix Cv(nbasis_dim, nvir_dim, kit); // virtual block of scf_vector
  RefSCMatrix Co(nbasis_dim, nocc_dim, kit); // occupied block of scf_vector
  for (p=0; p<nbasis; p++) {
    c_pq = scf_vector[p];
    for (q=0; q<noso; q++) {
      if (q<nocc) Co->set_element(p, q, *c_pq++);
      else Cv->set_element(p, q-nocc, *c_pq++);
      }
    }

  // Compute the density-like Dmat_matrix
  // (Cv*Pab_matrix*Cv.t() + Co*Pkj_matrix*Co.t())
  RefSymmSCMatrix Pab_matrix(nvir_dim,kit);
  RefSymmSCMatrix Pkj_matrix(nocc_dim,kit);
  RefSymmSCMatrix Dmat_matrix(nbasis_dim,kit);
  Pab_matrix->assign(Pab); // fill in elements of Pab_matrix from Pab
  free(Pab);
  Pkj_matrix->assign(Pkj); // fill in elements of Pkj_matrix from Pkj
  free(Pkj);
  Dmat_matrix.assign(0.0);
  Dmat_matrix.accumulate_transform(Cv,Pab_matrix);
  Dmat_matrix.accumulate_transform(Co,Pkj_matrix);
  // We now have the density-like matrix Dmat_matrix

  // Compute the G matrix

  RefSymmSCMatrix Gmat(nbasis_dim,kit);
  init_cs_gmat();
  tim.enter("make_gmat for Laj");
  make_cs_gmat_new(Gmat, Dmat_matrix);
  if (debug_ > 1) {
    Dmat_matrix.print("Dmat");
    Gmat.print("Gmat");
    }
  tim.exit("make_gmat for Laj");

  // Finish computation of Laj
  RefSCMatrix Laj_matrix(nocc_dim,nvir_dim,kit); // elements are ordered as j*nvir+a
  Laj_matrix->assign(Laj);
  if (debug_ > 1) Laj_matrix->print("Laj (first bit)");
  Laj_matrix = Laj_matrix - 2.0*Co.t()*Gmat*Cv;
  if (debug_ > 1) Laj_matrix->print("Laj (all of it)");
  Laj_matrix->convert(Laj);  // Put new Laj_matrix elements into Laj

  tim.exit("Laj");

  //////////////////////////////////////
  // Computation of Laj is now complete
  //////////////////////////////////////

  ////////////////////////////
  // Solve the CPHF equations
  ////////////////////////////
  RefSCMatrix Paj_matrix(nvir_dim, nocc_dim, kit);
  tim.enter("cphf");
  cs_cphf(scf_vector, Laj, evals, Paj_matrix);
  tim.exit("cphf");

  free(Laj);

  // Finish computation of Waj
  for (a=0; a<nvir; a++) {
    waj_ptr = &Waj[a];
    for (j=0; j<nocc; j++) {
      *waj_ptr -= evals[j]*Paj_matrix->get_element(a,j);
      waj_ptr += nvir;
      }
    }
  // Waj is now complete
  RefSCMatrix Waj_matrix(nocc_dim, nvir_dim, kit);
  Waj_matrix->assign(Waj); // Put elements of Waj into Waj_matrix
  // NB. Waj_matrix elements are ordered as j*nvir+a
  free(Waj);


  // Finish computation of Wkj
  tim.enter("Pkj and Wkj");
  // Compute Dmat_matrix =
  // Co*Pkj_matrix*Co.t() + Co*Paj_matrix.t()*Cv.t()
  // + Cv*Paj_matrix*Co.t() + Cv*Pab_matrix*Cv.t();
  Dmat_matrix.assign(0.0);
  Dmat_matrix.accumulate_symmetric_sum(Cv*Paj_matrix*Co.t());
  Dmat_matrix.accumulate_transform(Co,Pkj_matrix);
  Dmat_matrix.accumulate_transform(Cv,Pab_matrix);
  tim.enter("make_gmat for Wkj");
  make_cs_gmat_new(Gmat, Dmat_matrix);
  tim.exit("make_gmat for Wkj");
  done_cs_gmat();
  for (i=0; i<thr_->nthread(); i++) tbints_[i] = 0;
  delete[] tbints_; tbints_ = 0;
  RefSCMatrix Wkj_matrix(nocc_dim, nocc_dim, kit);
  Wkj_matrix->assign(Wkj);
  Wkj_matrix = Wkj_matrix - 2.0*Co.t()*Gmat*Co;
  free(Wkj);
  // Wkj is now complete - not as Wkj but as Wkj_matrix
  tim.exit("Pkj and Wkj");

  ////////////////////////////////////////////////////////////////
  // We now have the matrices Pkj_matrix, Paj_matrix, Pab_matrix,
  // Wkj_matrix, Waj_matrix, Wab_matrix and can compute the
  // remaining contributions to the gradient
  ///////////////////////////////////////////////////////////////

  // Compute the second order correction to
  // the density matrix and energy weighted
  // density matrix in the AO basis
  RefSCMatrix P2AO_matrix(nbasis_dim, nbasis_dim, kit);
  RefSCMatrix P2MO_matrix(noso_dim, noso_dim, kit);
  RefSCMatrix W2AO_matrix(nbasis_dim, nbasis_dim, kit);
  RefSCMatrix W2MO_matrix(noso_dim, noso_dim, kit);
  RefSCMatrix SCF_matrix(nbasis_dim, noso_dim, kit);
  for (i=0; i<nocc; i++) {
    for (j=0; j<nocc; j++) {
      P2MO_matrix->set_element(i,j,Pkj_matrix->get_element(i,j));
      W2MO_matrix->set_element(i,j,Wkj_matrix->get_element(i,j));
      }
    for (j=nocc; j<noso; j++) {
      P2MO_matrix->set_element(i,j,Paj_matrix->get_element(j-nocc,i));
      W2MO_matrix->set_element(i,j,Waj_matrix->get_element(i,j-nocc));
      }
    }
  for (i=nocc; i<noso; i++) {
    for (j=0; j<nocc; j++) {
      P2MO_matrix->set_element(i,j,Paj_matrix->get_element(i-nocc,j));
      W2MO_matrix->set_element(i,j,Waj_matrix->get_element(j,i-nocc));
      }
    for (j=nocc; j<noso; j++) {
      P2MO_matrix->set_element(i,j,Pab_matrix->get_element(i-nocc,j-nocc));
      W2MO_matrix->set_element(i,j,Wab_matrix->get_element(i-nocc,j-nocc));
      }
    }
  for (i=0; i<nbasis; i++) {
    for (j=0; j<nocc; j++) {
      SCF_matrix->set_element(i,j,Co->get_element(i,j));
      }
    for (j=nocc; j<noso; j++) {
      SCF_matrix->set_element(i,j,Cv->get_element(i,j-nocc));
      }
    }
  P2AO_matrix = SCF_matrix * P2MO_matrix * SCF_matrix.t();
  W2AO_matrix = SCF_matrix * W2MO_matrix * SCF_matrix.t();
//  P2AO_matrix = Co*(Pkj_matrix*Co.t() + Paj_matrix.t()*Cv.t()) +
//                Cv*(Paj_matrix*Co.t() + Pab_matrix*Cv.t());
//  W2AO_matrix = Co*(Wkj_matrix*Co.t() + Waj_matrix*Cv.t()) +
//                Cv*(Waj_matrix.t()*Co.t() + Wab_matrix*Cv.t());

  if (debug_ > 1) {
    SCF_matrix.print("SCF_matrix");
    P2MO_matrix.print("P2MO_matrix");
    W2MO_matrix.print("W2MO_matrix");
    P2AO_matrix.print("P2AO_matrix");
    W2AO_matrix.print("W2AO_matrix");
    }

  // Convert P2AO_matrix and W2AO_matrix to double*
  P2AO = new double[nbasis*nbasis];
  W2AO = new double[nbasis*nbasis];
  P2AO_matrix->convert(P2AO);
  W2AO_matrix->convert(W2AO);

  // Compute the HF density matrix and
  // energy weighted density matrix
  PHF = new double[nbasis*nbasis];
  WHF = new double[nbasis*nbasis];
  phf_ptr = PHF;
  whf_ptr = WHF;
  for (p=0; p<nbasis; p++) {
    for (q=0; q<nbasis; q++) {
      *phf_ptr++ = 0.0;
      *whf_ptr++ = 0.0;
      }
    }
  phf_ptr = PHF;
  whf_ptr = WHF;
  for (p=0; p<nbasis; p++) {
    for (q=0; q<nbasis; q++) {
      c_pi = scf_vector[p];
      c_qi = scf_vector[q];
      for (i=0; i<nocc; i++) {
        tmpval = 2* *c_pi++ * *c_qi++;
        *phf_ptr += tmpval;
        *whf_ptr -= evals[i] * tmpval;
        } // exit i loop
      phf_ptr++;
      whf_ptr++;
      }   // exit q loop
    }     // exit p loop

  // Compute the MP2 density and energy weighted density

  // PMP2 = PHF + P2AO; WMP2 = WHF + W2AO
  PMP2 = new double[nbasis*nbasis];
  WMP2 = new double[nbasis*nbasis];
  // Initialize PMP2 and WMP2
  pmp2_ptr = PMP2;
  wmp2_ptr = WMP2;
  for (p=0; p<nbasis; p++) {
    for (q=0; q<nbasis; q++) {
      *pmp2_ptr++ = 0.0;
      *wmp2_ptr++ = 0.0;
      }
    }
  pmp2_ptr = PMP2;
  wmp2_ptr = WMP2;
  p2ao_ptr = P2AO;
  w2ao_ptr = W2AO;
  phf_ptr = PHF;
  whf_ptr = WHF;
  for (p=0; p<nbasis; p++) {
    for (q=0; q<nbasis; q++) {
      *pmp2_ptr++ = *phf_ptr++ + *p2ao_ptr++;
      *wmp2_ptr++ = *whf_ptr++ + *w2ao_ptr++;
      }
    }
  delete[] W2AO;

  if (debug_ > 1) {
    RefSCMatrix tmpmat(ao_dimension(), ao_dimension(), kit);
    tmpmat->assign(PMP2);
    tmpmat.print("PMP2");
    tmpmat->assign(P2AO);
    tmpmat.print("P2AO");
    tmpmat->assign(PHF);
    tmpmat.print("PHF");
    tmpmat->assign(WMP2);
    tmpmat.print("WMP2");
    }

  ////////////////////////////////////////////////
  // Compute the contribution to the MP2 gradient
  // from the separable part of the 2PDM
  ////////////////////////////////////////////////

  zero_gradients(ginter, natom, 3);
  zero_gradients(hf_ginter, natom, 3);
  tim.enter("sep 2PDM contrib.");

  CSGradS2PDM** s2pdmthread = new CSGradS2PDM*[thr_->nthread()];
  for (i=0; i<thr_->nthread(); i++) {
      s2pdmthread[i] = new CSGradS2PDM(i, thr_->nthread(), me, nproc,
                                       lock, basis(), tbintder_[i],
                                       PHF, P2AO, tol, debug_, dynamic_);
      thr_->add_thread(i,s2pdmthread[i]);
#     if SINGLE_THREAD_S2PDM
      s2pdmthread[i]->run();
#     endif
    }
# if !SINGLE_THREAD_S2PDM
  thr_->start_threads();
  thr_->wait_threads();
# endif
  for (i=0; i<thr_->nthread(); i++) {
      s2pdmthread[i]->accum_mp2_contrib(ginter);
      s2pdmthread[i]->accum_hf_contrib(hf_ginter);
      delete s2pdmthread[i];
    }
  sum_gradients(msg_, ginter, molecule()->natom(), 3);
  sum_gradients(msg_, hf_ginter, molecule()->natom(), 3);
  delete[] s2pdmthread;

  tim.exit("sep 2PDM contrib.");
  delete[] P2AO;

  // The separable 2PDM contribution to the gradient has now been
  // accumulated in ginter on node 0; add it to the total gradients
  if (me == 0) {
    accum_gradients(gradient, ginter, natom, 3);
    accum_gradients(hf_gradient, hf_ginter, natom, 3);
    }
  // Print out the contribution to the gradient from sep. 2PDM
  if (debug_) {
    print_natom_3(hf_ginter,
                  "Contribution from separable 2PDM to HF gradient [au]:");
    print_natom_3(ginter,
                  "Contribution from separable 2PDM to MP2 gradient [au]:");
    }

  // Done with two-electron integrals
  tbint_ = 0;
  if (dograd || do_d1_) {
    delete[] tbintder_;
    tbintder_ = 0;
    }

  /////////////////////////////////////////////////////////////
  // Compute the one-electron contribution to the MP2 gradient
  /////////////////////////////////////////////////////////////

  zero_gradients(ginter, natom, 3);
  zero_gradients(hf_ginter, natom, 3);
  tim.enter("hcore contrib.");
  hcore_cs_grad(PHF, PMP2, hf_ginter, ginter);
  tim.exit("hcore contrib.");
  delete[] PHF;
  delete[] PMP2;
  // The hcore contribution to the gradient has now been accumulated
  // in ginter on node 0; add it to the total gradients
  if (me == 0) {
    accum_gradients(gradient, ginter, natom, 3);
    accum_gradients(hf_gradient, hf_ginter, natom, 3);
    }
  // Print out the contribution to the gradient from hcore
  if (debug_) {
    print_natom_3(hf_ginter, "Contribution to HF gradient from hcore [au]:");
    print_natom_3(ginter, "Contribution to MP2 gradient from hcore [au]:");
    }

  zero_gradients(ginter, natom, 3);
  zero_gradients(hf_ginter, natom, 3);
  tim.enter("overlap contrib.");
  overlap_cs_grad(WHF, WMP2, hf_ginter, ginter);
  delete[] WHF;
  tim.exit("overlap contrib.");
  delete[] WMP2;
  // The overlap contribution to the gradient has now been accumulated
  // in ginter on node 0; add it to the total gradients
  if (me == 0) {
      accum_gradients(gradient, ginter, natom, 3);
      accum_gradients(hf_gradient, hf_ginter, natom, 3);
    }

  // Print out the overlap contribution to the gradient
  if (debug_) {
    print_natom_3(hf_ginter, "Overlap contribution to HF gradient [au]:");
    print_natom_3(ginter,"Overlap contribution to MP2 gradient [au]:");
    }

  ////////////////////////////////////////////////////////
  // Compute the nuclear contribution to the MP2 gradient
  ////////////////////////////////////////////////////////

  if (me == 0) {
    nuclear_repulsion_energy_gradient(ginter);
    accum_gradients(gradient, ginter, natom, 3);
    accum_gradients(hf_gradient, ginter, natom, 3);

    }

  // Print out the nuclear contribution to the gradient
  if (debug_)
    print_natom_3(ginter,"Nuclear contribution to MP2 gradient [au]:");


  ////////////////////////////////////////////////////////
  // The computation of the MP2 gradient is now complete;
  // print out the gradient
  ////////////////////////////////////////////////////////
  if (debug_) {
    ExEnv::out0() << "Obtaining HF gradient" << endl;
    print_natom_3(ref()->gradient(),"HF gradient from HF");
    print_natom_3(hf_gradient,"Total HF gradient from MP2 [au]:");
    }
  print_natom_3(gradient,"Total MP2 gradient [au]:");

  msg_->bcast(gradient_dat, natom*3);
  RefSCVector gradientvec = matrixkit()->vector(moldim());
  gradientvec->assign(gradient_dat);
  set_gradient(gradientvec);

  msg_->bcast(hf_gradient_dat, natom*3);
  hf_gradient_ = matrixkit()->vector(moldim());
  hf_gradient_->assign(hf_gradient_dat);

  delete[] gradient;
  delete[] gradient_dat;

  delete[] hf_gradient;
  delete[] hf_gradient_dat;

  for (i=0; i<natom; i++) {
    delete[] ginter[i];
    delete[] hf_ginter[i];
    }
  delete[] ginter;
  delete[] hf_ginter;

  delete[] scf_vector;
  delete[] scf_vector_dat;
  delete[] evals;

  tim.exit("mp2-mem");
  }

///////////////////////////////////////////////////////////
// Compute the contribution to the MP2 gradient from hcore
///////////////////////////////////////////////////////////
void
MBPT2::hcore_cs_grad(double *PHF, double *PMP2,
                     double **hf_ginter, double **ginter)
{

  int i, j, k, l, m;
  int jj, kk;
  int jsize, ksize;
  int j_offset, k_offset;
  int jk_index;
  int index;
  int nshell;
  int nbasis;
  int nproc = msg_->n();
  int me = msg_->me();

  const double *oneebuf; // 1-electron buffer
  double tmpval1, tmpval2;
  double gxyz[3];
  double hf_gxyz[3];

  // Initialize 1e object
  Ref<OneBodyDerivInt> obintder_ = integral()->hcore_deriv();
  oneebuf = obintder_->buffer();

  nshell = basis()->nshell();
  nbasis = basis()->nbasis();

  ///////////////////////////////////////////////////////////////////////////////
  // Compute the kinetic and nuclear-electron energy contribution to the gradient
  ///////////////////////////////////////////////////////////////////////////////

  jk_index = 0;

  for (i=0; i<molecule()->natom(); i++) {
    for (j=0; j<nshell; j++) {
      jsize = basis()->shell(j).nfunction();
      j_offset = basis()->shell_to_function(j);

      for (k=0; k<=j; k++) {
        ksize = basis()->shell(k).nfunction();
        k_offset = basis()->shell_to_function(k);

        if (jk_index++%nproc == me) {
          obintder_->compute_shell(j,k,i);

          for (l=0; l<3; l++) { gxyz[l] = 0.0; hf_gxyz[l] = 0.0; }

          index = 0;

          for (jj=0; jj<jsize; jj++) {
            for (kk=0; kk<ksize; kk++) {
              tmpval1 = PMP2[(j_offset + jj)*nbasis + k_offset + kk];
              tmpval2 = PHF[(j_offset + jj)*nbasis + k_offset + kk];
              for (m=0; m<3; m++) {
                gxyz[m] += oneebuf[index] * tmpval1;
                hf_gxyz[m] += oneebuf[index] * tmpval2;
                index++;
                } // exit m loop
              }   // exit kk loop
            }     // exit jj loop

          if (j != k) {
            // off-diagonal blocks
            for (l=0; l<3; l++) {
              gxyz[l] *= 2.0;
              hf_gxyz[l] *= 2.0;
              }
            }

          for (l=0; l<3; l++) {
            ginter[i][l] += gxyz[l];
            hf_ginter[i][l] += hf_gxyz[l];
            }

          } // exit "if"
        }   // exit k loop
      }     // exit j loop
    }       // exit i loop

  /* Accumulate the nodes' intermediate gradients on node 0 */
  sum_gradients(msg_, ginter, molecule()->natom(), 3);
  sum_gradients(msg_, hf_ginter, molecule()->natom(), 3);
}


////////////////////////////////////////////////////
// Compute the overlap contribution to the gradient
////////////////////////////////////////////////////
void
MBPT2::overlap_cs_grad(double *WHF, double *WMP2,
                       double **hf_ginter, double **ginter)
{
  int i, j, k, l, m;
  int jj, kk;
  int jj_index, kk_index;
  int jsize, ksize;
  int j_offset, k_offset;
  int jk_index;
  int index;
  int nshell;
  int nbasis;

  const double *oneebuf; // 1-electron buffer
  double tmpval1, tmpval2;
  double hf_tmpval1, hf_tmpval2;
  double gxyz[3];
  double hf_gxyz[3];

  // Initialize 1e object
  Ref<OneBodyDerivInt> obintder_ = integral()->overlap_deriv();
  oneebuf = obintder_->buffer();

  nshell = basis()->nshell();
  nbasis = basis()->nbasis();
  int nproc = msg_->n();
  int me = msg_->me();

  for (i=0; i<molecule()->natom(); i++) {
    jk_index = 0;

    for (j=0; j<nshell; j++) {
      j_offset = basis()->shell_to_function(j);
      jsize = basis()->shell(j).nfunction();

      for (k=0; k<=j; k++) {
        k_offset = basis()->shell_to_function(k);
        ksize = basis()->shell(k).nfunction();

        if (jk_index++%nproc == me) {
          obintder_->compute_shell(j,k,i);

          for (l=0; l<3; l++) {
            hf_gxyz[l] = 0.0;
            gxyz[l] = 0.0;
            }
          index = 0;

          for (jj=0; jj<jsize; jj++) {
            jj_index = j_offset + jj;
            for (kk=0; kk<ksize; kk++) {
              kk_index = k_offset + kk;
              if (kk_index > jj_index) {
                index += 3;  // increment index since m-loop will be skipped
                break;       // skip to next jj value
                }
              // NB. WMP2 is not a symmetrix matrix
              tmpval1 = WMP2[jj_index*nbasis + kk_index];
              tmpval2 = WMP2[kk_index*nbasis + jj_index];
              hf_tmpval1 = WHF[jj_index*nbasis + kk_index];
              hf_tmpval2 = WHF[kk_index*nbasis + jj_index];

              for (m=0; m<3; m++) {
                if (jj_index != kk_index) {
                  gxyz[m] += oneebuf[index] * (tmpval1 + tmpval2);
                  hf_gxyz[m] += oneebuf[index] * (hf_tmpval1 + hf_tmpval2);
                  }
                else {
                  gxyz[m] += oneebuf[index] * tmpval1;
                  hf_gxyz[m] += oneebuf[index] * hf_tmpval1;
                  }
                index++;
                } // exit m loop
              }   // exit kk loop
            }     // exit jj loop

          for (l=0; l<3; l++) {
            ginter[i][l] += gxyz[l];
            hf_ginter[i][l] += hf_gxyz[l];
            }
          } // exit "if"
        }   // exit k loop
      }     // exit j loop
    }       // exit i loop

  /* Accumulate the nodes' intermediate gradients on node 0 */
  sum_gradients(msg_, ginter, molecule()->natom(), 3);
  sum_gradients(msg_, hf_ginter, molecule()->natom(), 3);
}

static void
sum_gradients(const Ref<MessageGrp>& msg, double **f, int n1, int n2)
{
  int i;

  if (msg->n() == 1) return;

  for (i=0; i<n1; i++) {
    msg->sum(f[i],n2);
    }
}

static void
zero_gradients(double **f, int n1, int n2)
{
  for (int i=0; i<n1; i++) {
    for (int j=0; j<3; j++) f[i][j] = 0.0;
    }
}

static void
accum_gradients(double **g, double **f, int n1, int n2)
{
  for (int i=0; i<n1; i++) {
    for (int j=0; j<3; j++) g[i][j] += f[i][j];
    }
}

/////////////////////////////////////
// Compute the batchsize
//
// Only arrays allocated before exiting the loop over
// i-batches are included here  - only these arrays
// affect the batch size.
/////////////////////////////////////
int
MBPT2::compute_cs_batchsize(size_t mem_static, int nocc_act)
{
  size_t mem_dyn;   // dynamic memory available
  distsize_t maxdyn;
  int ni;

  ///////////////////////////////////////
  // the largest memory requirement will
  // either occur just before the end of
  // the 1. q.b.t. (mem1) or just before
  // the end of the i-batch loop (mem2)
  ///////////////////////////////////////

  if (mem_alloc > mem_static) {
    mem_dyn = mem_alloc - mem_static;
    }
  else {
    mem_dyn = 0;
    }

  // First determine if calculation is possible at all (i.e., if ni=1 possible)

  ni = 1;
  maxdyn = compute_cs_dynamic_memory(ni, nocc_act);

  if (maxdyn > mem_dyn) {
    return 0;
    }

  ni = 2;
  while (ni<=nocc_act) {
    maxdyn = compute_cs_dynamic_memory(ni, nocc_act);
    if (maxdyn >= mem_dyn) {
      ni--;
      break;
      }
    ni++;
    }
  if (ni > nocc_act) ni = nocc_act;

  return ni;
}

/////////////////////////////////////
// Compute required (dynamic) memory
// for a given batch size
//
// Only arrays allocated before exiting the loop over
// i-batches are included here  - only these arrays
// affect the batch size.
/////////////////////////////////////
distsize_t
MBPT2::compute_cs_dynamic_memory(int ni, int nocc_act)
{
  int index;
  distsize_t mem1, mem2, mem3;
  distsize_t maxdyn;
  distsize_t tmp;
  int i, j;
  int nij;
  int nproc = msg_->n();

  ///////////////////////////////////////
  // the largest memory requirement will
  // either occur just before the end of
  // the 1. q.b.t. (mem1) or just before
  // the end of the i-batch loop (mem2)
  ///////////////////////////////////////

  // compute nij as nij on node 0, since nij on node 0 is >= nij on other nodes
  index = 0;
  nij = 0;
  for (i=0; i<ni; i++) {
    for (j=0; j<nocc; j++) {
      if (index++ % nproc == 0) nij++;
      }
    }
  mem1 = sizeof(double)*((distsize_t)nij*nbasis*nbasis
                         + nbasis*nvir);
  mem2 = sizeof(double)*((distsize_t)thr_->nthread()*ni*nbasis*nfuncmax*nfuncmax
                         + (distsize_t)nij*nbasis*nbasis
                         + ni*nbasis + nbasis*nfuncmax
                         + 2*nfuncmax*nfuncmax*nfuncmax*nfuncmax);
  mem3 = sizeof(double)*((distsize_t)ni*nbasis*nfuncmax*nfuncmax
                         + (distsize_t)nij*nbasis*nbasis
                         + 2*(2 + nbasis*nfuncmax));
  tmp = (mem2>mem3 ? mem2:mem3);
  maxdyn = (mem1>tmp ? mem1:tmp);

  return maxdyn;
}

////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
