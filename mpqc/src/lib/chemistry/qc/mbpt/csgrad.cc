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

#include <util/misc/formio.h>
#include <util/misc/timer.h>
#include <util/group/memory.h>
#include <util/group/message.h>
#include <util/class/class.h>
#include <util/state/state.h>
#include <math/scmat/matrix.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/mbpt/bzerofast.h>
#include <chemistry/qc/mbpt/mbpt.h>
#include <chemistry/qc/mbpt/util.h>

static void sum_gradients(const RefMessageGrp& msg, double **f, int n1, int n2);
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

  RefSCMatrixKit kit = SCMatrixKit::default_matrixkit();

  int nij;        // number of i,j pairs on a node (for e.g., mo_int)
  double *mo_int; // MO integrals of type (ov|ov)
                  // (and these integrals divided by
                  // orbital energy denominators)
  double *integral_iqjs; // half-transformed integrals
  double *iqjs_buf;

  int nocc_act, nvir_act;
  int i, j, k;
  int ii, bb;
  int x, y;
  int isize, jsize;
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
  int flags;
  int me;
  int nproc;
  int rest;
  int p_offset, q_offset, r_offset, s_offset;

  int aoint_computed = 0; 
  int aointder_computed = 0; 
  int derset, xyz;
  int natom = molecule()->natom();     // the number of atoms
  int int_index;
  int mem_static;    // static memory in bytes
  int qp, sr;
  int factor_pqrs;
  int ij_proc;          // the processor which has ij pair
  int ij_index;         // of the ij pairs on a proc, this ij pair is number ij_index
                        // (i.e., ij_index < nij)
  int ik_proc;          // the processor which has ik pair
  int ik_index;
  int ij_offset;
  int jloop, kloop;
  int ntri;

  int ni;

  double *evals;              // scf eigenvalues
  const double *intbuf;       // 2-electron AO integral buffer
  const double *intderbuf;    // 2-electron AO integral derivative buffer
  double *iajb_ptr, *ibja_ptr, *iakb_ptr, *ibka_ptr;
  double *iajc_ptr, *ibjc_ptr, *icjb_ptr, *icja_ptr;
  double *ijkb_ptr, *ibkj_ptr;
  double pqrs;
  double *c_sa, c_rj;
  double *c_qk, *c_pi, *c_qi, *c_sj;
  double *c_qx, *c_qa, *c_sb, *c_pa, *c_pq, *c_sy;
  double delta_ijab, delta_ijbc, delta_ijac;
  double ecorr_mp2 = 0.0;
  double escf;
  double emp2;
  double tol;                 // log2 of the erep tolerance
                              // (erep < 2^tol => discard)
  double dtol;                // non-log2 version of the above
  double *Wkj,*Wab,*Waj;      // occ-occ, vir-vir and vir-occ parts of 
                              // second order correction to MP2
                              // energy weighted density matrix
  double *Pkj,*Pab;           // occ-occ and vir-vir parts of second order
                              // correction to MP2 density matrix
  double *Laj;                // MP2 Lagrangian
  double *Lpi;                // contrib to MP2 Lagrangian partially in AO basis
  double *pkj_ptr, *pab_ptr;
  double *wkj_ptr, *wjk_ptr, *wab_ptr, *wba_ptr, *waj_ptr;
  double *laj_ptr, *lpi_ptr, *lqi_ptr;
  double *gamma_iajs, *gamma_iajs_tmp, *gamma_iqrs; 
                              // partially back-transformed non-sep 2PDM's
  double *gamma_iqjs_tmp;
  double *gamma_pqrs;
  double *gamma_iajs_ptr;
  double *gamma_iqjs_ptr, *gamma_irjq_ptr;
  double *gamma_iqrs_ptr, *gamma_iprs_ptr;
  double *gamma_iqsr_ptr;
  double *gamma_pqrs_ptr;
  double *gammabuf;           // buffer used for sending elements of gamma_iqjs
  double *mo_intbuf;          // buffer used for sending mo integrals
  double *grad_ptr1, *grad_ptr2;
  double tmpval;
  double *Dmat;
  double *P2AO, *W2AO;
  double *p2ao_ptr, *w2ao_ptr;
  double *PHF, *WHF;
  double *phf_ptr, *whf_ptr;
  double *PMP2, *WMP2;
  double *pmp2_ptr, *wmp2_ptr;

  double *integral_iqrs; // quarter transformed two-el integrals
  double *ixjs_tmp;      // three-quarter transformed two-el integrals
  double *integral_ixjs;  // all three-quarter transformed two-el integrals
  double *integral_iajy; // mo integrals (y = any MO)
  double *integral_ikja; // mo integrals
  double *iqjs_contrib;  // local contributions to integral_iqjs
  double *iqjr_contrib;  // local contributions to integral_iqjr
  double *integral_iqjs_ptr;
  double *iajy_ptr;
  double *ixjs_ptr;
  double *ikja_ptr;
  double *iajs_ptr, *ikjs_ptr;
  double *iqrs_ptr, *iprs_ptr;
  double *iqjs_ptr, *iqjr_ptr;
  const double *pqrs_ptr;

  double **gradient, *gradient_dat;  // The MP2 gradient
  double **ginter;    // Intermediates for the MP2 gradient

  BiggestContribs biggest_coefs(5,10);
  CharacterTable ct = molecule()->point_group()->char_table();

  int dograd = do_gradient();

  // this controls how often mem->catchup is called
  int catchup_ctr;
  const int catchup_mask = 3;

  tim_enter("mp2-mem");

  nfuncmax = basis()->max_nfunction_in_shell();

  DerivCenters der_centers;

  nshell = basis()->nshell();

  me = msg_->me();

  if (me == 0) {
    cout << endl << indent
         << "Entered memgrp based MP2 routine" << endl;
    }
  
  nproc = msg_->n();
  if (me == 0)
    cout << indent << scprintf("nproc = %i", nproc) << endl;

  tol = (int) (-10.0/log10(2.0));  // discard ereps smaller than 10^-10
  dtol = 1.0e-10;

  nocc = 0;
  for (i=0; i<nbasis; i++) {
    if (reference_->occupation(i) == 2.0) nocc++;
    }

  nocc_act = nocc - nfzc;
  nvir  = nbasis - nocc;
  nvir_act = nvir - nfzv;

  // Do a few preliminary tests to make sure the desired calculation
  // can be done (and appears to be meaningful!)

  if (nocc_act <= 0) {
    if (me == 0) {
      cerr << "There are no active occupied orbitals; program exiting" << endl;
      }
    abort();
    }

  if (nvir_act <= 0) {
    if (me == 0) {
      cerr << "There are no active virtual orbitals; program exiting" << endl;
      }
    abort();
    }

  ////////////////////////////////////////////////////////
  // Compute batch size ni for mp2 loops;
  //
  // The following arrays are kept throughout (all of type double):
  //   scf_vector, gradient, ginter, Pkj, Pab, Wkj, Wab, Waj, Laj, iqjs_buf
  // and memory allocated for these arrays is called mem_static
  //
  ////////////////////////////////////////////////////////
  if (me == 0) {
    mem_static = nbasis*nbasis; // scf vector
    mem_static += 2*nbasis*nfuncmax; // iqjs & iqjr
    if (dograd) {
      mem_static += 6*natom; // gradient & ginter
      mem_static += (nocc*(nocc+1))/2; // Pkj
      mem_static += (nvir*(nvir+1))/2; // Pab
      mem_static += nocc*nocc; // Wkj
      mem_static += nvir*nvir; // Wab
      mem_static += 2*nocc*nvir; // Waj & Laj
      }
    mem_static *= sizeof(double);
    ni = compute_cs_batchsize(mem_static, nocc_act); 
    }

  // Send value of ni and mem_static to other nodes
  msg_->bcast(ni);
  msg_->bcast(mem_static);

  // Compute the storage remaining for the integral routines
  distsize_t dyn_mem = compute_cs_dynamic_memory(ni,nocc_act);
  int mem_remaining;
  if (mem_alloc <= (dyn_mem + mem_static)) mem_remaining = 0;
  else mem_remaining = mem_alloc - (int)(dyn_mem + mem_static);

  cout << node0 << indent
       << "Memory available per node:      " << mem_alloc << " Bytes"
       << endl;
  cout << node0 << indent
       << "Static memory used per node:    " << mem_static << " Bytes"
       << endl;
  cout << node0 << indent
       << "Total memory used per node:     " << dyn_mem+mem_static << " Bytes"
       << endl;
  cout << node0 << indent
       << "Memory required for one pass:   "
       << compute_cs_dynamic_memory(nocc_act,nocc_act)+mem_static
       << " Bytes"
       << endl;
  cout << node0 << indent
       << "Minimum memory required:        "
       << compute_cs_dynamic_memory(1,nocc_act)+mem_static
       << " Bytes"
       << endl;
  cout << node0 << indent
       << "Batch size:                     " << ni
       << endl;

  if (ni == 0) {
    cerr << node0 << "Batch size is 0: more memory or processors are needed"
         << endl;
    abort();
    }

  if (ni == nocc_act) {
    npass = 1;
    rest = 0;
    }
  else {
    rest = nocc_act%ni;
    npass = (nocc_act - rest)/ni + 1;
    if (rest == 0) npass--;
    }

  if (me == 0) {
    cout << indent
         << scprintf(" npass  rest  nbasis  nshell  nfuncmax") << endl;
    cout << indent
         << scprintf("  %-4i   %-3i   %-5i    %-4i     %-3i",
                     npass,rest,nbasis,nshell,nfuncmax)
         << endl;
    cout << indent
         << scprintf(" nocc   nvir   nfzc   nfzv") << endl;
    cout << indent
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

  RefDiagSCMatrix occ;
  RefSCMatrix Scf_Vec;
  RefDiagSCMatrix evalmat;
  eigen(evalmat, Scf_Vec, occ);

  if (debug_ > 1) {
    evalmat.print("eigenvalues");
    Scf_Vec.print("eigenvectors");
    }

  double *scf_vector_dat = new double[nbasis*nbasis];
  Scf_Vec.t()->convert(scf_vector_dat);

  int *orbital_map = new int[nbasis];
  evals = new double[nbasis];
  double** scf_vector = new double*[nbasis];
  int idoc=0, ivir=0;
  const double epsilon = 1.0e-4;
  for (i=0; i<nbasis; i++) {
    if (occ(i) >= 2.0 - epsilon) {
      evals[idoc] = evalmat(i);
      scf_vector[idoc] = &scf_vector_dat[i*nbasis];
      orbital_map[idoc] = i;
      idoc++;
      }
    else {
      evals[ivir+nocc] = evalmat(i);
      scf_vector[ivir+nocc] = &scf_vector_dat[i*nbasis];
      orbital_map[ivir+nocc] = i;
      ivir++;
      }
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
    Laj            = (double*) malloc(nvir*nocc*sizeof(double));

    gradient_dat = new double[natom*3];
    gradient = new double*[natom];
    for (i=0; i<natom; i++) {
      gradient[i] = &gradient_dat[i*3];
      }

    ginter = new double*[natom];
    for (i=0; i<natom; i++) {
      ginter[i] = new double[3];
      for (xyz=0; xyz<3; xyz++) ginter[i][xyz] = 0;
      }

    //////////////////////////////
    // Initialize various arrays
    //////////////////////////////

    bzerofast(Pkj,nocc*(nocc+1)/2);
    bzerofast(Wkj,nocc*nocc);
    bzerofast(Pab,nvir*(nvir+1)/2);
    bzerofast(Wab,nvir*nvir);
    bzerofast(Waj,nvir*nocc);
    bzerofast(Laj,nvir*nocc);

    if (me == 0) zero_gradients(gradient, natom, 3);
    }
  else {
    Pkj = 0;
    Pab = 0;
    Wkj = 0;
    Wab = 0;
    Waj = 0;
    Laj = 0;

    gradient_dat = 0;
    gradient = 0;
    ginter = 0;
    }

  if (debug_ && me == 0) {
    for (j=0; j<nbasis; j++) {
      cout << indent
           << scprintf("eigenvalue[%3d] = %15.10lf", j, evals[j]);
      if (j < nfzc) cout << " (frozen docc)";
      else if (j < nocc_act + nfzc) cout << " (active docc)";
      else if (j < nvir_act + nocc_act + nfzc) cout << " (active uocc)";
      else cout << " (frozen uocc)";
      cout << endl;
      }
    }

  if (nproc > 1) iqjs_buf = new double[2 + nfuncmax*nbasis];

  /////////////////////////////////////
  //  Begin MP2 loops
  /////////////////////////////////////

  // debug print
  if (debug_ && me == 0) {
    cout << indent
         << scprintf("node %i, begin loop over i-batches",me) << endl;
    }
  // end of debug print

  // Initialize the integrals
  integral()->set_storage(mem_remaining);
  tbint_ = integral()->electron_repulsion();
  intbuf = tbint_->buffer();
  tbintder_ = integral()->electron_repulsion_deriv();
  intderbuf = tbintder_->buffer();

  int mem_integral_intermediates = integral()->storage_used();
  int mem_integral_storage = mem_remaining - mem_integral_intermediates;
  if (mem_integral_storage<0) mem_integral_storage = 0;
  tbint_->set_integral_storage(mem_integral_storage);

  cout << node0 << endl << indent
       << scprintf("Memory used for integral intermediates: %i Bytes",
                   mem_integral_intermediates)
       << endl;
  cout << node0 << indent
       << scprintf("Memory used for integral storage:       %i Bytes",
                   mem_integral_storage)
       << endl;

  if (mem.null()) {
      cerr << "MBPT2: memory group not initialized" << endl;
      abort();
    }

  mem->set_localsize(nijmax*nbasis*nbasis*sizeof(double));

  mem->lock(0);

  MemoryGrpBuf<double> membuf(mem);
  MemoryGrpBuf<double> membuf_remote(mem);

  tim_enter("mp2 passes");
  for (pass=0; pass<npass; pass++) {

    if (debug_) {
      cout << indent << me << " beginning pass " << pass << endl;
      }

    i_offset = pass*ni + nfzc;
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
      cout << indent << "node " << me << ", nij = " << nij << endl;
    // end of debug print

    mem->sync(); // This must be here or gamma non-sep will be wrong when running
                 // on multiple processors with more than one pass

    r_offset = 0;

    // Allocate and initialize some arrays
    // (done here to avoid having these arrays
    // overlap with arrays allocated later)

    // Allocate (and initialize) some arrays

    integral_iqrs = new double[ni*nbasis*nfuncmax*nfuncmax];
    iqjs_contrib  = new double[nbasis*nfuncmax];
    iqjr_contrib  = new double[nbasis*nfuncmax];

    integral_iqjs = membuf.writeonly_on_node(0, nij*nbasis*nbasis);

    bzerofast(integral_iqjs, nij*nbasis*nbasis);

    integral_iqjs = 0;
    membuf.release();
    mem->lock(1);
    mem->sync();

    index = 0;

    // debug print
    if (debug_ && me == 0) {
      cout << indent
           << scprintf("Begin loop over shells (erep, 1.+2. qt)") << endl;
      }
    // end of debug print

    for (S=0; S<nshell; S++) {
      ns = basis()->shell(S).nfunction();
      s_offset = basis()->shell_to_function(S);

      for (R=0; R<=S; R++) {
        nr = basis()->shell(R).nfunction();
        r_offset = basis()->shell_to_function(R);

        if (index++ % nproc == me) {

          bzerofast(integral_iqrs, ni*nbasis*nfuncmax*nfuncmax);

          for (Q=0; Q<nshell; Q++) {
            nq = basis()->shell(Q).nfunction();
            q_offset = basis()->shell_to_function(Q);
            for (P=0; P<=Q; P++) {
              np = basis()->shell(P).nfunction();
              p_offset = basis()->shell_to_function(P);

           // if (scf_erep_bound(P,Q,R,S) < tol) {
           //   continue;  // skip ereps less than tol
           //   }
              if (tbint_->log2_shell_bound(P,Q,R,S) < tol) {
                continue;  // skip ereps less than tol
                }

              aoint_computed++;

              tim_enter("erep");
              tbint_->compute_shell(P,Q,R,S);
              tim_exit("erep");

              mem->catchup();

              tim_enter("1. q.t.");
              // Begin first quarter transformation;
              // generate (iq|rs) for i active

              offset = nr*ns*nbasis;
              pqrs_ptr = intbuf;
              for (bf1 = 0; bf1 < np; bf1++) {
                p = p_offset + bf1;
                for (bf2 = 0; bf2 < nq; bf2++) {
                  q = q_offset + bf2;

                  if (q < p) {
                    pqrs_ptr = &intbuf[ns*nr*(bf2+1 + nq*bf1)];
                    continue; // skip to next q value
                    }

                  for (bf3 = 0; bf3 < nr; bf3++) {
                    r = r_offset + bf3;

                    for (bf4 = 0; bf4 < ns; bf4++) {
                      s = s_offset + bf4;

                      if (s < r) {
                        pqrs_ptr++;
                        continue; // skip to next bf4 value
                        }

                      if (fabs(*pqrs_ptr) > dtol) {
                        iprs_ptr = &integral_iqrs[bf4 + ns*(p + nbasis*bf3)];
                        iqrs_ptr = &integral_iqrs[bf4 + ns*(q + nbasis*bf3)];
                        c_qi = &scf_vector[q][i_offset];
                        c_pi = &scf_vector[p][i_offset];
                        tmpval = *pqrs_ptr;
                        for (i=0; i<ni; i++) {
                          *iprs_ptr += *c_qi++*tmpval;
                          iprs_ptr += offset;
                          if (p != q) {
                            *iqrs_ptr += *c_pi++*tmpval;
                            iqrs_ptr += offset;
                            }
                          } // exit i loop
                        }   // endif

                      pqrs_ptr++;
                      } // exit bf4 loop
                    }   // exit bf3 loop
                  }     // exit bf2 loop
                }       // exit bf1 loop
              // end of first quarter transformation
              tim_exit("1. q.t.");

              }           // exit P loop
            }             // exit Q loop

#if PRINT1Q
      {
        double *tmp = integral_iqrs;
        for (int i = 0; i<ni; i++) {
            for (int r = 0; r<nr; r++) {
                for (int q = 0; q<nbasis; q++) {
                    for (int s = 0; s<ns; s++) {
                        printf("1Q: (%d %d|%d %d) = %12.8f\n",
                               i,q,r+r_offset,s+s_offset,*tmp);
                        tmp++;
                      }
                  }
              }
          }
      }
#endif

          tim_enter("2. q.t.");
          // Begin second quarter transformation;
          // generate (iq|jr) for i active and j active or frozen
          for (i=0; i<ni; i++) {
            for (j=0; j<nocc; j++) {

              bzerofast(iqjs_contrib, nbasis*nfuncmax);
              bzerofast(iqjr_contrib, nbasis*nfuncmax);

              catchup_ctr = 0;
              for (bf1=0; bf1<ns; bf1++) {
                s = s_offset + bf1;
                c_sj = &scf_vector[s][j];
                iqjr_ptr = iqjr_contrib;
                for (bf2=0; bf2<nr; bf2++) {
                  r = r_offset + bf2;
                  if (r > s) {
                    break; // skip to next bf1 value
                    }
                  c_rj = scf_vector[r][j];
                  iqjs_ptr = &iqjs_contrib[bf1*nbasis];
                  iqrs_ptr = &integral_iqrs[bf1 + ns*nbasis*(bf2 + nr*i)];
                  for (q=0; q<nbasis; q++) {
                    *iqjs_ptr++ += c_rj * *iqrs_ptr;
                    if (r != s) *iqjr_ptr += *c_sj * *iqrs_ptr;
                    iqjr_ptr++;
                    iqrs_ptr += ns;
                    } // exit q loop
                  // every so often process outstanding messages
                  if ((catchup_ctr++ & catchup_mask) == 0) mem->catchup();
                  }   // exit bf2 loop
                }     // exit bf1 loop

              // We now have contributions to iqjs and iqjr for one pair i,j,
              // all q, r in R and s in S; send iqjs and iqjr to the node
              // (ij_proc) which is going to have this ij pair
              ij_proc =  (i*nocc + j)%nproc;
              ij_index = (i*nocc + j)/nproc;

              // Sum the iqjs_contrib to the appropriate place
              ij_offset = nbasis*(s_offset + nbasis*ij_index);
              mem->sum_reduction_on_node(iqjs_contrib,
                                         ij_offset, ns*nbasis, ij_proc);

              ij_offset = nbasis*(r_offset + nbasis*ij_index);
              mem->sum_reduction_on_node(iqjr_contrib,
                                         ij_offset, nr*nbasis, ij_proc);

              }     // exit j loop
            }       // exit i loop
          // end of second quarter transformation
          tim_exit("2. q.t.");

          }     // endif
        }       // exit R loop
      }         // exit S loop
    // debug print
    if (debug_ && me == 0) {
      cout << indent << "End of loop over shells" << endl;
      }
    // end of debug print

    mem->lock(0);
    mem->sync();  // Make sure iqjs is complete on each node before continuing

    integral_iqjs = membuf.readwrite_on_node(0, nij*nbasis*nbasis);

//     if (me == 0) {
//         const int maxint = 300;
//         int index = 0;
//         int nint = 0;
//         const double *tmpint;
//         for (i=0; i<ni; i++) {
//             for (j=0; j<nocc; j++) {
//                 ij_proc =  (i*nocc + j)%nproc; // ij_proc has this ij pair
//                 ij_index = (i*nocc + j)/nproc;
//                 ij_offset = ij_index*nbasis*nbasis;
//                 tmpint = membuf.readonly_on_node(ij_offset, 1, ij_proc);
//                 for (s=0; s<nbasis; s++) {
//                     for (q=0; q<nbasis; q++) {
//                         printf("(%d %d|%d %d) = %12.5f (node %d)\n",
//                                i,q,j,s, *tmpint, ij_proc);
//                         tmpint++;
//                         nint++;
//                         if (nint > maxint) break;
//                       }
//                     if (nint > maxint) break;
//                   }
//                 membuf.release();
//                 if (nint > maxint) break;
//               }
//             if (nint > maxint) break;
//           }
//       }
//     fflush(stdout);
//     mem->sync();

    delete[] integral_iqrs;
    delete[] iqjs_contrib;
    delete[] iqjr_contrib;

    // Allocate and initialize some arrays
    ixjs_tmp = new double[nbasis];

    // debug print
    if (debug_ && me == 0) {
      cout << indent << "Begin 3. qt" << endl;
      }
    // end of debug print

    tim_enter("3. q.t.");
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
              for (x=0; x<nbasis; x++) {
                *ixjs_ptr++ += *c_qx++ * tmpval;
                }
              }   // exit q loop

            // Put ixjs into integral_iqjs, while overwriting what was there;
            // i.e., integral_iqjs will now contain three-quarter transformed
            // integrals ixjs
            integral_iqjs_ptr = &integral_iqjs[nbasis*(s + nbasis*ij_index)];
            ixjs_ptr = ixjs_tmp;
            for (x=0; x<nbasis; x++) {
              *integral_iqjs_ptr++ = *ixjs_ptr++;
              }
            }   // exit s loop
          ij_index++;
          }     // endif
        }       // exit j loop
      }         // exit i loop
    // end of third quarter transformation
    tim_exit("3. q.t.");

    // debug print
    if (debug_ && me == 0) {
      cout << indent << "End of 3. qt" << endl;
      }
    // end of debug print

    delete[] ixjs_tmp;

    // The array of half-transformed integrals integral_iqjs has now
    // been overwritten by three-quarter transformed integrals ixjs;
    // rename the array integral_ixjs, where x = any MO
    integral_ixjs = integral_iqjs;

    integral_iajy = new double[nbasis];
    // in iajy: i act; a,j act or frz; y act or frz occ or virt.
    integral_ikja = new double[nvir_act];
    // in ikja: i,j act; k act or frz; a act.

    // debug print
    if (debug_ && me == 0) {
      cout << indent << "Begin 4. qt" << endl;
      }
    // end of debug print

    // Begin fourth quarter transformation
    // generating MO integrals (ov|ov), (ov|oo) and (oo|ov)
    tim_enter("4. q.t.");
    index = 0;
    ij_index = 0;
    for (i=0; i<ni; i++) {
      for (j=0; j<nocc; j++) {
        if (index++ % nproc == me) {

          for (a=0; a<nvir; a++) {
            bzerofast(integral_iajy, nbasis);
            iajs_ptr = &integral_ixjs[a+nocc + nbasis*nbasis*ij_index];
            for (s=0; s<nbasis; s++) {
              c_sy = scf_vector[s];
              iajy_ptr = integral_iajy;
              tmpval = *iajs_ptr;
              for (y=0; y<nbasis; y++) {
                *iajy_ptr++ += *c_sy++ * tmpval;
                } // exit y loop
              iajs_ptr += nbasis;
              }   // exit s loop
            // Put integral_iajy into ixjs for one i,a,j while
            // overwriting elements of ixjs
            iajs_ptr = &integral_ixjs[a+nocc + nbasis*nbasis*ij_index];
            iajy_ptr = integral_iajy;
            for (y=0; y<nbasis; y++) {
              *iajs_ptr = *iajy_ptr++;
              iajs_ptr += nbasis;
              } // exit y loop
            }   // exit a loop

          // this is only needed for gradients
          if (dograd && j >= nfzc) {
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
    tim_exit("4. q.t.");

    // debug print
    if (debug_ && me == 0) {
      cout << indent << "End of 4. qt" << endl;
      }
    // end of debug print

    // The array integral_ixjs has now been overwritten by MO integrals
    // iajy and ikja, so rename the array mo_int
    mo_int = integral_ixjs;

    delete[] integral_iajy;
    delete[] integral_ikja;

    // Divide the (ia|jb) MO integrals by the term 
    // evals[i]+evals[j]-evals[a]-evals[b]
    // and keep these integrals in mo_int;
    // do this only for active i, j, a, and b
    tim_enter("divide (ia|jb)'s");

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
               *iajb_ptr++ /= evals[ii]+evals[j]-evals[a+nocc]-evals[bb];
                } // exit a loop
              }   // exit b loop
            }     // endif
          ij_index++;
          }       // endif
        }         // exit j loop
      }           // exit i loop
    tim_exit("divide (ia|jb)'s");

    // We now have the fully transformed integrals (ia|jb)
    // divided by the proper orbital energy denominators
    // for one batch of i, all j<nocc, and all a<nvir and b<nvir,
    // where i, j, a, and b are all active;
    // compute contribution to the MP2 correlation energy
    // from these integrals 

    tim_enter("compute ecorr");

    index = 0;
    ij_index = 0;
    for (i=0; i<ni; i++) {
      for (j=0; j<nocc; j++) {
        if (index++ % nproc == me) {

          if (j>=nfzc) {
            for (b=0; b<nvir_act; b++) {
              iajb_ptr = &mo_int[nocc + nbasis*(b+nocc + nbasis*ij_index)];
              ibja_ptr = &mo_int[b+nocc + nbasis*(nocc + nbasis*ij_index)];
              for (a=0; a<nvir_act; a++) {
                delta_ijab = evals[i_offset+i]+evals[j]-evals[nocc+a]-evals[nocc+b];
                // only include determinants with unique coefficients
                if (a>=b && i_offset+i>=j) {
                  // aaaa or bbbb
                  if (a>b && i_offset+i>j) {
                    biggest_coefs.insert(*iajb_ptr - *ibja_ptr,
                                         i_offset+i,j,a,b,1111);
                    }
                  // abab or baba or abba or baab
                  biggest_coefs.insert(*iajb_ptr,i_offset+i,j,a,b,1212);
                  }

                ecorr_mp2 += *iajb_ptr*(2**iajb_ptr - *ibja_ptr)*delta_ijab;
                iajb_ptr++;
                ibja_ptr += nbasis;;
                } // exit a loop
              }   // exit b loop
            }     // endif

          ij_index++;
          }       // endif
        }         // exit j loop
      }           // exit i loop
    tim_exit("compute ecorr");

    // debug print
    if (debug_ && me == 0) {
      cout << indent << "End of ecorr" << endl;
      }
    // end of debug print

    if (npass > 1 && pass < npass - 1) {
      double passe = ecorr_mp2;
      msg_->sum(passe);
      cout << node0 << indent
           << "Partial correlation energy for pass " << pass << ":" << endl;
      cout << node0 << indent
           << scprintf("  restart_ecorr_memgrp   = %14.10f", passe)
           << endl;
      cout << node0 << indent
           << scprintf("  restart_orbital_memgrp = %d", ((pass+1) * ni))
           << endl;
      }

    // don't go beyond this point if only the energy is needed
    if (!dograd) continue;

    integral_iqjs = 0;
    membuf.release();
    mem->sync(); // Make sure MO integrals are complete on all nodes before continuing

    mo_int = (double*) membuf.readonly_on_node(0, nij*nbasis*nbasis);

    // Update the matrices Pkj and Wkj with
    // contributions from (occ vir|occ vir) integrals
    index = 0;
    ij_index = 0;
    tim_enter("Pkj and Wkj");
    for (i=0; i<ni; i++) {
      for (j=0; j<nocc; j++) {
        if (index++ % nproc == me) {

          if (j>=nfzc) {
            for (kloop=me; kloop<me+nocc; kloop++) {
              // stagger k's to minimize contention
              k = kloop%nocc;
              if (k<=j) pkj_ptr = &Pkj[j*(j+1)/2 + k];
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
                  tmpval = 2**iakb_ptr * (*ibja_ptr++ - 2 * *iajb_ptr);
                  iakb_ptr += nbasis;
                  iajb_ptr += nbasis;
                  if (k<nfzc && k<=j) *pkj_ptr += tmpval/(evals[k]-evals[j]);
                  else if (k<=j) *pkj_ptr += tmpval;
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
    tim_exit("Pkj and Wkj");

    // debug print
    if (debug_ && me == 0) {
      cout << indent << "End of Pkj and Wkj" << endl;
      }
    // end of debug print

    // Update the matrices Pab and Wab with
    // contributions from (occ vir|occ vir) integrals
    tim_enter("Pab and Wab");
    index = 0;
    ij_index = 0;
    for (i=0; i<ni; i++) {
      for (j=0; j<nocc; j++) {
        if (index++ % nproc == me) {
          if (j>=nfzc) {

            offset = nocc + nocc*nbasis + nbasis*nbasis*ij_index;
            for (a=0; a<nvir_act; a++) {
              pab_ptr = &Pab[a*(a+1)/2];
              for (b=0; b<=a; b++) {  // active-active part of Pab and Wab
                wab_ptr = &Wab[a*nvir + b];
                wba_ptr = &Wab[b*nvir + a];
                ibjc_ptr = &mo_int[offset + b];
                icjb_ptr = &mo_int[offset + b*nbasis];
                iajc_ptr = &mo_int[offset + a];
                icja_ptr = &mo_int[offset + a*nbasis];
                for (c=0; c<nvir_act; c++) {
                  *pab_ptr += 2**iajc_ptr * (2 * *ibjc_ptr - *icjb_ptr);
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
    tim_exit("Pab and Wab");

    // debug print
    if (debug_ && me == 0) {
      cout << indent << "End of Pab and Wab" << endl;
      }
    // end of debug print

    ///////////////////////////////////////
    // Update Waj and Laj with contrib. 
    // from (oo|ov) and (ov|oo) integrals;
    // here a is active and j is active or
    // frozen
    ///////////////////////////////////////
    tim_enter("Waj and Laj");

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
                waj_ptr = &Waj[j*nvir]; // order as j*nvir+a to make loops more efficient
                laj_ptr = &Laj[j*nvir];
                for (a=0; a<nvir_act; a++) {
                  tmpval = 2**ibka_ptr * *ijkb_ptr;
                  ibka_ptr += nbasis;
                  *waj_ptr++ += tmpval;
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
                waj_ptr = &Waj[j*nvir];
                laj_ptr = &Laj[j*nvir];
                for (a=0; a<nvir_act; a++) {
                  tmpval = 4 * *ibka_ptr * *ibkj_ptr;
                  ibka_ptr += nbasis;
                  *waj_ptr++ -= tmpval;
                  *laj_ptr++ += tmpval; // This term had the wrong sign in Frisch's paper
                  } // exit a loop
                }   // exit j loop
              }     // exit b loop
            }       // endif
          ik_index++;
          }       // endif
        }         // exit k loop
      }           // exit i loop

    tim_exit("Waj and Laj");
    /////////////////////////////
    // End of Waj and Laj update
    /////////////////////////////

    // debug print
    if (debug_ && me == 0) {
      cout << indent << "End of Paj and Waj" << endl;
      }
    // end of debug print

    mo_int = 0;
    membuf.release();

    mem->sync(); // Need to synchronize before deleting mo_intbuf

    mo_int = membuf.readwrite_on_node(0, nij*nbasis*nbasis);

    gamma_iajs_tmp = new double[nbasis*nvir_act];
    if (!gamma_iajs_tmp) {
      cout << indent << "Could not allocate gamma_iajs_tmp" << endl;
      abort();
      }

    // debug print
    if (debug_ && me == 0) {
      cout << indent << "Begin 1+2qbt\n" << endl;
      }
    // end of debug print

    ///////////////////////////////////////////////////////////
    // Perform first and second quarter back-transformation.
    // Each node produces gamma_iajs, and gamma_iqjs 
    // for a subset of i and j, all a and all s;
    // the back-transf. is done only for active i, j, a, and b
    ///////////////////////////////////////////////////////////

    // Begin first quarter back-transformation
    tim_enter("1. q.b.t.");
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
    tim_exit("1. q.b.t.");

    // debug print
    if (debug_ && me == 0) {
      cout << indent << "End 1. qbt" << endl;
      }
    // end of debug print

    mo_int = 0;
    membuf.release();

    mem->sync(); // Make sure all nodes are done with gamma_iajs_tmp before renaming

    delete[] gamma_iajs_tmp;

    // The array mo_int has now been overwritten by the quarter 
    // back-transformed non-sep 2PDM gamma_iajs, so rename
    gamma_iajs = membuf.readwrite_on_node(0, nij*nbasis*nbasis);

    gamma_iqjs_tmp = new double[nbasis];
    if (!gamma_iqjs_tmp) {
      cerr << "Could not allocate gamma_iqjs_tmp" << endl;
      abort();
      }

    if (debug_ && me == 0) {
      cout << indent << "Begin 2. qbt" << endl;
      }

    // Begin second quarter back-transformation
    // (gamma_iqjs elements ordered as i,j,s,q,
    // i.e., q varies fastest)
    tim_enter("2. q.b.t.");
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
    tim_exit("2. q.b.t.");
    // end of second quarter back-transformation

    if (debug_ && me == 0) {
      cout << indent << "End 2. qbt" << endl;
      }

    gamma_iajs = 0;
    membuf.release();
    
    mem->sync(); // Keep this here to make sure all nodes have gamma_iqjs
                 // before it is needed below, and that gamma_iajs is not
                 // deleted prematurely

    // The quarter back-transformed elements gamma_iajs have now been
    // overwritten by the half back-transformed elements gamma_iqjs

    delete[] gamma_iqjs_tmp;

    /////////////////////////////////////////////////
    // End of 1. and 2. quarter back-transformation
    /////////////////////////////////////////////////

    // Allocate various arrays
    gamma_iqrs = new double[ni*nbasis*nfuncmax*nfuncmax];
    if (!gamma_iqrs) {
      cerr << "Could not allocate gamma_iqrs" << endl;
      abort();
      }
    gamma_pqrs = new double[nfuncmax*nfuncmax*nfuncmax*nfuncmax];
    if (!gamma_pqrs) {
      cerr << "Could not allocate gamma_pqrs" << endl;
      abort();
      }

    Lpi = new double[nbasis*ni];
    bzerofast(Lpi,nbasis*ni);

    if (debug_ && me == 0) {
      cout << indent << "Begin 3. and 4. qbt" << endl;
      }

    ////////////////////////////////////////////////////////
    // Perform third and fourth quarter back-transformation
    // and compute contrib. to gradient from non-sep 2PDM
    ////////////////////////////////////////////////////////

    index = 0;

    for (S=0; S<nshell; S++) {
      ns = basis()->shell(S).nfunction();
      s_offset = basis()->shell_to_function(S);

      for (R=0; R<=S; R++) {
        nr = basis()->shell(R).nfunction();
        r_offset = basis()->shell_to_function(R);

        // If both PQRS and PQRS derivative are zero, skip to next R
        if (tbint_->log2_shell_bound(R,S) < tol
            && tbintder_->log2_shell_bound(R,S) < tol) continue;

        if (index++ % nproc == me) {

          tim_enter("3. q.b.t.");
          // Begin third quarter back-transformation.

          bzerofast(gamma_iqrs,ni*nbasis*nfuncmax*nfuncmax);

          for (i=0; i<ni; i++) {
            for (jloop=me; jloop<me+nocc_act; jloop++) {
              // stagger j's to minimize contention
              j = jloop%nocc_act + nfzc;  // j runs from nfzc to nocc
              ij_proc =  (i*nocc + j)%nproc; // ij_proc has this ij pair
              ij_index = (i*nocc + j)/nproc;

              offset = s_offset*nbasis + ij_index*nbasis*nbasis;
              // Send for elements gamma_iqjs, if necessary
              gammabuf = (double*) membuf_remote.readonly_on_node(offset,
                                                                  nbasis * ns,
                                                                  ij_proc);
              for (bf1=0; bf1<nr; bf1++) {
                c_rj = scf_vector[bf1 + r_offset][j];
                gamma_iqjs_ptr = gammabuf;
                for (bf2=0; bf2<ns; bf2++) {
                  gamma_iqrs_ptr = &gamma_iqrs[bf2 + ns*nbasis*(bf1 + nr*i)];
                  for (q=0; q<nbasis; q++) {
                    *gamma_iqrs_ptr += c_rj * *gamma_iqjs_ptr++;
                    gamma_iqrs_ptr += ns;
                    } // exit q loop
                  }   // exit bf2 loop
                }     // exit bf1 loop

              membuf_remote.release();

              offset = r_offset*nbasis + ij_index*nbasis*nbasis;
              // Send for elements gamma_irjq, if necessary
              gammabuf = (double*) membuf_remote.readonly_on_node(offset,
                                                                  nbasis*nr,
                                                                  ij_proc);
              for (bf1=0; bf1<ns; bf1++) {
                s = bf1 + s_offset;
                c_sj = &scf_vector[s][j];
                gamma_irjq_ptr = gammabuf;
                for (bf2=0; bf2<nr; bf2++) {
                  r = bf2 + r_offset;
                  if (r != s) {
                    gamma_iqsr_ptr = &gamma_iqrs[bf1 + ns*nbasis*(bf2 + nr*i)];
                    for (q=0; q<nbasis; q++) {
                      *gamma_iqsr_ptr += *c_sj * *gamma_irjq_ptr++;
                      gamma_iqsr_ptr += ns;
                      } // exit q loop
                    }   // endif
                  else gamma_irjq_ptr += nbasis;
                  }     // exit bf2 loop
                }       // exit bf1 loop

              membuf_remote.release();

              }         // exit j loop
            }           // exit i loop

          // end of third quarter back-transformation
          // we now have gamma_iqrs (symmetrized)
          // for i-batch, all q, s in S, r in R
          tim_exit("3. q.b.t.");

          // only do this if integral is nonzero
          if (tbint_->log2_shell_bound(R,S) >= tol) {

            // Compute contrib to Laj from (ov|vv) integrals
            // (done in AO basis to avoid generating (ov|vv);
            // here, generate Lpi for i-batch; later, transform
            // Lpi to get contribution to Laj
            tim_enter("(ov|vv) contrib to Laj");
            for (Q=0; Q<nshell; Q++) {
              nq = basis()->shell(Q).nfunction();
              q_offset = basis()->shell_to_function(Q);
              for (P=0; P<=Q; P++) {
                np = basis()->shell(P).nfunction();
                p_offset = basis()->shell_to_function(P);
             // if (scf_erep_bound(P,Q,R,S) < tol) {
             //   continue;  // skip ereps less than tol
             //   }
                if (tbint_->log2_shell_bound(P,Q,R,S) < tol) {
                  continue;  // skip ereps less than tol
                  }
                tim_enter("erep");
                tbint_->compute_shell(P,Q,R,S);
                tim_exit("erep");

                mem->catchup();

                offset = nr*ns*nbasis;
                int_index = 0;

                for (bf1 = 0; bf1 < np; bf1++) {
                  p = p_offset + bf1;
                  for (bf2 = 0; bf2 < nq; bf2++) {
                    q = q_offset + bf2;

                    if (q < p) {
                      int_index = ns*nr*(bf2+1 + nq*bf1);
                      continue; // skip to next q value
                      }

                    for (bf3 = 0; bf3 < nr; bf3++) {
                      r = r_offset + bf3;

                      for (bf4 = 0; bf4 < ns; bf4++) {

                        if (fabs(intbuf[int_index]) > dtol) {
                          s = s_offset + bf4;

                          if (s < r) {
                            int_index++;
                            continue; // skip to next bf4 value
                            }

                          gamma_iqrs_ptr = &gamma_iqrs[bf4 + ns*(q + nbasis*bf3)];
                          gamma_iprs_ptr = &gamma_iqrs[bf4 + ns*(p + nbasis*bf3)];
                          pqrs = intbuf[int_index];

                          lpi_ptr = &Lpi[p*ni];
                          lqi_ptr = &Lpi[q*ni];

                          for (i=0; i<ni; i++) {
                            *lpi_ptr++ -= pqrs**gamma_iqrs_ptr;
                            if (p != q) {
                              *lqi_ptr++ -= pqrs**gamma_iprs_ptr;
                              }
                            gamma_iqrs_ptr += offset;
                            gamma_iprs_ptr += offset;
                            } // exit i loop
                          // every so often process outstanding messages
                          if ((catchup_ctr++
                              & catchup_mask) == 0) mem->catchup();
                          }   // endif

                        int_index++;
                        }     // exit bf4 loop
                      }       // exit bf3 loop
                    }         // exit bf2 loop
                  }           // exit bf1 loop

                }             // exit P loop
              }               // exit Q loop
            tim_exit("(ov|vv) contrib to Laj");
            }                 // endif

          if (tbintder_->log2_shell_bound(R,S) >= tol) {

            for (Q=0; Q<=S; Q++) {
              nq = basis()->shell(Q).nfunction();
              q_offset = basis()->shell_to_function(Q);

              for (P=0; P<=(Q==S ? R:Q); P++) {
                np = basis()->shell(P).nfunction();
                p_offset = basis()->shell_to_function(P);

                // If integral derivative is less than threshold skip to next P
                if (tbintder_->log2_shell_bound(P,Q,R,S) < tol) continue;
                aointder_computed++;

                tim_enter("4. q.b.t.");
                bzerofast(gamma_pqrs,nfuncmax*nfuncmax*nfuncmax*nfuncmax);

                offset = nr*ns*nbasis;

                // Begin fourth quarter back-transformation
                gamma_pqrs_ptr = gamma_pqrs;
                for (bf1=0; bf1<np; bf1++) {
                  p = bf1 + p_offset;
                  for (bf2=0; bf2<nr; bf2++) {
                    for (bf3=0; bf3<nq; bf3++) {
                      q = bf3 + q_offset;
                      for (bf4=0; bf4<ns; bf4++) {
                        c_pi = &scf_vector[p][i_offset];
                        c_qi = &scf_vector[q][i_offset];
                        gamma_iqrs_ptr = &gamma_iqrs[bf4 + ns*(q + nbasis*bf2)];
                        gamma_iprs_ptr = &gamma_iqrs[bf4 + ns*(p + nbasis*bf2)];
                        tmpval = 0.0;
                        for (i=0; i<ni; i++) {
                          tmpval += *c_pi * *gamma_iqrs_ptr;
                          if (p!=q) tmpval += *c_qi * *gamma_iprs_ptr;
                          c_pi++;
                          c_qi++;
                          gamma_iqrs_ptr += offset;
                          gamma_iprs_ptr += offset;
                          } // exit i loop
                        *gamma_pqrs_ptr += tmpval;
                        gamma_pqrs_ptr++;
                        }   // exit bf4 loop
                      }     // exit bf3 loop
                    }       // exit bf2 loop
                  }         // exit bf1 loop
                // end of fourth quarter back-transformation
                tim_exit("4. q.b.t.");
                // (we now have the contribution from one i-batch to the
                // non-separable part of the 2PDM for one shell block PQRS)

                // Evaluate derivative integrals
                tim_enter("erep derivs");
                tbintder_->compute_shell(P,Q,R,S,der_centers);
                tim_exit("erep derivs");

                mem->catchup();

                // Compute contribution to gradient from non-sep 2PDM
                // (i.e., contract derivative integrals with gamma_pqrs)
                int_index = 0;
                tim_enter("non-sep 2PDM contrib.");
                for (derset=0; derset<der_centers.n(); derset++) {
                  for (xyz=0; xyz<3; xyz++) {
                    grad_ptr1 = &ginter[der_centers.atom(derset)][xyz];
                    grad_ptr2 = &ginter[der_centers.omitted_atom()][xyz];
                    for (bf1=0; bf1<np; bf1++) {
                      p = bf1 + p_offset;
                      for (bf2=0; bf2<nq; bf2++) {
                        q = bf2 + q_offset;
                        qp = q*(q+1)/2 + p;
                        for (bf3=0; bf3<nr; bf3++) {
                          r = bf3 + r_offset;
                          gamma_pqrs_ptr = &gamma_pqrs[ns*(bf2 + nq*(bf3 + nr*bf1))];
                          for (bf4=0; bf4<ns; bf4++) {
                            s = bf4 + s_offset;
                            sr = s*(s+1)/2 + r;
                            if (q == s && p == r) factor_pqrs = 1;
                            else factor_pqrs = 2;
                            tmpval = intderbuf[int_index]*factor_pqrs**gamma_pqrs_ptr;
                            gamma_pqrs_ptr++;
                            if (q>=p && s>=r && (P != R || Q != S || sr >= qp)) {
                              *grad_ptr1 += tmpval;
                              if (der_centers.has_omitted_center())
                                *grad_ptr2 -= tmpval;
                               }
                            int_index++;
                            } // exit bf4 loop
                          }   // exit bf3 loop
                        }     // exit bf2 loop
                      }       // exit bf1 loop
                    }         // exit xyz loop
                  }           // exit derset loop
                tim_exit("non-sep 2PDM contrib.");

                } // exit P loop
              }   // exit Q loop
            }     // endif
          }       // endif
        }         // exit R loop
      }           // exit S loop

    if (debug_ && me == 0) {
      cout << indent << "End 3. and 4. qbt" << endl;
      }

    mem->sync(); // Make sure all nodes are done before deleting arrays

    if (debug_ > 1) {
      RefSCDimension ni_dim(new SCDimension(ni));
      RefSCDimension nbasis_dim(new SCDimension(nbasis));
      RefSCMatrix Lpi_mat(nbasis_dim, ni_dim, kit);
      Lpi_mat->assign(Lpi);
      Lpi_mat.print("Lpi");
      }

    if (debug_ && me == 0) {
      cout << indent << "Back-transform Lpi" << endl;
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

    delete[] gamma_iqrs;
    delete[] gamma_pqrs;

    if (debug_ && me == 0) {
      cout << indent << "Done with pass " << pass << endl;
      }
    }           // exit loop over i-batches (pass)
  tim_exit("mp2 passes");

  mem->set_localsize(0);

  if (nproc > 1) delete[] iqjs_buf;

  // debug print
  if (debug_ && me == 0) {
    cout << indent << "Exited loop over i-batches" << endl;
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

  if (me == 0) {
    emp2 = escf + ecorr_mp2;

    if (biggest_coefs.ncontrib()) {
      cout << endl << indent
           << "Largest first order coefficients (unique determinants):"
           << endl;
      }
    for (i=0; i<biggest_coefs.ncontrib(); i++) {
      int i0 = orbital_map[biggest_coefs.indices(i)[0]];
      int i1 = orbital_map[biggest_coefs.indices(i)[1]];
      int i2 = orbital_map[biggest_coefs.indices(i)[2] + nocc];
      int i3 = orbital_map[biggest_coefs.indices(i)[3] + nocc];
      int spincase = biggest_coefs.indices(i)[4];
      cout << indent
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

    // Print out various energies etc.

    if (debug_) {
      cout << indent << "Number of shell quartets for which AO integrals\n"
           << indent << "(or integral derivatives) would have been computed\n"
           << indent << "without bounds checking: "
           << npass*nshell*nshell*(nshell+1)*(nshell+1)/2
           << endl;

      cout << indent << "Number of shell quartets for which AO integrals\n"
           << indent << "were computed: " << aoint_computed
           << endl;

      if (dograd) {
        cout << indent
             << "Number of shell quartets for which AO integral derivatives\n"
             << indent << "were computed: " << aointder_computed
             << endl;
        }
      }

    cout<<endl<<indent
        <<scprintf("RHF energy [au]:                   %17.12lf\n", escf);
    cout<<indent
        <<scprintf("MP2 correlation energy [au]:       %17.12lf\n", ecorr_mp2);
    cout<<indent
        <<scprintf("MP2 energy [au]:                   %17.12lf\n", emp2);
    cout.flush();
    }
  if (method_ && strcmp(method_,"mp")) {
    cout << node0 << indent
         << "MBPT2: bad method for closed shell case: " << method_
         << ", using mp" << endl;
    }
  set_energy(emp2);

  // quit here if only the energy is needed
  if (!dograd) {
    delete[] scf_vector;
    delete[] scf_vector_dat;
    delete[] evals;
    tim_exit("mp2-mem");
    return;
    }

  // Accumulate intermediate gradients on node 0
  sum_gradients(msg_, ginter, natom, 3);

  // Add intermediate gradients to the gradient on node 0
  if (me==0)
    accum_gradients(gradient, ginter, natom, 3);

  // Print out contribution to the gradient from non-sep. 2PDM
  if (debug_ && me == 0) {
    cout << indent
         << "Contribution to MP2 gradient from non-separable 2PDM [au]:"
         << endl;
    for (i=0; i<natom; i++) {
      cout << indent << scprintf("%15.10lf  %15.10lf  %15.10lf",
                       ginter[i][0], ginter[i][1], ginter[i][2])
           << endl;
      }
    }

  ////////////////////////////////////////////////////////
  // Add contributions from all nodes to various matrices
  ////////////////////////////////////////////////////////
  tmpint = (nvir > nocc ? nvir:nocc);
  double *tmpmat = new double[tmpint*tmpint];
  msg_->sum(Laj,nvir*nocc,tmpmat);
  msg_->sum(Pkj,nocc*(nocc+1)/2,tmpmat); // Pkj is now complete
  msg_->sum(Pab,nvir*(nvir+1)/2,tmpmat); // Pab is now complete
  msg_->sum(Wab,nvir*nvir,tmpmat);
  msg_->sum(Wkj,nocc*nocc,tmpmat);
  msg_->sum(Waj,nvir*nocc,tmpmat);
  delete[] tmpmat;

  RefSCDimension nocc_act_dim(new SCDimension(nocc_act));
  RefSCDimension nvir_act_dim(new SCDimension(nvir_act));
  RefSCDimension nocc_dim(new SCDimension(nocc));
  RefSCDimension nvir_dim(new SCDimension(nvir));
  RefSCDimension nbasis_dim(new SCDimension(nbasis));

  RefSCMatrix s2(nocc_act_dim, nvir_act_dim, kit);
  BiggestContribs biggest_s2(2,10);
  // compute the S2 norm
  double s2norm = 0.0;
  for (j=nfzc; j<nocc; j++) {
    laj_ptr = &Laj[j*nvir];
    for (a=0; a<nvir_act; a++) {
      tmpval = 0.5**laj_ptr++/(evals[a+nocc]-evals[j]);
      s2.set_element(j-nfzc,a,tmpval);
      biggest_s2.insert(tmpval,j,a);
      s2norm += tmpval*tmpval;
      }
    }
  // compute the S2 matrix 1-norm
  double s2onenorm = 0.0;
  RefSCElementSumAbs sumabs = new SCElementSumAbs;
  RefSCElementOp genop = sumabs.pointer();
  for (a=0; a < nvir_act; a++) {
    sumabs->init();
    s2.get_column(a).element_op(genop);
    if (s2onenorm < sumabs->result()) s2onenorm = sumabs->result();
    }
  // compute the S2 matrix inf-norm
  double s2infnorm = 0.0;
  for (j=0; j < nocc_act; j++) {
    sumabs->init();
    s2.get_row(j).element_op(genop);
    if (s2infnorm < sumabs->result()) s2infnorm = sumabs->result();
    }
  // compute the S2 matrix 2-norm
  RefSymmSCMatrix s2s2(nocc_act_dim,kit);
  s2s2.assign(0.0);
  s2s2.accumulate_symmetric_product(s2);
  s2 = 0;
  double s2twonorm = sqrt(s2s2.eigvals().get_element(nocc_act-1));
  s2s2 = 0;
  // compute the S2 norm
  s2norm = sqrt(s2norm/(2*nocc_act));
  // print the norms
  cout <<node0
       <<indent<<scprintf("S2(ov) matrix 1-norm   = %12.8f", s2onenorm) <<endl
       <<indent<<scprintf("S2 matrix 2-norm       = %12.8f", s2twonorm) <<endl
       <<indent<<scprintf("S2(ov) matrix inf-norm = %12.8f", s2infnorm) <<endl
       <<indent<<scprintf("S2 norm                = %12.8f", s2norm) <<endl;
  if (biggest_s2.ncontrib()) {
    cout << node0 << endl
         << indent << "Largest S2 values (unique determinants):" << endl;
    }
  for (i=0; i<biggest_s2.ncontrib(); i++) {
    int i0 = orbital_map[biggest_s2.indices(i)[0]];
    int i1 = orbital_map[biggest_s2.indices(i)[1] + nocc];
    cout << node0 << indent << scprintf("  %2d %12.8f %2d %3s -> %2d %3s",
                                        i+1, biggest_s2.val(i),
                                        symorb_num_[i0]+1,
                                        ct.gamma(symorb_irrep_[i0]).symbol(),
                                        symorb_num_[i1]+1,
                                        ct.gamma(symorb_irrep_[i1]).symbol())
         << endl;
    }


  // Finish computation of Wab
  tim_enter("Pab and Wab");
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
  tim_exit("Pab and Wab");
  RefSCMatrix Wab_matrix(nvir_dim, nvir_dim, kit);
  Wab_matrix->assign(Wab); // Put elements of Wab into Wab_matrix
  free(Wab);

  // Update Wkj with contribution from Pkj
  tim_enter("Pkj and Wkj");
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
  tim_exit("Pkj and Wkj");

  /////////////////////////////////
  // Finish the computation of Laj
  /////////////////////////////////

  tim_enter("Laj");

  RefSCMatrix Cv(nbasis_dim, nvir_dim, kit); // virtual block of scf_vector
  RefSCMatrix Co(nbasis_dim, nocc_dim, kit); // occupied block of scf_vector
  for (p=0; p<nbasis; p++) {
    c_pq = scf_vector[p];
    for (q=0; q<nbasis; q++) {
      if (q<nocc) Co->set_element(p, q, *c_pq++);
      else Cv->set_element(p, q-nocc, *c_pq++);
      }
    }

  // Compute the density-like matrix Dmat_matrix
  RefSymmSCMatrix Pab_matrix(nvir_dim,kit);
  RefSymmSCMatrix Pkj_matrix(nocc_dim,kit);
  RefSCMatrix Dmat_matrix(nbasis_dim,nbasis_dim,kit);
  Pab_matrix->assign(Pab); // fill in elements of Pab_matrix from Pab
  free(Pab);
  Pkj_matrix->assign(Pkj); // fill in elements of Pkj_matrix from Pkj
  free(Pkj);
  Dmat_matrix = Cv*Pab_matrix*Cv.t() + Co*Pkj_matrix*Co.t();
  // We now have the density-like matrix Dmat_matrix

  // Compute the G matrix
  Dmat = new double[nbasis*nbasis];
  Dmat_matrix->convert(Dmat); // convert Dmat_matrix to Dmat (double*)

  RefSymmSCMatrix Gmat(nbasis_dim,kit);
  init_cs_gmat();
  tim_enter("make_gmat for Laj");
  make_cs_gmat(Gmat, Dmat);
  if (debug_ > 1) {
    Dmat_matrix.print("Dmat");
    Gmat.print("Gmat");
    }
  tim_exit("make_gmat for Laj");

  // Finish computation of Laj
  RefSCMatrix Laj_matrix(nocc_dim,nvir_dim,kit); // elements are ordered as j*nvir+a
  Laj_matrix->assign(Laj);
  if (debug_ > 1) Laj_matrix->print("Laj (first bit)");
  Laj_matrix = Laj_matrix - 2*Co.t()*Gmat*Cv;
  if (debug_ > 1) Laj_matrix->print("Laj (all of it)");
  Laj_matrix->convert(Laj);  // Put new Laj_matrix elements into Laj

  tim_exit("Laj");

  //////////////////////////////////////
  // Computation of Laj is now complete
  //////////////////////////////////////

  ////////////////////////////
  // Solve the CPHF equations
  ////////////////////////////
  RefSCMatrix Paj_matrix(nvir_dim, nocc_dim, kit);
  tim_enter("cphf");
  cs_cphf(scf_vector, Laj, evals, Paj_matrix);
  tim_exit("cphf");

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
  tim_enter("Pkj and Wkj");
  Dmat_matrix = Co*(Pkj_matrix*Co.t() + Paj_matrix.t()*Cv.t()) +
                Cv*(Paj_matrix*Co.t() + Pab_matrix*Cv.t());
  Dmat_matrix->convert(Dmat); // convert Dmat_matrix to Dmat (double*)
  tim_enter("make_gmat for Wkj");
  make_cs_gmat(Gmat, Dmat);
  tim_exit("make_gmat for Wkj");
  done_cs_gmat();
  RefSCMatrix Wkj_matrix(nocc_dim, nocc_dim, kit);
  Wkj_matrix->assign(Wkj);
  Wkj_matrix = Wkj_matrix - 2*Co.t()*Gmat*Co;
  free(Wkj);
  // Wkj is now complete - not as Wkj but as Wkj_matrix
  tim_exit("Pkj and Wkj");

  delete[] Dmat;

//scf_done_bounds(); // must be called after last call to make_gmat

  ////////////////////////////////////////////////////////////////
  // We now have the matrices Pkj_matrix, Paj_matrix, Pab_matrix,
  // Wkj_matrix, Waj_matrix, Wab_matrix and can compute the 
  // remaining contributions to the gradient
  ///////////////////////////////////////////////////////////////

  // Compute the second order correction to 
  // the density matrix and energy weighted
  // density matrix in the AO basis
  RefSCMatrix P2AO_matrix(nbasis_dim, nbasis_dim, kit);
  RefSCMatrix P2MO_matrix(nbasis_dim, nbasis_dim, kit);
  RefSCMatrix W2AO_matrix(nbasis_dim, nbasis_dim, kit);
  RefSCMatrix W2MO_matrix(nbasis_dim, nbasis_dim, kit);
  RefSCMatrix SCF_matrix(nbasis_dim, nbasis_dim, kit);
  for (i=0; i<nocc; i++) {
    for (j=0; j<nocc; j++) {
      P2MO_matrix->set_element(i,j,Pkj_matrix->get_element(i,j));
      W2MO_matrix->set_element(i,j,Wkj_matrix->get_element(i,j));
      SCF_matrix->set_element(i,j,Co->get_element(i,j));
      }
    for (j=nocc; j<nbasis; j++) {
      P2MO_matrix->set_element(i,j,Paj_matrix->get_element(j-nocc,i));
      W2MO_matrix->set_element(i,j,Waj_matrix->get_element(i,j-nocc));
      SCF_matrix->set_element(i,j,Cv->get_element(i,j-nocc));
      }
    }
  for (i=nocc; i<nbasis; i++) {
    for (j=0; j<nocc; j++) {
      P2MO_matrix->set_element(i,j,Paj_matrix->get_element(i-nocc,j));
      W2MO_matrix->set_element(i,j,Waj_matrix->get_element(j,i-nocc));
      SCF_matrix->set_element(i,j,Co->get_element(i,j));
      }
    for (j=nocc; j<nbasis; j++) {
      P2MO_matrix->set_element(i,j,Pab_matrix->get_element(i-nocc,j-nocc));
      W2MO_matrix->set_element(i,j,Wab_matrix->get_element(i-nocc,j-nocc));
      SCF_matrix->set_element(i,j,Cv->get_element(i,j-nocc));
      }
    }
  P2AO_matrix = SCF_matrix * P2MO_matrix * SCF_matrix.t();
  W2AO_matrix = SCF_matrix * W2MO_matrix * SCF_matrix.t();
//  P2AO_matrix = Co*(Pkj_matrix*Co.t() + Paj_matrix.t()*Cv.t()) +
//                Cv*(Paj_matrix*Co.t() + Pab_matrix*Cv.t());
//  W2AO_matrix = Co*(Wkj_matrix*Co.t() + Waj_matrix*Cv.t()) +
//                Cv*(Waj_matrix.t()*Co.t() + Wab_matrix*Cv.t());

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
  delete[] WHF;
  delete[] W2AO;

  if (debug_ > 1) {
    RefSCMatrix tmpmat(basis()->basisdim(), basis()->basisdim(), kit);
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
  tim_enter("sep 2PDM contrib.");
  s2pdm_contrib(intderbuf, PHF, P2AO, ginter);
  tim_exit("sep 2PDM contrib.");
  delete[] PHF;
  delete[] P2AO;

  // The separable 2PDM contribution to the gradient has now been
  // accumulated in ginter on node 0; add it to the total gradients
  if (me == 0) {
    accum_gradients(gradient, ginter, natom, 3);
    }
  // Print out the contribution to the gradient from sep. 2PDM
  if (debug_ && me == 0) {
    cout <<indent
         << "Contribution from separable 2PDM to MP2 gradient [au]:" << endl;
    for (i=0; i<natom; i++) {
      cout << indent << scprintf("%15.10lf  %15.10lf  %15.10lf\n",
                                 ginter[i][0], ginter[i][1], ginter[i][2]);
      }
    }

  // Done with two-electron integrals
  tbint_ = 0;
  tbintder_ = 0;

  /////////////////////////////////////////////////////////////
  // Compute the one-electron contribution to the MP2 gradient
  /////////////////////////////////////////////////////////////

  zero_gradients(ginter, natom, 3);
  tim_enter("hcore contrib.");
  hcore_cs_grad(PMP2, ginter);
  tim_exit("hcore contrib.");
  delete[] PMP2;
  // The hcore contribution to the gradient has now been accumulated
  // in ginter on node 0; add it to the total gradients
  if (me == 0) {
    accum_gradients(gradient, ginter, natom, 3);
    }
  // Print out the contribution to the gradient from hcore
  if (debug_ && me == 0) {
    cout << indent << "Contribution to MP2 gradient from hcore [au]:" << endl;
    for (i=0; i<natom; i++) {
      cout << indent << scprintf("%15.10lf  %15.10lf  %15.10lf\n",
                                 ginter[i][0], ginter[i][1], ginter[i][2]);
      }
    }

  zero_gradients(ginter, natom, 3);
  tim_enter("overlap contrib.");
  overlap_cs_grad(WMP2, ginter);
  tim_exit("overlap contrib.");
  delete[] WMP2;
  // The overlap contribution to the gradient has now been accumulated
  // in ginter on node 0; add it to the total gradients
  if (me == 0) {
    accum_gradients(gradient, ginter, natom, 3);
    }
  // Print out the overlap contribution to the gradient
  if (debug_ && me == 0) {
    cout << indent << "Overlap contribution to MP2 gradient [au]:" << endl;
    for (i=0; i<natom; i++) {
      cout << indent << scprintf("%15.10lf  %15.10lf  %15.10lf\n",
                                 ginter[i][0], ginter[i][1], ginter[i][2]);
      }
    }

  ////////////////////////////////////////////////////////
  // Compute the nuclear contribution to the MP2 gradient
  ////////////////////////////////////////////////////////

  if (me == 0) {
    zero_gradients(ginter, natom, 3);
    for (i=0; i<natom; i++) {
      molecule()->nuclear_repulsion_1der(i,ginter[i]);
      }
    accum_gradients(gradient, ginter, natom, 3);

    // Print out the nuclear contribution to the gradient
    if (debug_) {
      cout << indent
           << scprintf("Nuclear contribution to MP2 gradient [au]:") << endl;
      for (i=0; i<natom; i++) {
        cout << indent << scprintf("%15.10lf  %15.10lf  %15.10lf\n",
                                   ginter[i][0], ginter[i][1], ginter[i][2]);
        }
      }
    }


  ////////////////////////////////////////////////////////
  // The computation of the MP2 gradient is now complete;
  // print out the gradient
  ////////////////////////////////////////////////////////
  if (me == 0) {
    cout << endl << indent << "Total MP2 gradient [au]:" << endl;
    for (i=0; i<natom; i++) {
      cout << indent
           << scprintf("%15.10lf  %15.10lf  %15.10lf\n",
                       gradient[i][0], gradient[i][1], gradient[i][2]);
      }
    cout.flush();
    }

  msg_->bcast(gradient_dat, natom*3);
  RefSCVector gradientvec = matrixkit()->vector(moldim());
  gradientvec->assign(gradient_dat);
  set_gradient(gradientvec);

  delete[] gradient;
  delete[] gradient_dat;

  for (i=0; i<natom; i++) {
    delete[] ginter[i];
    }
  delete[] ginter;

  delete[] scf_vector;
  delete[] scf_vector_dat;
  delete[] evals;
  delete[] orbital_map;

  tim_exit("mp2-mem");
  }

///////////////////////////////////////////////////////////
// Compute the contribution to the MP2 gradient from hcore
///////////////////////////////////////////////////////////
void 
MBPT2::hcore_cs_grad(double *PMP2, double **ginter)
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
  double tmpval1;
  double gxyz[3];

  // Initialize 1e object
  RefOneBodyDerivInt obintder_ = integral()->hcore_deriv();
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

          for (l=0; l<3; l++) gxyz[l] = 0.0;

          index = 0;

          for (jj=0; jj<jsize; jj++) {
            for (kk=0; kk<ksize; kk++) {
              tmpval1 = PMP2[(j_offset + jj)*nbasis + k_offset + kk];
              for (m=0; m<3; m++) {
                gxyz[m] += oneebuf[index] * tmpval1;
                index++;
                } // exit m loop
              }   // exit kk loop
            }     // exit jj loop

          if (j != k) {
            for (l=0; l<3; l++) gxyz[l] *= 2.0; // off-diagonal blocks
            }

          for (l=0; l<3; l++) ginter[i][l] += gxyz[l];

          } // exit "if"
        }   // exit k loop
      }     // exit j loop
    }       // exit i loop

  /* Accumulate the nodes' intermediate gradients on node 0 */
  sum_gradients(msg_, ginter, molecule()->natom(), 3);
}


////////////////////////////////////////////////////
// Compute the overlap contribution to the gradient
////////////////////////////////////////////////////
void 
MBPT2::overlap_cs_grad(double *WMP2, double **ginter)
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
  double gxyz[3];

  // Initialize 1e object
  RefOneBodyDerivInt obintder_ = integral()->overlap_deriv();
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

          for (l=0; l<3; l++) gxyz[l] = 0.0;
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

              for (m=0; m<3; m++) {
                if (jj_index != kk_index) {
                  gxyz[m] += oneebuf[index] * (tmpval1 + tmpval2);
                  }
                else {
                  gxyz[m] += oneebuf[index] * tmpval1;
                  }
                index++;
                } // exit m loop
              }   // exit kk loop
            }     // exit jj loop

          for (l=0; l<3; l++) ginter[i][l] += gxyz[l];
          } // exit "if"
        }   // exit k loop
      }     // exit j loop
    }       // exit i loop
  
  /* Accumulate the nodes' intermediate gradients on node 0 */
  sum_gradients(msg_, ginter, molecule()->natom(), 3);
}


//////////////////////////////////////////////////////////////
// Compute (in the AO basis) the contribution to the gradient 
// from the separable part of the two particle density matrix
//////////////////////////////////////////////////////////////
void
MBPT2::s2pdm_contrib(const double *intderbuf, double *PHF,
                     double *P2AO, double **ginter)
{               

  int P, Q, R, S;
  int QP, SR;
  int p, q, r;
  int np, nq, nr, ns;
  int p_offset, q_offset, r_offset, s_offset;
  int bf1, bf2, bf3, bf4;
  int index;
  int xyz;
  int derset;
  int nshell = basis()->nshell();
  int nbasis = basis()->nbasis();
  int flags;
  int nproc = msg_->n();
  int me = msg_->me();

  double *grad_ptr1, *grad_ptr2;
  double tmpval;
  double *phf_pq, *phf_pr, *phf_ps, *phf_qr, *phf_qs, *phf_rs;
  double *p2ao_pq, *p2ao_pr, *p2ao_ps, *p2ao_qr, *p2ao_qs, *p2ao_rs;
  double k_QP, k_SR, k_QPSR; // factors needed since we loop over nonredund
                             // shell quartets but do redund integrals within
                             // shell quartets when applicable
  double gamma_factor; // factor multiplying integrals; needed because we
                       // loop over nonredund shell quarters but do redund
                       // integrals within shell quartets when applicable
  double *gammasym_pqrs; // symmetrized sep. 2PDM
  double *gammasym_ptr;
  const double *integral_ptr;

  int tol = (int) (-10.0/log10(2.0));  // discard erep derivatives smaller than 10^-10

  gammasym_pqrs = new double[nfuncmax*nfuncmax*nfuncmax*nfuncmax];

  DerivCenters der_centers;

  index = 0;

  tim_enter("PQRS loop");

  for (Q=0; Q<nshell; Q++) {
    nq = basis()->shell(Q).nfunction();
    q_offset = basis()->shell_to_function(Q);

    for (S=0; S<=Q; S++) {
      ns = basis()->shell(S).nfunction();
      s_offset = basis()->shell_to_function(S);

      for (R=0; R<=S; R++) {
        nr = basis()->shell(R).nfunction();
        r_offset = basis()->shell_to_function(R);
        k_SR = (R == S ? 0.5 : 1.0);
        SR = S*(S+1)/2 + R;

        for (P=0; P<=(S==Q ? R:Q); P++) {
          // If integral derivative is 0, skip to next P
          if (tbintder_->log2_shell_bound(P,Q,R,S) < tol) continue;

          index++;

          if (index%nproc == me) {
            np = basis()->shell(P).nfunction();
            p_offset = basis()->shell_to_function(P);
            k_QP = (P == Q ? 0.5 : 1.0);
            QP = Q*(Q+1)/2 + P;
            k_QPSR = (QP == SR ? 0.5 : 1.0);
            gamma_factor = k_QP*k_SR*k_QPSR;

            // Evaluate derivative integrals
            tbintder_->compute_shell(P,Q,R,S,der_centers);

            //////////////////////////////
            // Symmetrize sep. 2PDM
            //////////////////////////////
            gammasym_ptr = gammasym_pqrs;
            for (bf1=0; bf1<np; bf1++) {
              p = p_offset + bf1;
              phf_pq =  &PHF [p*nbasis + q_offset];
              p2ao_pq = &P2AO[p*nbasis + q_offset];

              for (bf2=0; bf2<nq; bf2++) {
                q = q_offset + bf2;
                phf_pr =  &PHF [p*nbasis + r_offset];
                p2ao_pr = &P2AO[p*nbasis + r_offset];
                phf_qr =  &PHF [q*nbasis + r_offset];
                p2ao_qr = &P2AO[q*nbasis + r_offset];

                for (bf3=0; bf3<nr; bf3++) {
                  r = r_offset + bf3;
                  phf_ps =  &PHF [p*nbasis + s_offset];
                  phf_qs =  &PHF [q*nbasis + s_offset];
                  phf_rs =  &PHF [r*nbasis + s_offset];
                  p2ao_ps = &P2AO[p*nbasis + s_offset];
                  p2ao_qs = &P2AO[q*nbasis + s_offset];
                  p2ao_rs = &P2AO[r*nbasis + s_offset];

                  for (bf4=0; bf4<ns; bf4++) {

                    *gammasym_ptr++ = gamma_factor*(
                                          4**phf_pq*(*phf_rs + *p2ao_rs)
                                        + 4**phf_rs**p2ao_pq
                                        - *phf_qs*(*phf_pr + *p2ao_pr)
                                        - *phf_qr*(*phf_ps + *p2ao_ps)
                                        - *phf_ps**p2ao_qr
                                        - *phf_pr**p2ao_qs);

                    phf_ps++;
                    phf_qs++;
                    phf_rs++;
                    p2ao_ps++;
                    p2ao_qs++;
                    p2ao_rs++;
                    } // exit bf4 loop
                  phf_pr++;
                  p2ao_pr++;
                  phf_qr++;
                  p2ao_qr++;
                  }   // exit bf3 loop
                phf_pq++;
                p2ao_pq++;
                }     // exit bf2 loop
              }       // exit bf1 loop

            ///////////////////////////////////////////////////////////
            // Contract symmetrized sep 2PDM with integral derivatives
            ///////////////////////////////////////////////////////////
            integral_ptr = intderbuf;
            for (derset=0; derset<der_centers.n(); derset++) {

              for (xyz=0; xyz<3; xyz++) {
                grad_ptr1 = &ginter[der_centers.atom(derset)][xyz];
                if (der_centers.has_omitted_center())
                  grad_ptr2 = &ginter[der_centers.omitted_atom()][xyz];

                gammasym_ptr = gammasym_pqrs;
                for (bf1=0; bf1<np; bf1++) {
                  for (bf2=0; bf2<nq; bf2++) {
                    for (bf3=0; bf3<nr; bf3++) {
                      for (bf4=0; bf4<ns; bf4++) {
                        tmpval = *integral_ptr++ * *gammasym_ptr++;
                        *grad_ptr1 += tmpval;
                        *grad_ptr2 -= tmpval;
                        } // exit bf4 loop
                      }   // exit bf3 loop
                    }     // exit bf2 loop
                  }       // exit bf1 loop
                }         // exit xyz loop
              }           // exit derset loop


            } // exit "if"
          }   // exit P loop
        }     // exit R loop
      }       // exit S loop
    }         // exit Q loop

  tim_exit("PQRS loop");

  delete[] gammasym_pqrs;

  // Accumulate intermediate gradients on node 0
  sum_gradients(msg_, ginter, molecule()->natom(), 3);

}

static void
sum_gradients(const RefMessageGrp& msg, double **f, int n1, int n2)
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
MBPT2::compute_cs_batchsize(int mem_static, int nocc_act)
{
  int index;
  int mem1, mem2, mem3;
  int mem_dyn;   // dynamic memory available
  distsize_t maxdyn;
  int tmp;
  int i, j;
  int ni;
  int nproc = msg_->n();

  ///////////////////////////////////////
  // the largest memory requirement will
  // either occur just before the end of
  // the 1. q.b.t. (mem1) or just before
  // the end of the i-batch loop (mem2)
  ///////////////////////////////////////

  mem_dyn = mem_alloc - mem_static;

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
  mem1 = sizeof(double)*(nij*(distsize_t)nbasis*(distsize_t)nbasis
                         + nbasis*(distsize_t)nvir);
  mem2 = sizeof(double)*(ni*nbasis*nfuncmax*(distsize_t)nfuncmax
                         + nij*(distsize_t)nbasis*(distsize_t)nbasis
                         + ni*(distsize_t)nbasis + nbasis*(distsize_t)nfuncmax
                         + nfuncmax*nfuncmax*nfuncmax*(distsize_t)nfuncmax);
  mem3 = sizeof(double)*(ni*nbasis*nfuncmax*(distsize_t)nfuncmax
                         + nij*(distsize_t)nbasis*(distsize_t)nbasis
                         + 2*(2 + nbasis*(distsize_t)nfuncmax));
  tmp = (mem2>mem3 ? mem2:mem3);
  maxdyn = (mem1>tmp ? mem1:tmp);

  return maxdyn;
}

////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
