//
// hsosv2lb.cc
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

#include <math.h>

#include <util/misc/regtime.h>
#include <util/misc/formio.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/mbpt/mbpt.h>
#include <chemistry/qc/mbpt/bzerofast.h>

using namespace std;
using namespace sc;

static void iqs(int *item,int *index,int left,int right);
static void iquicksort(int *item,int *index,int n);
static void findprocminmax(int *nbf, int nproc,
                           int *procmin, int *procmax, int *minbf, int *maxbf);
static void findshellmax(int *myshellsizes, int nRshell, int *shellmax, 
                         int *shellmaxindex);
static void expandintarray(int *&a, int dim);

void
MBPT2::compute_hsos_v2_lb()
{
  int i, j, k, l;
  int s1, s2;
  int a, b;
  int isocc, asocc;   // indices running over singly occupied orbitals
  int nfuncmax = basis()->max_nfunction_in_shell();
  int nvir;
  int nshell;
  int shellmax;
  int shellmaxindex;
  int nocc=0,ndocc=0,nsocc=0;
  int i_offset; 
  int npass, pass;
  int ni;
  int np, nq, nr, ns; 
  int P, Q, R, S;
  int p, q, r, s;
  int bf1, bf2, bf3, bf4;
  int bf3_offset;
  int nbfmoved;
  int nbfav;         // average number of r basis functions per node
  int minbf, maxbf;  // max/min number of (r) basis functions on a node
  int index;
  int compute_index;
  int col_index;
  int tmp_index;
  int dim_ij;
  int docc_index, socc_index, vir_index;
  int me;
  int nproc;
  int procmin, procmax; // processor with most/fewest basis functions
  int rest;
  int r_offset;
  int min;
  int iproc;
  int nRshell;
  int imyshell;
  int *myshells;       // the R indices processed by node me
  int *myshellsizes;   // sizes of the shells (after split) on node me 
  int *split_info;     // on each node: offset for each shell; -1 if shell not split 
  int *shellsize;      // size of each shell
  int *sorted_shells;  // sorted shell indices: large shells->small shells
  int *nbf;            // number of basis functions processed by each node
  int *proc;           // element k: processor which will process shell k
  int aoint_computed = 0; 
  double A, B, C, ni_top, max, ni_double; // variables used to compute ni
  double *evals_open;    // reordered scf eigenvalues
  const double *intbuf;        // 2-electron AO integral buffer
  double *trans_int1;    // partially transformed integrals
  double *trans_int2;    // partially transformed integrals
  double *trans_int3;    // partially transformed integrals
  double *trans_int4;    // fully transformed integrals
  double *trans_int4_tmp; // scratch array
  double *mo_int_do_so_vir=0;//mo integral (is|sa); i:d.o.,s:s.o.,a:vir
  double *mo_int_tmp=0;  // scratch array used in global summations
  double *socc_sum=0;    // sum of 2-el integrals involving only s.o.'s
  double *socc_sum_tmp=0;// scratch array
  double *iqrs, *iprs;
  double *iars_ptr;
  double iars;
  double iajr;
  double *iajr_ptr;
  double *iajb;
  double pqrs;
  double *c_qa;
  double *c_rb, *c_pi, *c_qi, *c_sj;
  double delta_ijab;
  double delta;
  double contrib1, contrib2;
  double ecorr_opt2=0,ecorr_opt1=0;
  double ecorr_zapt2;
  double ecorr_opt2_contrib=0, ecorr_zapt2_contrib=0;
  double escf;
  double eopt2,eopt1,ezapt2;
  double tol;          // log2 of the erep tolerance (erep < 2^tol => discard)

  me = msg_->me();
 
  ExEnv::out0() << indent << "Just entered OPT2 program (opt2v2lb)" << endl;

  tol = (int) (-10.0/log10(2.0));  // discard ereps smaller than 10^-10

  nproc = msg_->n();
  ExEnv::out0() << indent << "nproc = " << nproc << endl;

  ndocc = nsocc = 0;
  const double epsilon = 1.0e-4;
  for (i=0; i<oso_dimension()->n(); i++) {
    if      (reference_->occupation(i) >= 2.0 - epsilon) ndocc++;
    else if (reference_->occupation(i) >= 1.0 - epsilon) nsocc++;
    }

  // Do a few preliminary tests to make sure the desired calculation
  // can be done (and appears to be meaningful!)

  if (ndocc == 0 && nsocc == 0) {
    ExEnv::err0() << "There are no occupied orbitals; program exiting" << endl;
    abort();
    }

  if (nfzc > ndocc) {
    ExEnv::err0()
         << "The number of frozen core orbitals exceeds the number" << endl
         << "of doubly occupied orbitals; program exiting" << endl;
    abort();
    }

  if (nfzv > noso - ndocc - nsocc) {
    ExEnv::err0()
         << "The number of frozen virtual orbitals exceeds the number" << endl
         << "of unoccupied orbitals; program exiting" << endl;
    abort();
    }

  ndocc = ndocc - nfzc;
  // nvir = # of unocc. orb. + # of s.o. orb. - # of frozen virt. orb.
  nvir  = noso - ndocc - nfzc - nfzv; 
  // nocc = # of d.o. orb. + # of s.o. orb - # of frozen d.o. orb.
  nocc  = ndocc + nsocc;
  nshell = basis()->nshell();

  // Allocate storage for some arrays used for keeping track of which R
  // indices are processed by each node                                
  shellsize = (int*) malloc(nshell*sizeof(int));
  sorted_shells = (int*) malloc(nshell*sizeof(int));
  nbf = (int*) malloc(nproc*sizeof(int));
  proc = (int*) malloc(nshell*sizeof(int));


  ///////////////////////////////////////////////////////
  // Begin distributing R shells between nodes so all
  // nodes get ca. the same number of r basis functions
  ///////////////////////////////////////////////////////

  // Compute the size of each shell
  for (i=0; i<nshell; i++) {
    shellsize[i] = basis()->shell(i).nfunction();
    }

  // Do an index sort (large -> small) of shellsize to form sorted_shells
  iquicksort(shellsize,sorted_shells,nshell);

  // Initialize nbf
  for (i=0; i<nproc; i++) nbf[i] = 0;

  for (i=0; i<nshell; i++) {
    min = nbf[0];
    iproc = 0;
    for (j=1; j<nproc; j++) {
      if (nbf[j] < min) {
        iproc = j;
        min = nbf[j];
        }
      }
    proc[sorted_shells[i]] = iproc;
    nbf[iproc] += shellsize[sorted_shells[i]];
    }
  if (me == 0) {
    ExEnv::out0() << indent << "Distribution of basis functions between nodes:" << endl;
    for (i=0; i<nproc; i++) {
      if (i%12 == 0) ExEnv::out0() << indent;
      ExEnv::out0() << scprintf(" %4i",nbf[i]);
      if ((i+1)%12 == 0) ExEnv::out0() << endl;
      }
    ExEnv::out0() << endl;
    }

  // Determine which shells are to be processed by node me
  nRshell = 0;
  for (i=0; i<nshell; i++) {
    if (proc[i] == me) nRshell++;
    }
  myshells = (int*) malloc(nRshell*sizeof(int));
  imyshell = 0;
  for (i=0; i<nshell; i++) {
    if (proc[i] == me) {
      myshells[imyshell] = i;
      imyshell++;
      }
    }

  /////////////////////////////////////////////////////////////
  // End of preliminary distribution of R shells between nodes
  /////////////////////////////////////////////////////////////

  // Compute the average number of basis functions per node
  nbfav = nbasis/nproc;
  if (nbasis%nproc) nbfav++;

  myshellsizes = (int*) malloc(nRshell*sizeof(int));
  split_info   = (int*) malloc(nRshell*sizeof(int));
  for (j=0; j<nRshell; j++) {
    myshellsizes[j] = basis()->shell(myshells[j]).nfunction();
    split_info[j] = -1;
    }

  // Find the processor with the most/fewest basis functions
  findprocminmax(nbf,nproc,&procmin,&procmax,&minbf,&maxbf);
  if (maxbf > nbfav) {
    ExEnv::out0() << indent << "Redistributing basis functions" << endl;
    }

  while (maxbf > nbfav) {
    msg_->sync();
    if (me == procmax) {

      findshellmax(myshellsizes, nRshell, &shellmax, &shellmaxindex);
      nbfmoved = 0;
      while (maxbf>nbfav && minbf<nbfav && shellmax>1) {
        shellmax--;
        nbfmoved++;
        maxbf--;
        minbf++;
        }
      myshellsizes[shellmaxindex] = shellmax;
      if (split_info[shellmaxindex] == -1) split_info[shellmaxindex] = 0;
      shellmax += nbfmoved;

      // Send nbfmoved from procmax to all other nodes
      msg_->bcast(nbfmoved,procmax);

      // Send variables to node procmin
      msg_->send(procmin,&myshells[shellmaxindex],1);
      msg_->send(procmin,&shellmax,1);

      }
    else {
      // Receive nbfmoved from procmax
      msg_->bcast(nbfmoved,procmax);
      }

    nbf[procmax] -= nbfmoved;

    if (me == procmin) {
      expandintarray(myshellsizes,nRshell);
      expandintarray(myshells,nRshell);
      expandintarray(split_info,nRshell);
      nRshell++;
      myshellsizes[nRshell-1] = nbfmoved;
      msg_->recv(procmax,&myshells[nRshell-1],1);
      msg_->recv(procmax,&split_info[nRshell-1],1);
      split_info[nRshell-1] -= myshellsizes[nRshell-1];
      }

    nbf[procmin] += nbfmoved;
    msg_->sync();
    findprocminmax(nbf,nproc,&procmin,&procmax,&minbf,&maxbf);

    }

  if (me == 0) {
    ExEnv::out0() << indent
         << "New distribution of basis functions between nodes:" << endl;
    for (i=0; i<nproc; i++) {
      if (i%12 == 0) ExEnv::out0() << indent;
      ExEnv::out0() << scprintf(" %4i",nbf[i]);
      if ((i+1)%12 == 0) ExEnv::out0() << endl;
      }
    ExEnv::out0() << endl;
    }


  //////////////////////////////////////////////////////////
  // End of distribution of R shells and r basis functions
  //////////////////////////////////////////////////////////

  // Compute batch size ni for opt2 loops;
  // need to store the following arrays of type double : trans_int1-4,
  // trans_int4_tmp, scf_vector, evals_open, socc_sum, socc_sum_tmp,
  // mo_int_do_so_vir, mo_int_tmp,
  // and the following arrays of type int: myshells, shellsize,
  // sorted_shells, nbf, and proc
    A = -0.5*sizeof(double)*nbf[me]*nvir;
    B = sizeof(double)*(nfuncmax*nfuncmax*nbasis + nvir + nocc*nbf[me]*nvir
                        + nbf[me]*nvir*0.5);
    C = sizeof(double)*(2*nvir*nvir + (nbasis+1)*(nvir+nocc) + 2*nsocc
                        + 2*ndocc*nsocc*(nvir-nsocc))
        + sizeof(int)*(3*nshell + nproc + nRshell);
    ni_top = -B/(2*A);
    max = A*ni_top*ni_top + B*ni_top +C;
    if (max <= mem_alloc) {
      ni = nocc;
      }
    else {
      ni_double = (-B + sqrt((double)(B*B - 4*A*(C-mem_alloc))))/(2*A);
      ni = (int) ni_double;
      if (ni > nocc) ni = nocc;
      max = mem_alloc;
      }

  size_t mem_remaining = mem_alloc - (size_t)max;

  // Set ni equal to the smallest batch size for any node
  msg_->min(ni);
  msg_->bcast(ni);

  if (ni < nsocc) {
    ExEnv::err0() << "Not enough memory allocated" << endl;
    abort();
    }

  if (ni < 1) {     // this applies only to a closed shell case
    ExEnv::err0() << "Not enough memory allocated" << endl;
    abort();
    }

  ExEnv::out0() << indent << "Computed batchsize: " << ni << endl;

  if (nocc == ni) {
    npass = 1;
    rest = 0;
    }
  else {
    rest = nocc%ni;
    npass = (nocc - rest)/ni + 1;
    if (rest == 0) npass--;
    }

  if (me == 0) {
    ExEnv::out0() << indent << " npass  rest  nbasis  nshell  nfuncmax"
                     "  ndocc  nsocc  nvir  nfzc  nfzv" << endl;
    ExEnv::out0() << indent
         << scprintf("  %-4i   %-3i   %-5i    %-4i     %-3i"
                     "     %-3i    %-3i    %-3i    %-3i   %-3i\n",
             npass,rest,nbasis,nshell,nfuncmax,ndocc,nsocc,nvir,nfzc,nfzv);
    ExEnv::out0() << indent
         << scprintf("Using %i bytes of memory",mem_alloc) << endl;
    }

  //////////////////////
  // Test that ni is OK
  //////////////////////
  if (me == 0) {
    ExEnv::out0() << indent
         << scprintf("Memory allocated: %i", mem_alloc) << endl;
    ExEnv::out0() << indent
         << scprintf("Memory used     : %lf", A*ni*ni+B*ni+C) << endl;
    if (A*ni*ni + B*ni +C > mem_alloc) {
      ExEnv::err0() << "Problems with memory allocation: "
           << "Using more memory than allocated" << endl;
      abort();
      }
    }

  //////////////////////////////////////////////////////////////////
  // The scf vector might be distributed between the nodes,
  // but for OPT2 each node needs its own copy of the vector;
  // therefore, put a copy of the scf vector on each node;
  // while doing this, duplicate columns corresponding to singly
  // occupied orbitals and order columns as [socc docc socc unocc]
  // Also rearrange scf eigenvalues as [socc docc socc unocc]
  // want socc first to get the socc's in the first batch
  // (need socc's to compute energy denominators - see
  // socc_sum comment below)
  /////////////////////////////////////////////////////////
  evals_open = (double*) malloc((noso+nsocc-nfzc-nfzv)*sizeof(double));

  RefDiagSCMatrix occ;
  RefDiagSCMatrix evals;
  RefSCMatrix Scf_Vec;
  eigen(evals, Scf_Vec, occ);

  if (debug_) {
    evals.print("eigenvalues");
    Scf_Vec.print("eigenvectors");
    }

  double *scf_vectort_dat = new double[nbasis*noso];
  Scf_Vec->convert(scf_vectort_dat);

  double** scf_vectort = new double*[nocc + nvir];

  int idoc = 0, ivir = 0, isoc = 0;
  for (i=nfzc; i<noso-nfzv; i++) {
    if (occ(i) >= 2.0 - epsilon) {
      evals_open[idoc+nsocc] = evals(i);
      scf_vectort[idoc+nsocc] = &scf_vectort_dat[i*nbasis];
      idoc++;
      }
    else if (occ(i) >= 1.0 - epsilon) {
      evals_open[isoc] = evals(i);
      scf_vectort[isoc] = &scf_vectort_dat[i*nbasis];
      evals_open[isoc+nocc] = evals(i);
      scf_vectort[isoc+nocc] = &scf_vectort_dat[i*nbasis];
      isoc++;
      }
    else {
      if (ivir < nvir) {
        evals_open[ivir+nocc+nsocc] = evals(i);
        scf_vectort[ivir+nocc+nsocc] = &scf_vectort_dat[i*nbasis];
        }
      ivir++;
      }
    }
  // need the transpose of the vector
  double **scf_vector = new double*[nbasis];
  double *scf_vector_dat = new double[(nocc+nvir)*nbasis];
  for (i=0; i<nbasis; i++) {
    scf_vector[i] = &scf_vector_dat[(nocc+nvir)*i];
    for (j=0; j<nocc+nvir; j++) {
      scf_vector[i][j] = scf_vectort[j][i];
      }
    }
  delete[] scf_vectort;
  delete[] scf_vectort_dat;

  ////////////////////////////////////////
  // Allocate storage for various arrays
  ////////////////////////////////////////

  dim_ij = nocc*ni - ni*(ni - 1)/2;

  trans_int1 = (double*) malloc(nfuncmax*nfuncmax*nbasis*ni*sizeof(double));
  trans_int2 = (double*) malloc(nvir*ni*sizeof(double));
  trans_int3 = (double*) malloc(nbf[me]*nvir*dim_ij*sizeof(double));
  trans_int4 = (double*) malloc(nvir*nvir*sizeof(double));
  trans_int4_tmp = (double*) malloc(nvir*nvir*sizeof(double));
  if (nsocc) socc_sum = (double*) malloc(nsocc*sizeof(double));
  if (nsocc) socc_sum_tmp = (double*) malloc(nsocc*sizeof(double));
  if (nsocc) mo_int_do_so_vir = 
                     (double*) malloc(ndocc*nsocc*(nvir-nsocc)*sizeof(double));
  if (nsocc) mo_int_tmp = 
                     (double*) malloc(ndocc*nsocc*(nvir-nsocc)*sizeof(double));

  if (nsocc) bzerofast(mo_int_do_so_vir,ndocc*nsocc*(nvir-nsocc));

  // create the integrals object
  integral()->set_storage(mem_remaining);
  tbint_ = integral()->electron_repulsion();
  intbuf = tbint_->buffer();

 /////////////////////////////////////
 //  Begin opt2 loops
 /////////////////////////////////////


  Timer tim;

  for (pass=0; pass<npass; pass++) {
    i_offset = pass*ni;  
    if ((pass == npass - 1) && (rest != 0)) ni = rest;

    r_offset = 0;
    bzerofast(trans_int3,nbf[me]*nvir*dim_ij);

    tim.enter("RS loop");

    for (imyshell=0; imyshell<nRshell; imyshell++) {

      R = myshells[imyshell];
      nr = myshellsizes[imyshell];

      for (S = 0; S < nshell; S++) {
        ns = basis()->shell(S).nfunction();
        tim.enter("bzerofast trans_int1");
        bzerofast(trans_int1,nfuncmax*nfuncmax*nbasis*ni);
        tim.exit("bzerofast trans_int1");

        tim.enter("PQ loop");

        for (P = 0; P < nshell; P++) {
          np = basis()->shell(P).nfunction();

          for (Q = 0; Q <= P; Q++) {
            if (tbint_->log2_shell_bound(P,Q,R,S) < tol) {
              continue;                          // skip ereps less than tol
              }

            aoint_computed++;

            nq = basis()->shell(Q).nfunction();

            tim.enter("erep");
            tbint_->compute_shell(P,Q,R,S);
            tim.exit("erep");

            tim.enter("1. quart. tr.");

            for (bf1 = 0; bf1 < np; bf1++) {
              p = basis()->shell_to_function(P) + bf1;
 
              for (bf2 = 0; bf2 < nq; bf2++) {
                q = basis()->shell_to_function(Q) + bf2;
                if (q > p) {
                  // if q > p: want to skip the loops over bf3-4
                  // and larger bf2 values, so increment bf1 by 1
                  // ("break")
                  break;
                  }

                for (bf3 = 0; bf3 < nr; bf3++) {
                  bf3_offset = 0;
                  if (split_info[imyshell] != -1) bf3_offset = split_info[imyshell];

                  for (bf4 = 0; bf4 < ns; bf4++) {

                    index = bf4 + ns*(bf3+bf3_offset + 
                               basis()->shell(R).nfunction()*(bf2 + nq*bf1));

                    if (fabs(intbuf[index]) > 1.0e-15) {
                      pqrs = intbuf[index];

                      iqrs = &trans_int1[((bf4*nr + bf3)*nbasis + q)*ni];
                      iprs = &trans_int1[((bf4*nr + bf3)*nbasis + p)*ni];
                    
                      if (p == q) pqrs *= 0.5;
                      col_index = i_offset;
                      c_pi = &scf_vector[p][col_index];
                      c_qi = &scf_vector[q][col_index];

                      for (i=ni; i; i--) {
                        *iqrs++ += pqrs * *c_pi++;
                        *iprs++ += pqrs * *c_qi++;
                        }
                      }
                    }   // exit bf4 loop
                  }     // exit bf3 loop
                }       // exit bf2 loop
              }         // exit bf1 loop
            tim.exit("1. quart. tr.");
            }           // exit Q loop
          }             // exit P loop
        tim.exit("PQ loop");

        // Begin second and third quarter transformations

        for (bf3 = 0; bf3 < nr; bf3++) {
          r = r_offset + bf3;

          for (bf4 = 0; bf4 < ns; bf4++) {
            s = basis()->shell_to_function(S) + bf4;

            tim.enter("bzerofast trans_int2");
            bzerofast(trans_int2,nvir*ni);
            tim.exit("bzerofast trans_int2");

            tim.enter("2. quart. tr.");

            for (q = 0; q < nbasis; q++) {
              iars_ptr = trans_int2;
              iqrs = &trans_int1[((bf4*nr + bf3)*nbasis + q)*ni];
              c_qa = &scf_vector[q][nocc];

              for (a = 0; a < nvir; a++) {

                for (i=ni; i; i--) {
                  *iars_ptr++ += *c_qa * *iqrs++;
                  }

                iqrs -= ni;
                c_qa++;
                }
              }             // exit q loop
            tim.exit("2. quart. tr.");

            // Begin third quarter transformation

            tim.enter("3. quart. tr.");

            for (i=0; i<ni; i++) {
              tmp_index = i*(i+1)/2 + i*i_offset;

              for (a=0; a<nvir; a++) {
                iars = trans_int2[a*ni + i];
                c_sj = scf_vector[s];
                iajr_ptr = &trans_int3[tmp_index + dim_ij*(a + nvir*r)];

                for (j=0; j<=i+i_offset; j++) {
                  *iajr_ptr++ += *c_sj++ * iars;
                  }
                }
              }   // exit i loop
            tim.exit("3. quart. tr.");

            } // exit bf4 loop
          }   // exit bf3 loop

        }         // exit S loop
      r_offset += nr;
      }           // exit R loop
    tim.exit("RS loop");

    // Begin fourth quarter transformation;
    // first tansform integrals with only s.o. indices;
    // these integrals are needed to compute the denominators
    // in the various terms contributing to the correlation energy
    // and must all be computed in the first pass;
    // the integrals are summed into the array socc_sum:
    // socc_sum[isocc] = sum over asocc of (isocc asocc|asocc isocc)
    // (isocc, asocc = s.o. and the sum over asocc runs over all s.o.'s)
    // the individual integrals are not saved here, only the sums are kept

    if (pass == 0) {
      tim.enter("4. quart. tr.");
      if (nsocc) bzerofast(socc_sum,nsocc);
      for (isocc=0; isocc<nsocc; isocc++) {

        index = 0;
        for (i=0; i<nRshell; i++) {
          for (j=0; j<myshellsizes[i]; j++) {
            r = basis()->shell_to_function(myshells[i]) + j;
            if (split_info[i] != -1) r += split_info[i];

            for (asocc=0; asocc<nsocc; asocc++) {
              socc_sum[isocc] += scf_vector[r][nocc+asocc]*
                             trans_int3[isocc*(isocc+1)/2 + isocc*i_offset
                                       + isocc + dim_ij*(asocc + nvir*index)];
              }
            index++;
            }
          }
        }       // exit i loop

      tim.exit("4. quart. tr.");

      // Sum socc_sum contributions from each node (only if nsocc > 0
      // since gop1 will fail if nsocc = 0)
      if (nsocc > 0) {
        tim.enter("global sum socc_sum");
        msg_->sum(socc_sum,nsocc,socc_sum_tmp);
        tim.exit("global sum socc_sum");
        }

      } 

    // Now we have all the sums of integrals involving s.o.'s (socc_sum);
    // begin fourth quarter transformation for all integrals (including
    // integrals with only s.o. indices); use restriction j <= (i_offset+i)
    // to save flops

    compute_index = 0;

    for (i=0; i<ni; i++) {

      for (j=0; j <= (i_offset+i); j++) {

        tim.enter("4. quart. tr.");

        bzerofast(trans_int4,nvir*nvir);

        index = 0;
        for (k=0; k<nRshell; k++) {
          for (l=0; l<myshellsizes[k]; l++) {
            r = basis()->shell_to_function(myshells[k]) + l;
            if (split_info[k] != -1) r += split_info[k];

            for (a=0; a<nvir; a++) {
              iajb = &trans_int4[a*nvir];
              iajr = trans_int3[i*(i+1)/2 + i*i_offset + j + dim_ij*(a+nvir*index)];
              c_rb = &scf_vector[r][nocc];

              for (b=0; b<nvir; b++) {
                *iajb++ += *c_rb++ * iajr;
                }
              }
            index++;
            }
          }   // end of k loop

        tim.exit("4. quart. tr.");

        tim.enter("global sum trans_int4");
        msg_->sum(trans_int4,nvir*nvir,trans_int4_tmp);
        tim.exit("global sum trans_int4");

        // We now have the fully transformed integrals (ia|jb)
        // for one i, one j (j <= i_offset+i), and all a and b;
        // compute contribution to the OPT1 and OPT2 correlation
        // energies; use restriction b <= a to save flops

        tim.enter("compute ecorr");

        for (a=0; a<nvir; a++) {
          for (b=0; b<=a; b++) {
            compute_index++;
            if (compute_index%nproc != me) continue;

            docc_index = ((i_offset+i) >= nsocc && (i_offset+i) < nocc) 
                        + (j >= nsocc && j < nocc);
            socc_index = ((i_offset+i)<nsocc)+(j<nsocc)+(a<nsocc)+(b<nsocc);
            vir_index = (a >= nsocc) + (b >= nsocc);

            if (socc_index >= 3) continue; // skip to next b value
 
            delta_ijab = evals_open[i_offset+i] + evals_open[j] 
                       - evals_open[nocc+a] - evals_open[nocc+b];
            
            // Determine integral type and compute energy contribution
            if (docc_index == 2 && vir_index == 2) {
              if (i_offset+i == j && a == b) {
                contrib1 = trans_int4[a*nvir + b]*trans_int4[a*nvir + b];
                ecorr_opt2 += contrib1/delta_ijab;
                ecorr_opt1 += contrib1/delta_ijab;
                }
              else if (i_offset+i == j || a == b) {
                contrib1 = trans_int4[a*nvir + b]*trans_int4[a*nvir + b];
                ecorr_opt2 += 2*contrib1/delta_ijab;
                ecorr_opt1 += 2*contrib1/delta_ijab;
                }
              else {
                contrib1 = trans_int4[a*nvir + b];
                contrib2 = trans_int4[b*nvir + a];
                ecorr_opt2 += 4*(contrib1*contrib1 + contrib2*contrib2
                             - contrib1*contrib2)/delta_ijab;
                ecorr_opt1 += 4*(contrib1*contrib1 + contrib2*contrib2
                             - contrib1*contrib2)/delta_ijab;
                }
              }
            else if (docc_index == 2 && socc_index == 2) {
              contrib1 = (trans_int4[a*nvir + b] - trans_int4[b*nvir + a])*
                         (trans_int4[a*nvir + b] - trans_int4[b*nvir + a]);
              ecorr_opt2 += contrib1/
                           (delta_ijab - 0.5*(socc_sum[a]+socc_sum[b]));
              ecorr_opt1 += contrib1/delta_ijab;
              }
            else if (socc_index == 2 && vir_index == 2) {
              contrib1 = (trans_int4[a*nvir + b] - trans_int4[b*nvir + a])*
                         (trans_int4[a*nvir + b] - trans_int4[b*nvir + a]);
              ecorr_opt2 += contrib1/
                          (delta_ijab - 0.5*(socc_sum[i_offset+i]+socc_sum[j]));
              ecorr_opt1 += contrib1/delta_ijab;
              }
            else if (docc_index == 2 && socc_index == 1 && vir_index == 1) {
              if (i_offset+i == j) {
                contrib1 = trans_int4[a*nvir + b]*trans_int4[a*nvir + b];
                ecorr_opt2 += contrib1/(delta_ijab - 0.5*socc_sum[b]);
                ecorr_opt1 += contrib1/delta_ijab;
                }
              else {
                contrib1 = trans_int4[a*nvir + b];
                contrib2 = trans_int4[b*nvir + a];
                ecorr_opt2 += 2*(contrib1*contrib1 + contrib2*contrib2 
                            - contrib1*contrib2)/(delta_ijab - 0.5*socc_sum[b]);
                ecorr_opt1 += 2*(contrib1*contrib1 + contrib2*contrib2 
                             - contrib1*contrib2)/delta_ijab;
                }
              }
            else if (docc_index == 1 && socc_index == 2 && vir_index == 1) {
              contrib1 = trans_int4[b*nvir+a]*trans_int4[b*nvir+a];
              if (j == b) {
                // To compute the energy contribution from an integral of the
                // type (is1|s1a) (i=d.o., s1=s.o., a=unocc.), we need the
                // (is|sa) integrals for all s=s.o.; these integrals are
                // therefore stored here in the array mo_int_do_so_vir, and
                // the energy contribution is computed after exiting the loop
                // over i-batches (pass)
                mo_int_do_so_vir[a-nsocc + (nvir-nsocc)*
                                (i_offset+i-nsocc + ndocc*b)] =
                                trans_int4[b*nvir + a];
                ecorr_opt2_contrib += 1.5*contrib1/delta_ijab;
                ecorr_opt1         += 1.5*contrib1/delta_ijab;
                ecorr_zapt2_contrib += contrib1/
                                (delta_ijab - 0.5*(socc_sum[j]+socc_sum[b]))
                           + 0.5*contrib1/delta_ijab;
                }
              else {
                ecorr_opt2 += contrib1/
                             (delta_ijab - 0.5*(socc_sum[j] + socc_sum[b]));
                ecorr_opt1 += contrib1/delta_ijab;
                }
              }
            else if (docc_index == 1 && socc_index == 1 && vir_index == 2) {
              if (a == b) {
                contrib1 = trans_int4[a*nvir + b]*trans_int4[a*nvir + b];
                ecorr_opt2 += contrib1/(delta_ijab - 0.5*socc_sum[j]);
                ecorr_opt1 += contrib1/delta_ijab;
                }
              else {
                contrib1 = trans_int4[a*nvir + b];
                contrib2 = trans_int4[b*nvir + a];
                ecorr_opt2 += 2*(contrib1*contrib1 + contrib2*contrib2
                            - contrib1*contrib2)/(delta_ijab - 0.5*socc_sum[j]);
                ecorr_opt1 += 2*(contrib1*contrib1 + contrib2*contrib2
                             - contrib1*contrib2)/delta_ijab;
                }
              }
            }   // exit b loop
          }     // exit a loop
        tim.exit("compute ecorr");
        }       // exit j loop
      }         // exit i loop

    if (nsocc == 0 && npass > 1 && pass < npass - 1) {
      double passe = ecorr_opt2;
      msg_->sum(passe);
      ExEnv::out0() << indent
           << "Partial correlation energy for pass " << pass << ":" << endl;
      ExEnv::out0() << indent
           << scprintf("  restart_ecorr        = %14.10f", passe)
           << endl;
      ExEnv::out0() << indent
           << scprintf("  restart_orbital_v2lb = %d", ((pass+1) * ni))
           << endl;
      }
    }           // exit loop over i-batches (pass)



  // Compute contribution from excitations of the type is1 -> s1a where
  // i=d.o., s1=s.o. and a=unocc; single excitations of the type i -> a,
  // where i and a have the same spin, contribute to this term;
  // (Brillouin's theorem not satisfied for ROHF wave functions);
  // do this only if nsocc > 0 since gop1 will fail otherwise

  tim.enter("compute ecorr");

  if (nsocc > 0) {
    tim.enter("global sum mo_int_do_so_vir");
    msg_->sum(mo_int_do_so_vir,ndocc*nsocc*(nvir-nsocc),mo_int_tmp);
    tim.exit("global sum mo_int_do_so_vir");
    }

  // Add extra contribution for triplet and higher spin multiplicities
  // contribution = sum over s1 and s2<s1 of (is1|s1a)*(is2|s2a)/delta

  if (me == 0 && nsocc) {
    for (i=0; i<ndocc; i++) {

      for (a=0; a<nvir-nsocc; a++) {
        delta = evals_open[nsocc+i] - evals_open[nocc+nsocc+a];

        for (s1=0; s1<nsocc; s1++) {

          for (s2=0; s2<s1; s2++) {
            contrib1 = mo_int_do_so_vir[a + (nvir-nsocc)*(i + ndocc*s1)]*
                  mo_int_do_so_vir[a + (nvir-nsocc)*(i + ndocc*s2)]/delta;
            ecorr_opt2 += contrib1;
            ecorr_opt1 += contrib1;
            }
          }
        }     // exit a loop
      }       // exit i loop
    }

  tim.exit("compute ecorr");

  ecorr_zapt2 = ecorr_opt2 + ecorr_zapt2_contrib;
  ecorr_opt2 += ecorr_opt2_contrib;
  msg_->sum(ecorr_opt1);
  msg_->sum(ecorr_opt2);
  msg_->sum(ecorr_zapt2);
  msg_->sum(aoint_computed);

  escf = reference_->energy();
  hf_energy_ = escf;

  if (me == 0) {
    eopt2 = escf + ecorr_opt2;
    eopt1 = escf + ecorr_opt1;
    ezapt2 = escf + ecorr_zapt2;

    // Print out various energies etc.
    ExEnv::out0() << indent
         << "Number of shell quartets for which AO integrals would" << endl
         << indent << "have been computed without bounds checking: "
         << npass*nshell*nshell*(nshell+1)*(nshell+1)/2 << endl;
    ExEnv::out0() << indent
         << "Number of shell quartets for which AO integrals" << endl
         << indent << "were computed: " << aoint_computed << endl;
            
    ExEnv::out0() << indent
         << scprintf("ROHF energy [au]:                  %17.12lf\n", escf);
    ExEnv::out0() << indent
         << scprintf("OPT1 energy [au]:                  %17.12lf\n", eopt1);
    ExEnv::out0() << indent
         << scprintf("OPT2 second order correction [au]: %17.12lf\n",
                     ecorr_opt2);
    ExEnv::out0() << indent
         << scprintf("OPT2 energy [au]:                  %17.12lf\n", eopt2);
    ExEnv::out0() << indent
         << scprintf("ZAPT2 correlation energy [au]:     %17.12lf\n",
                     ecorr_zapt2);
    ExEnv::out0() << indent
         << scprintf("ZAPT2 energy [au]:                 %17.12lf\n", ezapt2);
    ExEnv::out0().flush();
    }

  msg_->bcast(eopt1);
  msg_->bcast(eopt2);
  msg_->bcast(ezapt2);

  if (method_ == "opt1") {
    set_energy(eopt1);
    set_actual_value_accuracy(reference_->actual_value_accuracy()
                              *ref_to_mp2_acc);
    }
  else if (method_ == "opt2") {
    set_energy(eopt2);
    set_actual_value_accuracy(reference_->actual_value_accuracy()
                              *ref_to_mp2_acc);
    }
  else if (nsocc == 0 && method_ == "mp") {
    set_energy(ezapt2);
    set_actual_value_accuracy(reference_->actual_value_accuracy()
                              *ref_to_mp2_acc);
    }
  else {
    if (method_ != "zapt") {
      ExEnv::out0() << indent
           << "MBPT2: bad method: " << method_ << ", using zapt" << endl;
      }
    set_energy(ezapt2);
    set_actual_value_accuracy(reference_->actual_value_accuracy()
                              *ref_to_mp2_acc);
    }

  free(trans_int1);
  free(trans_int2);
  free(trans_int3);
  free(trans_int4);
  free(trans_int4_tmp);
  if (nsocc) free(socc_sum);
  if (nsocc) free(socc_sum_tmp);
  if (nsocc) free(mo_int_do_so_vir);
  if (nsocc) free(mo_int_tmp);
  free(evals_open);
  free(myshells);
  free(shellsize);
  free (myshellsizes);
  free (split_info);
  free(sorted_shells);
  free(nbf);
  free(proc);

  delete[] scf_vector;
  delete[] scf_vector_dat;

  }

/////////////////////////////////////////////////////////////////
// Function iquicksort performs a quick sort (larger -> smaller) 
// of the integer data in item by the integer indices in index;
// data in item remain unchanged
/////////////////////////////////////////////////////////////////
static void
iquicksort(int *item,int *index,int n)
{
  int i;
  if (n<=0) return;
  for (i=0; i<n; i++) {
    index[i] = i;
    }
  iqs(item,index,0,n-1);
  }

static void
iqs(int *item,int *index,int left,int right)
{
  register int i,j;
  int x,y;
 
  i=left; j=right;
  x=item[index[(left+right)/2]];
 
  do {
    while(item[index[i]]>x && i<right) i++;
    while(x>item[index[j]] && j>left) j--;
 
    if (i<=j) {
      if (item[index[i]] != item[index[j]]) {
        y=index[i];
        index[i]=index[j];
        index[j]=y;
        }
      i++; j--;
      }
    } while(i<=j);
       
  if (left<j) iqs(item,index,left,j);
  if (i<right) iqs(item,index,i,right);
}

////////////////////////////////////////////////////////////////////
// Function findprocminmax finds the processor with the most/fewest
// basis functions and the corresponding number of basis functions
////////////////////////////////////////////////////////////////////
static void
findprocminmax(int *nbf, int nproc,
               int *procmin, int *procmax, int *minbf, int *maxbf)
{
  int i;

  *procmax = *procmin = 0;
  *maxbf = nbf[0];
  *minbf = nbf[0];

  for (i=1; i<nproc; i++) {
    if (nbf[i] > *maxbf) {
      *maxbf = nbf[i];
      *procmax = i;
      }
    if (nbf[i] < *minbf) {
      *minbf = nbf[i];
      *procmin = i;
      }
    }
}

/////////////////////////////////////////////////////////////////
// Function findshellmax finds the largest shell on a processor
/////////////////////////////////////////////////////////////////
static void
findshellmax(int *myshellsizes, int nRshell, int *shellmax, int *shellmaxindex)
{
  int i;

  *shellmax = myshellsizes[0];
  *shellmaxindex = 0;

  for (i=1; i<nRshell; i++) {
    if (myshellsizes[i] > *shellmax) {
      *shellmax = myshellsizes[i];
      *shellmaxindex = i;
      }
    }
}

//////////////////////////////////////////////////////////////
// Function expand_array expands the dimension of an array of
// doubles by 1; 
// NB: THE ARRAY MUST HAVE BEEN ALLOCATED WITH MALLOC
//////////////////////////////////////////////////////////////
static void 
expandintarray(int *&a, int olddim)
{
  int i;
  int *tmp;

  tmp = (int*) malloc((olddim+1)*sizeof(int));

  for (i=0; i<olddim; i++) {
    tmp[i] = a[i];
    }
  tmp[olddim] = 0;

  free(a);

  a = tmp;
}

////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
