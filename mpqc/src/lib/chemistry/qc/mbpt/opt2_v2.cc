//
// opt2_v2.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Ida Nielsen <ibniels@ca.sandia.gov>
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

typedef int dmt_matrix;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <util/group/picl.h>
#include <util/class/class.h>
#include <util/state/state.h>
#include <math/array/math_lib.h>
#include <math/dmt/libdmt.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <chemistry/qc/dmtsym/sym_dmt.h>
#include <chemistry/qc/dmtscf/scf_dmt.h>

extern "C" {
#include <chemistry/qc/dmtqc/libdmtqc.h>
#include <util/misc/libmisc.h>
}

extern "C" {
 int int_find_nfuncmax(centers_t*);
 int scf_init_bounds(centers_t*,double*);
 int scf_erep_bound(int,int,int,int);
 void scf_done_bounds();
}

static void iqs(int *item,int *index,int left,int right);
static void iquicksort(int *item,int *index,int n);

#include <chemistry/qc/mbpt/opt2.h>
#include <chemistry/qc/mbpt/bzerofast.h>

int
mbpt_opt2_v2(centers_t *centers, scf_struct_t *scf_info, dmt_matrix Scf_Vec,
             double_vector_t *_evals, int nfzc, int nfzv, int mem_alloc,
             FILE* outfile)
{

  int_initialize_offsets2(centers,centers,centers,centers);

  int i, j, k;
  int s1, s2;
  int a, b;
  int isocc, asocc;   /* indices running over singly occupied orbitals */
  int nfuncmax = int_find_nfuncmax(centers);
  int nbasis = centers->nfunc;
  int nvir;
  int nshell;
  int nocc=0,ndocc=0,nsocc=0;
  int i_offset; 
  int npass, pass;
  long ni;
  int np, nq, nr, ns; 
  int P, Q, R, S;
  int p, q, r, s;
  int bf1, bf2, bf3, bf4;
  int index;
  int compute_index;
  int col_index;
  int tmp_index;
  int dim_ij;
  int docc_index, socc_index, vir_index;
  int flags;
  int me;
  int nproc;
  int rest;
  int r_offset;
  int sum;
  int min;
  int iproc;
  int nRshell;
  int imyshell;
  int *myshells;       /* the R indices processed by node me */
  int *shellsize;      /* size of each shell                 */
  int *sorted_shells;  /* sorted shell indices: large shells->small shells */
  int *nbf;            /* number of basis functions processed by each node */
  int *proc;           /* element k: processor which will process shell k  */
  int aoint_computed = 0; 
  double A, B, C, ni_top, max, ni_double; /* variables used to compute ni  */
  double *evals = _evals->d;  /* scf eigenvalues (passed in)               */
  double *evals_open;    /* reordered scf eigenvalues                      */
  double *intbuf;        /* 2-electron AO integral buffer                  */
  double *trans_int1;    /* partially transformed integrals                */
  double *trans_int2;    /* partially transformed integrals                */
  double *trans_int3;    /* partially transformed integrals                */
  double *trans_int4;    /* fully transformed integrals                    */
  double *trans_int4_tmp; /* scratch array                                 */ 
  double *mo_int_do_so_vir; /*mo integral (is|sa); i:d.o.,s:s.o.,a:vir     */
  double *mo_int_tmp;    /* scratch array used in gop1                     */
  double *socc_sum;      /* sum of 2-el integrals involving only s.o.'s    */
  double *socc_sum_tmp;  /* scratch array                                  */ 
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
  double eopt2,eopt1;
  double tol;          /* log2 of the erep tolerance (erep < 2^tol => discard) */
  loop_t *loop;

  me = mynode0();

  if (me == 0) {
    fprintf(outfile,"Just entered OPT2 program (opt2_v2)\n");
    }

  flags = INT_EREP|INT_NOSTRB|INT_NOSTR1|INT_NOSTR2;
  intbuf = int_initialize_erep(flags,0,centers,centers,centers,centers);

  scf_init_bounds(centers,intbuf);

  tol = (int) (-10.0/log10(2.0));  /* discard ereps smaller than 10^-10 */

  nproc = numnodes0();
  if (me == 0) fprintf(outfile,"nproc = %i\n", nproc);

  ndocc = scf_info->nclosed;
  nsocc = scf_info->nopen;

  /* do a few preliminary tests to make sure the desired calculation *
   * can be done (and appears to be meaningful!)                     */

  if (ndocc == 0 && nsocc == 0) {
    if (me == 0) {
      fprintf(outfile,"There are no occupied orbitals; program exiting\n");
      }
    abort();
    }

  if (nfzc > ndocc) {
    if (me == 0) {
      fprintf(outfile,"The number of frozen core orbitals exceeds the number\n"
                      "of doubly occupied orbitals; program exiting\n");
      }
    abort();
    }

  if (nfzv > nbasis - ndocc - nsocc) {
    if (me == 0) {
      fprintf(outfile,"The number of frozen virtual orbitals exceeds the number\n"
                      "of unoccupied orbitals; program exiting\n");
      }
    abort();
    }

  ndocc = ndocc - nfzc;
  /* nvir = # of unocc. orb. + # of s.o. orb. - # of frozen virt. orb. */
  nvir  = nbasis - ndocc - nfzc - nfzv; 
  /* nocc = # of d.o. orb. + # of s.o. orb - # of frozen d.o. orb. */
  nocc  = ndocc + nsocc;
  nshell = centers->nshell;

  /* allocate storage for some arrays used for keeping track of which R   *
   * indices are processed by each node                                   */
  shellsize = (int*) malloc(nshell*sizeof(int));
  sorted_shells = (int*) malloc(nshell*sizeof(int));
  nbf = (int*) malloc(nproc*sizeof(int));
  proc = (int*) malloc(nshell*sizeof(int));


  /******************************************************
  * Begin distributing R shells between nodes so each   *
  * node gets ca. the same number of r basis functions  *
  ******************************************************/

  /* compute size of each shell */
  for (i=0; i<nshell; i++) {
    shellsize[i] = INT_SH_NFUNC((centers),i);
    }

  /* do an index sort (large -> small) of shellsize to form sorted_shells */
  iquicksort(shellsize,sorted_shells,nshell);

  /* initialize nbf */
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
    fprintf(outfile,"Distribution of basis functions between nodes:\n");
    for (i=0; i<nproc; i++) {
      fprintf(outfile," %4i",nbf[i]);
      if ((i+1)%12 == 0) fprintf(outfile,"\n");
      }
    fprintf(outfile,"\n");
    }

  /* determine which shells are to be processed by node me */
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

  /************************************************
  * End of distribution of R shells between nodes *
  ************************************************/


  /* compute batch size ni for opt2 loops;                                *
   * need to store the following arrays of type double : trans_int1-4,    *
   * trans_int4_tmp, scf_vector, evals_open, socc_sum, socc_sum_tmp,      *
   * mo_int_do_so_vir, mo_int_tmp,                                        *
   * and the following arrays of type int: myshells, shellsize,           *
   * sorted_shells, nbf, and proc                                         */
//  ni = (mem_alloc 
//        - sizeof(double)*(2*nvir*nvir + nbasis*(nvir+nocc) + (nocc+nvir) 
//                         + 2*nsocc + 2*ndocc*nsocc*(nvir-nsocc))
//        - sizeof(int)*(3*nshell + nproc + nRshell)) /
//        (sizeof(double)*(nfuncmax*nfuncmax*nbasis + nvir + nbf[me]*nvir*nocc));
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
      }

  /* set ni equal to the smallest batch size for any node */
  gmin0(&ni,1,2,mtype_get(),0);
  bcast0(&ni,sizeof(int),mtype_get(),0);
  /* gilow(&ni,1,&tmpint);  paragon specific */

  if (ni < nsocc) {
    if (me == 0) fprintf(outfile,"Not enough memory allocated\n");
    abort();
    }

  if (ni < 1) {     /* this applies only to a closed shell case */
    if (me == 0) fprintf(outfile,"Not enough memory allocated\n");
    abort();
    }

  if (me == 0) fprintf(outfile,"Computed batchsize: %i\n",ni);

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
    fprintf(outfile," npass  rest  nbasis  nshell  nfuncmax"
                    "  ndocc  nsocc  nvir  nfzc  nfzv\n");
    fprintf(outfile,"  %-4i   %-3i   %-5i    %-4i     %-3i"
                    "     %-3i    %-3i    %-3i    %-3i   %-3i\n",
            npass,rest,nbasis,nshell,nfuncmax,ndocc,nsocc,nvir,nfzc,nfzv);
    fprintf(outfile,"Using %i bytes of memory\n",mem_alloc);
    }


  /* rearrange scf eigenvalues as [socc docc socc unocc]    *
   * want socc first to get the socc's in the first batch   *
   * (need socc's to compute energy denominators - see      *
   * socc_sum comment below)                                */

  evals_open = (double*) malloc((nbasis+nsocc-nfzc-nfzv)*sizeof(double));
  for (i=0; i < nsocc; i++) evals_open[i] = evals[i+nfzc+ndocc];
  for (i=nsocc; i < nocc; i++) evals_open[i] = evals[i-nsocc+nfzc];
  for (i=nocc; i < nvir+nocc; i++) evals_open[i] = evals[i+nfzc-nsocc];


  /* the scf vector is distributed between the nodes, but for OPT2  *
   * each node needs its own copy of the vector;                    *
   * therefore, put a copy of the scf vector on each node;          * 
   * while doing this, duplicate columns corresponding to singly    *
   * occupied orbitals and order columns as [socc docc socc unocc]  */

  double** scf_vector = new double*[nbasis];
  for (i=0; i<nbasis; i++) {
      scf_vector[i] = new double[nocc+nvir];
    }

  loop = dmt_ngl_create("%mr",Scf_Vec);
  while(dmt_ngl_next(loop)) {
    int iind,isize,jsize;
    double *col;

    dmt_ngl_create_inner(loop,0);
    while(dmt_ngl_next_inner_m(loop,&iind,&isize,&k,&jsize,&col)) {
      if (k >= nfzc && k < ndocc+nfzc) {
        for (i=0; i < nbasis; i++) scf_vector[i][k-nfzc+nsocc] = col[i];
        }
      if (k >= ndocc+nfzc && k < ndocc+nfzc+nsocc) {
        for (i=0; i < nbasis; i++) {
          scf_vector[i][k-ndocc-nfzc] = col[i];
          scf_vector[i][k-nfzc+nsocc] = col[i];
          }
        }
      if (k >= ndocc+nsocc+nfzc && k < nbasis-nfzv) {
        for (i=0; i < nbasis; i++) scf_vector[i][k-nfzc+nsocc] = col[i];
        }
      }
    }

  dmt_ngl_kill(loop);



  /* allocate storage for various arrays */

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


/**************************************************************************
 *   begin opt2 loops                                                     *
 **************************************************************************/


  for (pass=0; pass<npass; pass++) {
    i_offset = pass*ni;  
    if ((pass == npass - 1) && (rest != 0)) ni = rest;

    r_offset = 0;
    bzerofast(trans_int3,nbf[me]*nvir*dim_ij);

    tim_enter("RS loop");

    for (imyshell=0; imyshell<nRshell; imyshell++) {

      R = myshells[imyshell];
      nr = INT_SH_NFUNC((centers),R);

      for (S = 0; S < nshell; S++) {
        ns = INT_SH_NFUNC((centers),S);
        tim_enter("bzerofast trans_int1");
        bzerofast(trans_int1,nfuncmax*nfuncmax*nbasis*ni);
        tim_exit("bzerofast trans_int1");

        tim_enter("PQ loop");

        for (P = 0; P < nshell; P++) {
          np = INT_SH_NFUNC((centers),P);

          for (Q = 0; Q <= P; Q++) {
            if (scf_erep_bound(P,Q,R,S) < tol) {
              continue;                          /* skip ereps less than tol */
              }

            aoint_computed++;

            nq = INT_SH_NFUNC((centers),Q);

            tim_enter("erep");
            int_erep(INT_EREP|INT_NOBCHK|INT_NOPERM|INT_REDUND,&P,&Q,&R,&S);
            tim_exit("erep");

            tim_enter("1. quart. tr.");

            index = 0;

            for (bf1 = 0; bf1 < np; bf1++) {
              p = centers->func_num[P] + bf1;
 
              for (bf2 = 0; bf2 < nq; bf2++) {
                q = centers->func_num[Q] + bf2;
                if (q > p) {
                  /* if q > p: want to skip the loops over bf3-4  */
                  /* and larger bf2 values, so increment bf1 by 1 */
                  /* ("break") and adjust the value of index      */
                  index = (bf1 + 1) * nq * nr * ns;
                  break;
                  }

                for (bf3 = 0; bf3 < nr; bf3++) {

                  for (bf4 = 0; bf4 < ns; bf4++,index++) {

                    if (INT_NONZERO(intbuf[index])) {
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
                    }   /* exit bf4 loop */
                  }     /* exit bf3 loop */
                }       /* exit bf2 loop */
              }         /* exit bf1 loop */
            tim_exit("1. quart. tr.");
            }           /* exit Q loop */
          }             /* exit P loop */
        tim_exit("PQ loop");

        /* begin second and third quarter transformations */

        for (bf3 = 0; bf3 < nr; bf3++) {
          r = r_offset + bf3;

          for (bf4 = 0; bf4 < ns; bf4++) {
            s = centers->func_num[S] + bf4;

            tim_enter("bzerofast trans_int2");
            bzerofast(trans_int2,nvir*ni);
            tim_exit("bzerofast trans_int2");

            tim_enter("2. quart. tr.");

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
              }             /* exit q loop */
            tim_exit("2. quart. tr.");

            /* begin third quarter transformation */

            tim_enter("3. quart. tr.");

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
              }   /* exit i loop */
            tim_exit("3. quart. tr.");

            } /* exit bf4 loop */
          }   /* exit bf3 loop */

        }         /* exit S loop */
      r_offset += nr;
      }           /* exit R loop */
    tim_exit("RS loop");


    /* begin fourth quarter transformation;                                *
     * first tansform integrals with only s.o. indices;                    *
     * these integrals are needed to compute the denominators              *
     * in the various terms contributing to the correlation energy         *
     * and must all be computed in the first pass;                         *
     * the integrals are summed into the array socc_sum:                   *
     * socc_sum[isocc] = sum over asocc of (isocc asocc|asocc isocc)       *
     * (isocc, asocc = s.o. and the sum over asocc runs over all s.o.'s)   *
     * the individual integrals are not saved here, only the sums are kept */


    if (pass == 0) {
      tim_enter("4. quart. tr.");
      if (nsocc) bzerofast(socc_sum,nsocc);
      for (isocc=0; isocc<nsocc; isocc++) {

        for (index=0; index<nbf[me]; index++) {
          i = 0;
          sum = INT_SH_NFUNC((centers),myshells[i]);
          while (sum <= index) {
            i++;
            sum += INT_SH_NFUNC((centers),myshells[i]);
            }
          sum -= INT_SH_NFUNC((centers),myshells[i]);
          r = centers->func_num[myshells[i]] + index - sum;

          for (asocc=0; asocc<nsocc; asocc++) {
            socc_sum[isocc] += scf_vector[r][nocc+asocc]*
                           trans_int3[isocc*(isocc+1)/2 + isocc*i_offset
                                     + isocc + dim_ij*(asocc + nvir*index)];
            }
          }
        }       /* exit isocc loop */

      tim_exit("4. quart. tr.");

      /* sum socc_sum contributions from each node (only if nsocc > 0 *
       * since gop1 will fail if nsocc = 0)                           */
      if (nsocc > 0) {
        tim_enter("gop1 socc_sum");
        gop1(socc_sum,nsocc,socc_sum_tmp,'+',3);
        tim_exit("gop1 socc_sum");
        }

      } 

    /* now we have all the sums of integrals involving s.o.'s (socc_sum);   *
     * begin fourth quarter transformation for all integrals (including     *
     * integrals with only s.o. indices); use restriction j <= (i_offset+i) *
     * to save flops                                                        */

    compute_index = 0;

    for (i=0; i<ni; i++) {

      for (j=0; j <= (i_offset+i); j++) {

        tim_enter("4. quart. tr.");

        bzerofast(trans_int4,nvir*nvir);

        for (index=0; index<nbf[me]; index++) {
          k = 0;
          sum = INT_SH_NFUNC((centers),myshells[k]);
          while (sum <= index) {
            k++;
            sum += INT_SH_NFUNC((centers),myshells[k]);
            }
          sum -= INT_SH_NFUNC((centers),myshells[k]);
          r = centers->func_num[myshells[k]] + index - sum;

          for (a=0; a<nvir; a++) {
            iajb = &trans_int4[a*nvir];
            iajr = trans_int3[i*(i+1)/2 + i*i_offset + j + dim_ij*(a+nvir*index)];
            c_rb = &scf_vector[r][nocc];

            for (b=0; b<nvir; b++) {
              *iajb++ += *c_rb++ * iajr;
              }
            }
          }

        tim_exit("4. quart. tr.");

        tim_enter("gop1 trans_int4");
        gop1(trans_int4,nvir*nvir,trans_int4_tmp,'+',3);
        tim_exit("gop1 trans_int4");

        /* we now have the fully transformed integrals (ia|jb)      *
         * for one i, one j (j <= i_offset+i), and all a and b;     *
         * compute contribution to the OPT1 and OPT2 correlation    *
         * energies; use restriction b <= a to save flops           */

        tim_enter("compute ecorr");

        for (a=0; a<nvir; a++) {
          for (b=0; b<=a; b++) {
            compute_index++;
            if (compute_index%nproc != me) continue;

            docc_index = ((i_offset+i) >= nsocc && (i_offset+i) < nocc) 
                        + (j >= nsocc && j < nocc);
            socc_index = ((i_offset+i)<nsocc)+(j<nsocc)+(a<nsocc)+(b<nsocc);
            vir_index = (a >= nsocc) + (b >= nsocc);

            if (socc_index >= 3) continue; /* skip to next b value */
 
            delta_ijab = evals_open[i_offset+i] + evals_open[j] 
                       - evals_open[nocc+a] - evals_open[nocc+b];
            
            /* determine integral type and compute energy contribution */
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
                /* to compute the energy contribution from an integral of the *
                 * type (is1|s1a) (i=d.o., s1=s.o., a=unocc.), we need the    *
                 * (is|sa) integrals for all s=s.o.; these integrals are      *
                 * therefore stored here in the array mo_int_do_so_vir, and   *
                 * the energy contribution is computed after exiting the loop *
                 * over i-batches (pass)                                      */
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
            }   /* exit b loop */
          }     /* exit a loop */
        tim_exit("compute ecorr");
        }       /* exit j loop */
      }         /* exit i loop */
    }           /* exit loop over i-batches (pass) */



  /* compute contribution from excitations of the type is1 -> s1a where   *
   * i=d.o., s1=s.o. and a=unocc; single excitations of the type i -> a,  *
   * where i and a have the same spin, contribute to this term;           *
   * (Brillouin's theorem not satisfied for ROHF wave functions);         *
   * do this only if nsocc > 0 since gop1 will fail otherwise             */

  tim_enter("compute ecorr");

  if (nsocc > 0) {
    tim_enter("gop1 mo_int_do_so_vir");
    gop1(mo_int_do_so_vir,ndocc*nsocc*(nvir-nsocc),mo_int_tmp,'+',3);
    tim_exit("gop1 mo_int_do_so_vir");
    }

  /* add extra contribution for triplet and higher spin multiplicities *
   * contribution = sum over s1 and s2<s1 of (is1|s1a)*(is2|s2a)/delta */

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
        }     /* exit a loop */
      }       /* exit i loop */
    }

  tim_exit("compute ecorr");

  ecorr_zapt2 = ecorr_opt2 + ecorr_zapt2_contrib;
  ecorr_opt2 += ecorr_opt2_contrib;
  gsum0(&ecorr_opt1,1,5,mtype_get(),0);
  gsum0(&ecorr_opt2,1,5,mtype_get(),0);
  gsum0(&ecorr_zapt2,1,5,mtype_get(),0);
  gsum0(&aoint_computed,1,2,mtype_get(),0);

  if (me == 0) {
    escf = scf_info->e_elec + scf_info->nuc_rep;
    eopt2 = escf + ecorr_opt2;
    eopt1 = escf + ecorr_opt1;

    /* print out various energies etc.*/

    fprintf(outfile,"Number of shell quartets for which AO integrals would\n"
                    "have been computed without bounds checking: %i\n",
                     npass*nshell*nshell*(nshell+1)*(nshell+1)/2);
    fprintf(outfile,"Number of shell quartets for which AO integrals\n"
                    "were computed: %i\n",aoint_computed);

    fprintf(outfile,"ROHF energy [au]:                  %13.8lf\n", escf);
    fprintf(outfile,"OPT1 energy [au]:                  %13.8lf\n", eopt1);
    fprintf(outfile,"OPT2 second order correction [au]: %13.8lf\n", ecorr_opt2);
    fprintf(outfile,"OPT2 energy [au]:                  %13.8lf\n", eopt2);
    fprintf(outfile,"ZAPT2 correlation energy [au]:     %13.8lf\n", ecorr_zapt2);
    fprintf(outfile,"ZAPT2 energy [au]:                 %13.8lf\n", 
               escf + ecorr_zapt2);
    fflush(outfile);
    }

  int_done_erep();
  int_done_offsets2(centers,centers,centers,centers);
  scf_done_bounds();

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
  free(sorted_shells);
  free(nbf);
  free(proc);

  for (i=0; i<nbasis; i++) {
      delete[] scf_vector[i];
    }
  delete[] scf_vector;

  return(0);

  }


/* Do a quick sort (larger -> smaller) of the integer data in item *
 * by the integer indices in index;                                *
 * data in item remain unchanged                                  */

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


////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ-CONDENSED")
// End:
