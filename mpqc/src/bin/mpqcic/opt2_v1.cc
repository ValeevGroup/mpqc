
typedef int dmt_matrix;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <tmpl.h>

extern "C" {
#ifdef PARAGON
#include <nx.h>
#endif
}

#include <util/group/picl.h>
#include <tmpl.h>
#include <math/array/math_lib.h>
#include <math/dmt/libdmt.h>
#include <util/class/class.h>
#include <util/state/state.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <chemistry/qc/dmtsym/sym_dmt.h>
#include <chemistry/qc/dmtscf/scf_dmt.h>

extern "C" {
#include <chemistry/qc/dmtqc/libdmtqc.h>
#include <util/misc/libmisc.h>
}

extern "C" {
 int scf_erep_bound(int,int,int,int);
 int int_find_nfuncmax(centers_t*);
 int scf_init_bounds(centers_t*,double*);
 void scf_done_bounds();
 int gcollect(double*,int*,double*);
}

#include "opt2.h"
#include "bzerofast.h"

int
opt2_v1(centers_t *centers, scf_struct_t *scf_info, dmt_matrix Scf_Vec,
     double_vector_t *_evals, int nfzc, int nfzv, int mem_alloc, FILE* outfile)
{

  int_initialize_offsets2(centers,centers,centers,centers);

  int i, j, k;
  int s1, s2;
  int a, b;
  int isocc, asocc;   /* indices running over singly occupied orbitals */
  int nfuncmax = int_find_nfuncmax(centers);
  int nbasis = centers->nfunc;
  int nvir;
  int nocc=0;
  int ndocc=0,nsocc=0;
  int i_offset; 
  int npass, pass;
  int ni;             /* batch size */
  int np, nq, nr, ns; 
  int P, Q, R, S;
  int p, q, r, s;
  int bf1, bf2, bf3, bf4;
  int index;
  int col_index;
  int docc_index, socc_index, vir_index;
  int flags;
  int me;
  int nproc;
  int rest;
  int a_rest;
  int a_number;              /* number of a-values processed by each node */
  int a_offset;
  int *a_vector;           /* each node's # of iajb integrals for one i,j */
  int shell_index;
  int compute_index;
  int tmp_index;
  int dim_ij;
  int aoint_computed = 0;
  int nshell;
  double A, B, C, ni_top, max, ni_double; /* Variables used to compute ni  */
  double *evals = _evals->d;  /* scf eigenvalues (passed in)               */
  double *evals_open;    /* reordered scf eigenvalues                      */
  double *intbuf;        /* 2-electron AO integral buffer                  */
  double *trans_int1;    /* partially transformed integrals                */
  double *trans_int2;    /* partially transformed integrals                */
  double *trans_int3;    /* partially transformed integrals                */
  double *trans_int4_node;/* each node's subset of fully transf. integrals */
  double *trans_int4;    /* fully transformed integrals                    */
  double *mo_int_do_so_vir; /*mo integral (is|sa); i:d.o.,s:s.o.,a:vir     */
  double *mo_int_tmp;    /* scratch array used in gop1                     */
  double *socc_sum;      /* sum of 2-el integrals involving only s.o.'s    */
  double *iqrs, *iprs;
  double *iars_ptr, *iajs_ptr, *iajr_ptr;
  double iajr;
  double iars;
  double *iajb;
  double pqrs;
  double *c_qa;
  double *c_rb, *c_rj, *c_sj, *c_pi, *c_qi;
  double delta_ijab;
  double delta;
  double contrib1, contrib2;
  double ecorr_opt2=0,ecorr_opt1=0;
  double ecorr_zapt2;
  double ecorr_opt2_contrib=0, ecorr_zapt2_contrib=0;
  double escf;
  double eopt2,eopt1;
  double tol;     /* log2 of the erep tolerance (erep < 2^tol => discard) */
  loop_t *loop;


  me = mynode0();

  if (me == 0) {
    fprintf(outfile,"Just entered OPT2 program (opt2_v1)\n");
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


  /* compute number of a-values (a_number) processed by each node */

  a_number = nvir/nproc; 
  a_rest = nvir%nproc;
  if (me < a_rest) a_number++;

  if (a_number < nsocc) { 
    if (me == 0) fprintf(outfile,"not enough memory allocated\n");
    /* must have all socc's on node 0 for computation of socc_sum*/
    abort();
    }

  if (me < a_rest) a_offset = me*a_number; /* a_offset for each node */
  else a_offset = a_rest*(a_number + 1) + (me - a_rest)*a_number;

  /* fill in elements of a_vector for gcollect */

  a_vector = (int*) malloc(nproc*sizeof(int));
  if (!a_vector) {
    fprintf(stderr,"could not allocate storage for a_vector\n");
    abort();
    }
  for (i=0; i<nproc; i++) {
    a_vector[i] = nvir*(nvir/nproc)*sizeof(double);
    }
  for (i=0; i<a_rest; i++) {
    a_vector[i] += nvir*sizeof(double); /* first a_rest nodes hold an extra a */
    }


  /* compute batch size ni for opt2 loops                                 *
   * need to store the following arrays: trans_int1-4, trans_int4_node,   *
   * scf_vector, evals_open, socc_sum, mo_int_do_so_vir, mo_int_tmp and   *
   * a_vector;                                                            *
   * since a_number is not the same on all nodes, use node 0's a_number   *
   * (which is >= all other a_numbers) and broadcast ni afterwords        */

  if (me == 0) {
//    ni = (mem_alloc - sizeof(double)*(nvir*nvir + nbasis*(nvir+nocc)
//                                     + (nocc+nvir) + nsocc 
//                                     + 2*ndocc*nsocc*(nvir-nsocc) 
//                                     + nvir*a_number)
//                    - sizeof(int)*nproc)/
//          (sizeof(double)*(2*nfuncmax*nfuncmax*nbasis + nbasis*a_number*nocc));
    A = -0.5*sizeof(double)*nbasis*a_number;
    B = sizeof(double)*(2*nfuncmax*nfuncmax + (nocc+0.5)*a_number)*nbasis;
    C = sizeof(double)*(nvir*nvir + (nbasis+1)*(nocc+nvir) + nsocc 
                        + 2*ndocc*nsocc*(nvir-nsocc) + nvir*a_number)
        + sizeof(int)*nproc;
    ni_top = -B/(2*A);  /* ni value for which memory requirement is max. */
    max = A*ni_top*ni_top + B*ni_top + C;
    if (max <= mem_alloc) {
      ni = nocc;
      npass = 1;
      rest = 0;
      }
    else {
      ni_double = (-B + sqrt((double)(B*B - 4*A*(C-mem_alloc))))/(2*A);
      ni = (int) ni_double;
      if (ni > nocc) ni = nocc;
      rest = nocc%ni;
      npass = (nocc - rest)/ni + 1;
      if (rest == 0) npass--;
      }
    }
  bcast0(&ni,sizeof(int),mtype_get(),0);
  bcast0(&npass,sizeof(int),mtype_get(),0);
  bcast0(&rest,sizeof(int),mtype_get(),0);

  if (ni < nsocc) {
    if (me == 0) fprintf(outfile,"Not enough memory allocated\n");
    abort();
    }

  if (ni < 1) {     /* this applies only to a closed shell case */
    if (me == 0) fprintf(outfile,"Not enough memory allocated\n");
    abort();
    }

  if (me == 0) {
    fprintf(outfile,"computed batchsize: %i\n", ni);
    }

  nshell = centers->nshell;
  if (me == 0) {
    fprintf(outfile," npass  rest  nbasis  nshell  nfuncmax"
                    "  ndocc  nsocc  nvir  nfzc  nfzv\n");
    fprintf(outfile,"   %-4i   %-3i   %-5i   %-4i     %-3i"
                    "     %-3i     %-3i   %-3i    %-3i   %-3i\n",
            npass,rest,nbasis,nshell,nfuncmax,ndocc,nsocc,nvir,nfzc,nfzv);
    fprintf(outfile,"Using %i bytes of memory\n", mem_alloc);
    }


  /* rearrange scf eigenvalues as [socc docc socc unocc]    *
   * want socc first to get the socc's in the first batch   *
   * (need socc's to compute energy denominators - see      *
   * socc_sum comment below)                                */

  evals_open = (double*) malloc((nbasis+nsocc-nfzc-nfzv)*sizeof(double));
  if (!evals_open) {
    fprintf(stderr,"could not allocate storage for evals_open\n");
    abort();
    }
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

  /* allocate storage for integral arrays */

  dim_ij = nocc*ni - ni*(ni-1)/2;

  trans_int1 = (double*) malloc(nfuncmax*nfuncmax*nbasis*ni*sizeof(double));
  trans_int2 = (double*) malloc(nfuncmax*nfuncmax*nbasis*ni*sizeof(double));
  trans_int3 = (double*) malloc(nbasis*a_number*dim_ij*sizeof(double));
  trans_int4_node= (double*) malloc(nvir*a_number*sizeof(double));
  trans_int4 = (double*) malloc(nvir*nvir*sizeof(double));
  if (!(trans_int1 && trans_int2 && trans_int3 && trans_int4_node && trans_int4)){
    fprintf(stderr,"could not allocate storage for integral arrays\n");
    abort();
    }
  if (nsocc) socc_sum  = (double*) malloc(nsocc*sizeof(double));
  if (nsocc) mo_int_do_so_vir = 
                   (double*) malloc(ndocc*nsocc*(nvir-nsocc)*sizeof(double));
  if (nsocc) mo_int_tmp = 
                   (double*) malloc(ndocc*nsocc*(nvir-nsocc)*sizeof(double));

  if (nsocc) bzerofast(mo_int_do_so_vir,ndocc*nsocc*(nvir-nsocc));

/**************************************************************************
*    begin opt2 loops                                                     *
***************************************************************************/

  for (pass=0; pass<npass; pass++) {
    i_offset= pass*ni;  
    if ((pass == npass - 1) && (rest != 0)) ni = rest;
    bzerofast(trans_int3,nbasis*a_number*dim_ij);

    shell_index = 0;

    tim_enter("RS loop");
    for (R = 0; R < centers->nshell; R++) {
      nr = INT_SH_NFUNC((centers),R);

      for (S = 0; S <= R; S++) {
        ns = INT_SH_NFUNC((centers),S);
        tim_enter("bzerofast trans_int1");
        bzerofast(trans_int1,nfuncmax*nfuncmax*nbasis*ni);
        tim_exit("bzerofast trans_int1");

        tim_enter("PQ loop");
        for (P = 0; P < centers->nshell; P++) {
          np = INT_SH_NFUNC((centers),P);

          for (Q = 0; Q <= P; Q++) {
            shell_index++;
            if (shell_index%nproc != me) continue; 

            if (scf_erep_bound(P,Q,R,S) < tol) {
              continue;                           /* skip ereps less than tol */
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
                    if (R==S && bf4>bf3) {
                      index = ((bf1*nq + bf2)*nr + (bf3+1))*ns;
                      break; 
                      }

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

        tim_enter("gop1 int");
        gop1(trans_int1,nr*ns*nbasis*ni,trans_int2,'+',3);
        tim_exit("gop1 int");

        /* begin second quarter transformation */

        tim_enter("bzerofast trans_int2");
        bzerofast(trans_int2,nfuncmax*nfuncmax*nbasis*ni);
        tim_exit("bzerofast trans_int2");

        tim_enter("2. quart. tr.");

        for (bf3 = 0; bf3 < nr; bf3++) {

          for (bf4 = 0; bf4 < ns; bf4++) {
            if (R == S && bf4 > bf3) continue;

            for (q = 0; q < nbasis; q++) {
              c_qa = &scf_vector[q][nocc + a_offset];
              iqrs = &trans_int1[((bf4*nr + bf3)*nbasis + q)*ni];
              iars_ptr = &trans_int2[((bf4*nr + bf3)*a_number)*ni];

              for (a = 0; a < a_number; a++) {

                for (i=ni; i; i--) {
                  *iars_ptr++ += *c_qa * *iqrs++;
                  }

                iqrs -= ni;
                c_qa++;
                }
              }
            }
          }
        tim_exit("2. quart. tr.");

        /* begin third quarter transformation */
        tim_enter("3. quart. tr.");


        for (bf3 = 0; bf3<nr; bf3++) {
          r = centers->func_num[R] + bf3;

          for (bf4 = 0; bf4 <= (R == S ? bf3:(ns-1)); bf4++) {
            s = centers->func_num[S] + bf4;

            for (i=0; i<ni; i++) {
              tmp_index = i*(i+1)/2 + i*i_offset;

              for (a=0; a<a_number; a++) {
                iars = trans_int2[((bf4*nr + bf3)*a_number + a)*ni + i];
                if (r == s) iars *= 0.5;
                iajs_ptr = &trans_int3[tmp_index + dim_ij*(a + a_number*s)];
                iajr_ptr = &trans_int3[tmp_index + dim_ij*(a + a_number*r)];
                c_rj = scf_vector[r];
                c_sj = scf_vector[s];

                for (j=0; j<=i+i_offset; j++) {
                  *iajs_ptr++ += *c_rj++ * iars;
                  *iajr_ptr++ += *c_sj++ * iars;
                  }
                }
              }
            }     /* exit bf4 loop */
          }       /* exit bf3 loop */
        tim_exit("3. quart. tr.");
        }         /* exit S loop */
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

    tim_enter("4. quart. tr.");
    if (pass == 0 && me == 0) {
      if (nsocc) bzerofast(socc_sum,nsocc);
      for (isocc=0; isocc<nsocc; isocc++) {

        for (r=0; r<nbasis; r++) {

          for (asocc=0; asocc<nsocc; asocc++) {
            socc_sum[isocc] += scf_vector[r][nocc+asocc]*
                                trans_int3[isocc*(isocc+1)/2 + isocc*i_offset 
                                          + isocc + dim_ij*(asocc + a_number*r)];
            }
          }
        }
      }

    tim_enter("bcast0 socc_sum");
    if (nsocc) bcast0(socc_sum,nsocc*sizeof(double),mtype_get(),0);
    tim_exit("bcast0 socc_sum");

    tim_exit("4. quart. tr.");

    /* now we have all the sums of integrals involving s.o.'s (socc_sum);   *
     * begin fourth quarter transformation for all integrals (including     *
     * integrals with only s.o. indices); use restriction j <= (i_offset+i) *
     * to save flops                                                        */

    compute_index = 0;

    for (i=0; i<ni; i++) {

      for (j=0; j <= (i_offset+i); j++) {

       tim_enter("4. quart. tr.");

        bzerofast(trans_int4_node,nvir*a_number);

        for (r=0; r<nbasis; r++) {

          for (a=0; a<a_number; a++) {
            iajb = &trans_int4_node[a*nvir];
            c_rb = &scf_vector[r][nocc];
            iajr = trans_int3[i*(i+1)/2 + i*i_offset + j + dim_ij*(a+a_number*r)];

            for (b=0; b<nvir; b++) {
              *iajb++ += *c_rb++ * iajr;
              }
            }
          }

        tim_exit("4. quart. tr.");

        /* collect each node's part of fully transf. int. into trans_int4 */
        tim_enter("gcollect");
        gcollect(trans_int4_node,a_vector,trans_int4);
        tim_exit("gcollect");


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
                /* to compute the total energy contribution from an integral  *
                 * of the type (is1|s1a) (i=d.o., s1=s.o., a=unocc.), we need *
                 * the (is|sa) integrals for all s=s.o.; these integrals are  *
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
   * do this only is nsocc > 0 since gop1 will fail otherwise             */

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
                     npass*nshell*nshell*(nshell+1)*(nshell+1)/4);
    fprintf(outfile,"Number of shell quartets for which AO integrals\n"
                    "were computed: %i\n",aoint_computed);
    fprintf(outfile,"ROHF energy [au]:                  %13.8lf\n", escf);
    fprintf(outfile,"OPT1 energy [au]:                  %13.8lf\n", eopt1);
    fprintf(outfile,"OPT2 second order correction [au]: %13.8lf\n", ecorr_opt2);
    fprintf(outfile,"OPT2 energy [au]:                  %13.8lf\n", eopt2);
    fprintf(outfile,"ZAPT2 correlation energy [au]:     %13.8lf\n", ecorr_zapt2);
    fprintf(outfile,"ZAPT2 energy [au]:                 %13.8lf\n", 
                     escf + ecorr_zapt2);
    }

  int_done_erep();
  int_done_offsets2(centers,centers,centers,centers);
  scf_done_bounds();

  free(trans_int1);
  free(trans_int2);
  free(trans_int3);
  free(trans_int4_node);
  free(trans_int4);
  free(a_vector);
  if (nsocc) free(socc_sum);
  if (nsocc) free(mo_int_do_so_vir);
  if (nsocc) free(mo_int_tmp);
  free(evals_open);

  for (i=0; i<nbasis; i++) {
      delete[] scf_vector[i];
    }
  delete[] scf_vector;

  return(0);

  }
