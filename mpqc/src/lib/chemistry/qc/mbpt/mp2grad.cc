typedef int dmt_matrix;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <tmpl.h>

#include <util/group/picl.h>
#include <util/group/memory.h>
#include <util/group/message.h>
#include <util/class/class.h>
#include <util/state/state.h>
#include <math/array/math_lib.h>
#include <math/dmt/libdmt.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <chemistry/qc/dmtsym/sym_dmt.h>
#include <chemistry/qc/dmtscf/scf_dmt.h>
#include <math/scmat/matrix.h>
#include <math/scmat/local.h>

extern "C" {
#include <chemistry/qc/dmtqc/libdmtqc.h>
#include <util/misc/libmisc.h>
}

#include <chemistry/qc/mbpt/gmat.h>
#include <chemistry/qc/mbpt/ffo.h>
#include <chemistry/qc/mbpt/cphf.h>
#include <chemistry/qc/mbpt/mp2grad.h>

static int
mp2grad(centers_t *centers, scf_struct_t *scf_info, dmt_matrix Scf_Vec,
        double_vector_t *_evals, int nfzc, int nfzv, int mem,FILE* outfile,
        RefMessageGrp msg, double_matrix_t *);

static void
s2pdm_contrib(double *intbuf, centers_t *centers, double *PHF, double *P2AO,
              double_matrix_t *ginter, int nproc, int me);

static void
hcore_grad(centers_t *centers, double *PMP2, double_matrix_t *ginter,
           int nproc, int me);

static void
overlap_grad(centers_t *centers, double *WMP2, double_matrix_t *ginter,
               int nproc, int me);

static void
sum_gradients(double_matrix_t *f, int nproc);

static int
compute_batchsize(int mem_alloc, int mem_static, int nproc, FILE* outfile);

extern "C" {
 int malloc_chain_check(int);
 int int_find_nfuncmax(centers_t*);
 int scf_erep_bound(int,int,int,int);
 int int_erep_4bound(int,int,int,int);
 int bzerofast(double *, int);
 void scf_done_bounds();
}

static int messagebuf[2];
static int nocc, nvir, nbasis;
static int send_offset;
static int nfuncmax;
static int nij;        // number of i,j pairs on a node (for e.g., mo_int)
static double *mo_int; // MO integrals of type (ov|ov)
                       // (and these integrals divided by
                       // orbital energy denominators)
static double *integral_iqjs;
static double *iqjs_buf;
static double *gamma_iqjs;

#define PRINT1Q 0

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

static int
mp2grad(centers_t *centers, scf_struct_t *scf_info, dmt_matrix Scf_Vec,
     double_vector_t *_evals, int nfzc, int nfzv, int mem_alloc, FILE* outfile,
        RefMessageGrp msg, double_matrix_t *gradientt)
{

  // New version of MP2 gradient program which uses the full
  // permutational symmetry of the two-electron integral derivatives

  int_initialize_offsets2(centers,centers,centers,centers);

  int i, j, k, l, m;
  int y;
  int isize, jsize;
  int s1, s2;
  int a, b, c;
  int nshell;
  int offset, offset2;
  int ik_offset;
  int i_offset, j_offset; 
  int npass, pass;
  long tmpint;
  int np, nq, nr, ns; 
  int P, Q, R, S;
  int p, q, r, s;
  int bf1, bf2, bf3, bf4;
  int index;
  int compute_index;
  int flags;
  int me;
  int nproc;
  int rest;
  int node;
  int p_offset, q_offset, r_offset, s_offset;
  int sum;
  int nfunc;
                        //  not split 
  int aoint_computed = 0; 
  int aointder_computed = 0; 
  int derset, xyz;
  int natom = centers->n;     // the number of atoms
  int int_index;
  int mem_static, mem_dyn;    // static and dynamic memory in bytes
  int mem_used;
  int qp, sr;
  int factor_pqrs;
  int factor_pq, factor_rs;
  int ij_proc;          // the processor which has ij pair
  int ij_index;         // of the ij pairs on a proc, this ij pair is number ij_index
                        // (i.e., ij_index < nij)
  int ik_proc;          // the processor which has ik pair
  int ik_index;
  int sendbuf[2];
  int ij_offset;
  int jloop, kloop;
  int ntri;

  long ni;

  double A, B, C;             // variables used to compute ni
  double  maxdyn, maxdyntmp;  // variables used to compute ni
  double *evals = _evals->d;  // scf eigenvalues (passed in)
  double *intbuf;             // 2-electron AO integral buffer
  double *iqrs, *iprs;
  double *iars_ptr;
  double iars;
  double iajr;
  double *iajr_ptr;
  double *iajb_ptr, *ibja_ptr, *iakb_ptr, *ibka_ptr;
  double *iajc_ptr, *ibjc_ptr, *icjb_ptr, *icja_ptr;
  double *ijkb_ptr, *ijkr_ptr, *ibkr_ptr, *ibkj_ptr;
  double *ibac_ptr, *icab_ptr, *ibar_ptr;
  double pqrs;
  double *c_sa, c_qj, c_pj, c_rj;
  double *c_qk, *c_rb, *c_rk, *c_sk, *c_pi, *c_qi, *c_sj;
  double *c_qa, *c_rc, *c_sb, *c_pa, *c_pq, *c_sy;
  double delta_ijab, delta_ikab, delta_ijbc, delta_ijac;
  double ecorr_mp2 = 0.0;
  double escf;
  double emp2;
  double tol;                 // log2 of the erep tolerance
                              // (erep < 2^tol => discard)
  double *Wkj,*Wab,*Waj;      // occ-occ, vir-vir and vir-occ parts of 
                              // second order correction to MP2
                              // energy weighted density matrix
  double *Pkj,*Pab;           // occ-occ and vir-vir parts of second order
                              // correction to MP2 density matrix
  double *Laj;                // MP2 Lagrangian
  double *Lpi;                // contrib to MP2 Lagrangian partially in AO basis
  double *pkj_ptr, pjk, *pab_ptr, *pbc_ptr, pcb;
  double *wkj_ptr, *wjk_ptr, *wab_ptr, *wba_ptr, *waj_ptr;
  double *laj_ptr, *lpi_ptr, *lqi_ptr;
  double *gamma_iajs, *gamma_iajs_tmp, *gamma_iqrs; 
                              // partially back-transformed non-sep 2PDM's
  double *gamma_iqjs_tmp;
  double *gamma_pqrs;
  double *gamma_ipqs, *gamma_iqps;
  double *gamma_iajs_ptr;
  double *gamma_ipqr_ptr, *gamma_iqpr_ptr;
  double *gamma_ipjs_ptr, *gamma_iqjs_ptr, *gamma_irjq_ptr;
  double *gamma_iqrs_ptr, *gamma_iprs_ptr;
  double *gamma_iqsr_ptr;
  double *gamma_pqrs_ptr;
  double *gammabuf;           // buffer used for sending elements of gamma_iqjs
  double *mo_intbuf;          // buffer used for sending mo integrals
  double *grad_ptr1, *grad_ptr2;
  double tmpval, tmpval1, tmpval2;
  double *Dmat;
  double *dmat_ptr;
  double *P2AO, *W2AO;
  double *p2ao_ptr, *w2ao_ptr;
  double *p2ao_pq, *p2ao_ps;
  double *PHF, *WHF;
  double *phf_ptr, *whf_ptr;
  double *phf_pq, *phf_rs, *phf_ps, *phf_qr;
  double *PMP2, *WMP2;
  double *pmp2_ptr, *wmp2_ptr;

  double *integral_iqrs; // quarter transformed two-el integrals
  double *integral_iajs; // three-quarter transformed two-el integrals
  double *integral_ikjs; // three-quarter transformed two-el integrals
  double *integral_ixjs;  // all three-quarter transformed two-el integrals
  double *integral_iajy; // mo integrals (y = any MO)
  double *integral_ikja; // mo integrals
  double *iqjs_contrib;  // local contributions to integral_iqjs
  double *iqjr_contrib;  // local contributions to integral_iqjr
  double *integral_iqjs_ptr, *integral_iqjr_ptr;
  double *iajy_ptr;
  double *iajk_ptr, *ikja_ptr;
  double *iajs_ptr, *ikjs_ptr;
  double *iqrs_ptr, *iprs_ptr;
  double *iqjs_ptr, *iqjr_ptr;
  double *pqrs_ptr;

  double_vector_t repder; // Intermediates for nuc. rep. contrib. to grad.
  double_matrix_t *gradient;  // The MP2 gradient
  double_matrix_t *ginter;    // Intermediates for the MP2 gradient

  nbasis = centers->nfunc;
  nfuncmax = int_find_nfuncmax(centers);

  der_centers_t der_centers;

  loop_t *loop;

  nshell = centers->nshell;

  me = mynode0();

  if (me == 0) gradient = new double_matrix_t;
  ginter   = new double_matrix_t;

  if (me == 0) {
    fprintf(outfile,"Entered MP2 program (mp2grad)\n");
    fflush(outfile);
    }
  
  nproc = numnodes0();
  if (me == 0) fprintf(outfile,"nproc = %i\n", nproc);

  // Set flags for erep and derivative erep integral evaluation
  flags = INT_EREP|INT_REDUND|INT_NOPERM|INT_NOSTRB;

  // Initialize erep integral buffer
  intbuf = int_initialize_erep(flags,1,centers,centers,centers,centers);

  // Initialize integral (derivative) bounds (in parallel)
  int_init_bounds_1der_nocomp();
  index = 0;
  for (R=0; R<nshell; R++) {
    for (S=0; S<=R; S++) {
      if (index++ % nproc == me) int_bounds_1der_comp(R,S);
      }
    }
  gop0_sc(&int_Q,1,'M',mtype_get());
  gop0_sc(&int_R,1,'M',mtype_get());
  ntri = nshell*(nshell+1)/2;
  gop0_sc(int_Qvec,ntri,'+',mtype_get());
  gop0_sc(int_Rvec,ntri,'+',mtype_get());


//scf_init_bounds(centers,intbuf);

  tol = (int) (-10.0/log10(2.0));  // discard ereps smaller than 10^-10

  nocc = scf_info->nclosed;

  // Do a few preliminary tests to make sure the desired calculation
  // can be done (and appears to be meaningful!)

  if (nocc == 0) {
    if (me == 0) {
      fprintf(outfile,"There are no occupied orbitals; program exiting\n");
      }
    abort();
    }


  if (nfzc != 0) {
    if (me == 0) {
      fprintf(outfile,"The number of frozen core orbitals is nonzero;\n"
                      "no orbitals can be frozen currently; program exits\n");
      }
    abort();
    }

  if (nfzv != 0) {
    if (me == 0) {
      fprintf(outfile,"The number of frozen virtual orbitals is nonzero;\n"
                      "no orbitals can be frozen currently; program exits\n");
      }
    abort();
    }

  nvir  = nbasis - nocc - nfzc - nfzv; 

  if (nvir == 0) {
    if (me == 0) {
      fprintf(outfile,"There are no virtual orbitals; program exiting\n");
      }
    abort();
    }

  nocc = nocc - nfzc;

  ////////////////////////////////////////////////////////
  // Compute batch size ni for mp2 loops;
  //
  // The following arrays are kept throughout (all of type double):
  //   scf_vector, gradient, ginter, Pkj, Pab, Wkj, Wab, Waj, Laj, iqjs_buf
  // and memory allocated for these arrays is called mem_static
  //
  ////////////////////////////////////////////////////////
  if (me == 0) {
    mem_static = sizeof(double)*(nbasis*nbasis + 6*natom + nocc*(nocc+1)/2
                               + nvir*(nvir+1)/2 + nocc*nocc + nvir*nvir
                               + 2*nocc*nvir + 2+nbasis*nfuncmax);
    ni = compute_batchsize(mem_alloc, mem_static, nproc, outfile); 
    }

  // Send value of ni to other nodes
  bcast0(&ni,sizeof(int),mtype_get(),0);

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
    fprintf(outfile," npass  rest  nbasis  nshell  nfuncmax  nocc  nvir\n");
    fprintf(outfile,"  %-4i   %-3i   %-5i    %-4i     %-3i     %-3i    %-3i\n",
            npass,rest,nbasis,nshell,nfuncmax,nocc,nvir);
    }

  ////////////////////////////////////////////////
  // The scf vector is distributed between nodes;
  // put a copy of the scf vector on each node;
  ////////////////////////////////////////////////

  double** scf_vector = new double*[nbasis];
  for (i=0; i<nbasis; i++) {
      scf_vector[i] = new double[nocc+nvir];
    }

  loop = dmt_ngl_create("%mr",Scf_Vec);
  while(dmt_ngl_next(loop)) {
    int iind;
    double *col;

    dmt_ngl_create_inner(loop,0);
    while(dmt_ngl_next_inner_m(loop,&iind,&isize,&k,&jsize,&col)) {
      if (k >= nfzc && k < nocc+nfzc) {
        for (i=0; i < nbasis; i++) scf_vector[i][k-nfzc] = col[i];
        }
      if (k >= nocc+nfzc && k < nbasis-nfzv) {
        for (i=0; i < nbasis; i++) scf_vector[i][k-nfzc] = col[i];
        }
      }
    }

  dmt_ngl_kill(loop);

  //////////////////////////////////////////////////////////////
  // Allocate storage for various arrays
  // (Pkj and Pab are symmetric, so store only lower triangle)
  //////////////////////////////////////////////////////////////

  Pkj            = (double*) malloc((nocc*(nocc+1)/2)*sizeof(double));
  Pab            = (double*) malloc((nvir*(nvir+1)/2)*sizeof(double));
  Wkj            = (double*) malloc(nocc*nocc*sizeof(double));
  Wab            = (double*) malloc(nvir*nvir*sizeof(double));
  Waj            = (double*) malloc(nvir*nocc*sizeof(double));
  Laj            = (double*) malloc(nvir*nocc*sizeof(double));


  if (me == 0) allocbn_double_matrix(gradient, "n1 n2", natom, 3);
  allocbn_double_matrix(ginter, "n1 n2", natom, 3);

  //////////////////////////////
  // Initialize various arrays
  //////////////////////////////

  bzerofast(Pkj,nocc*(nocc+1)/2);
  bzerofast(Wkj,nocc*nocc);
  bzerofast(Pab,nvir*(nvir+1)/2);
  bzerofast(Wab,nvir*nvir);
  bzerofast(Waj,nvir*nocc);
  bzerofast(Laj,nvir*nocc);

  if (me == 0) fill_double_matrix(gradient, 0, gradient->n1, 0, gradient->n2, 0.0);
  fill_double_matrix(ginter, 0, ginter->n1, 0, ginter->n2, 0.0);

  // Debug print
  if (me == 0) {
    for (j=0; j<nbasis; j++) {
      fprintf(stdout,"evals,j: %15.10lf %i\n", evals[j],j);
      }
    }
  fflush(outfile);
  // End of debug print

  if (nproc > 1) iqjs_buf = new double[2 + nfuncmax*nbasis];

  /////////////////////////////////////
  //  Begin MP2 loops
  /////////////////////////////////////

  // debug print
  if (me == 0) {
    fprintf(stdout,"node %i, begin loop over i-batches\n",me);
    fflush(outfile);
    }
  // end of debug print


  int nijmax = 0;
  index = 0;
  for (i=0; i<ni; i++) {
      for (j=0; j<nocc; j++) {
          if (index++ % nproc == me) nijmax++;
        }
    }

  RefMemoryGrp mem
      = MemoryGrp::create_memorygrp(nijmax*nbasis*nbasis*sizeof(double));

  mem->lock(0);

  MemoryGrpBuf<double> membuf(mem);
  MemoryGrpBuf<double> membuf_remote(mem);

  // find the start of memory for debugging purposes:
  const double *memstart = membuf.readonly_on_node(0,1,0);
  membuf.release();

  for (pass=0; pass<npass; pass++) {

    i_offset = pass*ni;
    if ((pass == npass - 1) && (rest != 0)) ni = rest;

    // Compute number of of i,j pairs on each node for
    // mo_int, gamma_iajs and more ?
    index = 0;
    nij = 0;
    for (i=0; i<ni; i++) {
      for (j=0; j<nocc; j++) {
        if (index++ % nproc == me) nij++;
        }
      }

    // debug print
    fprintf(stdout,"node %i, nij = %i\n", me, nij);
    fflush(outfile);
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
    if (me == 0) {
      fprintf(stdout,"Begin loop over shells (erep, 1.+2. qt)\n");
      fflush(outfile);
      }
    // end of debug print

    for (S=0; S<nshell; S++) {
      ns = INT_SH_NFUNC((centers),S);
      s_offset = centers->func_num[S];

      for (R=0; R<=S; R++) {
        nr = INT_SH_NFUNC((centers),R);
        r_offset = centers->func_num[R];

        if (index++ % nproc == me) {

          bzerofast(integral_iqrs, ni*nbasis*nfuncmax*nfuncmax);

          for (Q=0; Q<nshell; Q++) {
            nq = INT_SH_NFUNC((centers),Q);
            q_offset = centers->func_num[Q];
            for (P=0; P<=Q; P++) {
              np = INT_SH_NFUNC((centers),P);
              p_offset = centers->func_num[P];

           // if (scf_erep_bound(P,Q,R,S) < tol) {
           //   continue;  // skip ereps less than tol
           //   }
              if (int_erep_4bound(P,Q,R,S) < tol) {
                continue;  // skip ereps less than tol
                }

              aoint_computed++;

              tim_enter("erep");
              int_erep(INT_EREP|INT_NOBCHK|INT_NOPERM|INT_REDUND,&P,&Q,&R,&S);
              tim_exit("erep");

              tim_enter("1. q.t.");
              // Begin first quarter transformation

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
                        iprs_ptr++;
                        iqrs_ptr++;
                        continue; // skip to next bf4 value
                        }

                      if (INT_NONZERO(*pqrs_ptr)) {
                        iprs_ptr = &integral_iqrs[bf4 + ns*(p + nbasis*bf3)];
                        iqrs_ptr = &integral_iqrs[bf4 + ns*(q + nbasis*bf3)];
                        c_qi = &scf_vector[q][i_offset];
                        c_pi = &scf_vector[p][i_offset];
                        for (i=0; i<ni; i++) {
                          *iprs_ptr += *c_qi++**pqrs_ptr;
                          iprs_ptr += offset;
                          if (p != q) {
                            *iqrs_ptr += *c_pi++**pqrs_ptr;
                            iqrs_ptr += offset;
                            }
                          } // exit i loop
                        }   // endif

                      pqrs_ptr++;
                      iprs_ptr++;
                      iqrs_ptr++;
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
          // Begin second quarter transformation
          for (i=0; i<ni; i++) {
            for (j=0; j<nocc; j++) {

              bzerofast(iqjs_contrib, nbasis*nfuncmax);
              bzerofast(iqjr_contrib, nbasis*nfuncmax);

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
    if (me == 0) {
      fprintf(stdout,"End of loop over shells\n");
      fflush(outfile);
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
    integral_iajs = new double[nvir];
    integral_ikjs = new double[nocc];

    // debug print
    if (me == 0) {
    fprintf(stdout,"Begin 3. qt\n");
    fflush(outfile);
      }
    // end of debug print

    tim_enter("3. q.t.");
    // Begin third quarter transformation
    index = 0;
    ij_index = 0;
    for (i=0; i<ni; i++) {
      for (j=0; j<nocc; j++) {
        if (index++ % nproc == me) {

          for (s=0; s<nbasis; s++) {

            bzerofast(integral_iajs, nvir);
            bzerofast(integral_ikjs, nocc);
            for (q=0; q<nbasis; q++) {
              integral_iqjs_ptr = &integral_iqjs[q + nbasis*(s + nbasis*ij_index)];
              iajs_ptr = integral_iajs;
              c_qa = &scf_vector[q][nocc];
              for (a=0; a<nvir; a++) {
                *iajs_ptr++ += *c_qa++ * *integral_iqjs_ptr;
                }
              ikjs_ptr = integral_ikjs;
              c_qk = scf_vector[q];
              for (k=0; k<nocc; k++) {
                *ikjs_ptr++ += *c_qk++ **integral_iqjs_ptr;
                }
              }   // exit q loop

            // Put iajs and ikjs into integral_iqjs, while overwriting what was there
            // i.e., integral_iqjs will now contain three-quarter transformed integrals
            // iajs and ikjs
            integral_iqjs_ptr = &integral_iqjs[nocc + nbasis*(s + nbasis*ij_index)];
            iajs_ptr = integral_iajs;
            for (a=0; a<nvir; a++) {
              *integral_iqjs_ptr++ = *iajs_ptr++;
              }
            integral_iqjs_ptr = &integral_iqjs[nbasis*(s + nbasis*ij_index)];
            ikjs_ptr = integral_ikjs;
            for (k=0; k<nocc; k++) {
              *integral_iqjs_ptr++ = *ikjs_ptr++;
              }
            }   // exit s loop
          ij_index++;
          }     // endif
        }       // exit j loop
      }         // exit i loop
    // end of third quarter transformation
    tim_exit("3. q.t.");

    // debug print
    if (me == 0) {
      fprintf(stdout,"End of 3. qt\n");
      fflush(outfile);
      }
    // end of debug print

    delete[] integral_iajs;
    delete[] integral_ikjs;

    // The array of half-transformed integrals integral_iqjs has now
    // been overwritten by three-quarter transformed integrals iajs,
    // and ikjs; rename the array integral_ixjs, where x = any MO
    integral_ixjs = integral_iqjs;

    integral_iajy = new double[nbasis];
    integral_ikja = new double[nvir];

    // debug print
    if (me == 0) {
      fprintf(stdout,"Begin 4. qt\n");
      fflush(outfile);
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
              for (y=0; y<nbasis; y++) {
                *iajy_ptr++ += *c_sy++ * *iajs_ptr;
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

          for (k=0; k<nocc; k++) {
            bzerofast(integral_ikja, nvir);
            ikjs_ptr = &integral_ixjs[k + nbasis*nbasis*ij_index];
            for (s=0; s<nbasis; s++) {
              c_sa = &scf_vector[s][nocc];
              ikja_ptr = integral_ikja;
              for (a=0; a<nvir; a++) {
                *ikja_ptr++ += *c_sa++ * *ikjs_ptr;
                } // exit a loop 
              ikjs_ptr += nbasis;
              }   // exit s loop 
            // Put integral_ikja into ixjs for one i,k,j while
            // overwriting elements of ixjs
            ikjs_ptr = &integral_ixjs[k + nbasis*(nocc + nbasis*ij_index)];
            ikja_ptr = integral_ikja;
            for (a=0; a<nvir; a++) {
              *ikjs_ptr = *ikja_ptr++;
              ikjs_ptr += nbasis;
              } // exit a loop 
            }     // exit k loop 

        ij_index++;
        }   // endif
      }     // exit j loop
    }       // exit i loop
    // end of fourth quarter transformation
    tim_exit("4. q.t.");

    // debug print
    if (me == 0) {
      fprintf(stdout,"End of 4. qt\n");
      fflush(outfile);
      }
    // end of debug print

    // The array integral_ixjs has now been overwritten by MO integrals
    // iajy and ikja, so rename the array mo_int
    mo_int = integral_ixjs;

    delete[] integral_iajy;
    delete[] integral_ikja;

    // Divide the (ia|jb) MO integrals by the term 
    // evals[i]+evals[j]-evals[a]-evals[b]
    // and keep these integrals in mo_int
    tim_enter("divide (ia|jb)'s");

    index = 0;
    ij_index = 0;
    for (i=0; i<ni; i++) {
      for (j=0; j<nocc; j++) {
        if (index++ % nproc == me) {
          for (b=0; b<nvir; b++) {
            iajb_ptr = &mo_int[nocc + nbasis*(b+nocc + nbasis*ij_index)];
            for (a=0; a<nvir; a++) {
             *iajb_ptr++ /= evals[i+i_offset]+evals[j]-evals[a+nocc]-evals[b+nocc];
              } // exit a loop
            }   // exit b loop
          ij_index++;
          }     // endif
        }       // exit j loop
      }         // exit i loop
    tim_exit("divide (ia|jb)'s");

    // We now have the fully transformed integrals (ia|jb)
    // (divided by the proper orbital energy denominators)
    // for one batch of i, all j<nocc, and all a<nvir and b<nvir;
    // compute contribution to the MP2 correlation energy
    // from these integrals 

    tim_enter("compute ecorr");

    index = 0;
    ij_index = 0;
    for (i=0; i<ni; i++) {
      for (j=0; j<nocc; j++) {
        if (index++ % nproc == me) {

          for (b=0; b<nvir; b++) {
            iajb_ptr = &mo_int[nocc + nbasis*(b+nocc + nbasis*ij_index)];
            ibja_ptr = &mo_int[b+nocc + nbasis*(nocc + nbasis*ij_index)];
            for (a=0; a<nvir; a++) {
              delta_ijab = evals[i_offset+i]+evals[j]-evals[nocc+a]-evals[nocc+b];
              ecorr_mp2 += *iajb_ptr*(2**iajb_ptr - *ibja_ptr)*delta_ijab;
              iajb_ptr++;
              ibja_ptr += nbasis;;
              } // exit a loop
            }   // exit b loop

          ij_index++;
          }     // endif
        }       // exit j loop
      }         // exit i loop
    tim_exit("compute ecorr");

    // debug print
    if (me == 0) {
    fprintf(stdout,"End of ecorr\n");
    fflush(outfile);
      }
    // end of debug print

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
          for (kloop=me; kloop<me+nocc; kloop++) {
            // stagger k's to minimize contention
            k = kloop%nocc;
            if (k>=j) pkj_ptr = &Pkj[k*(k+1)/2 + j];
            wjk_ptr = &Wkj[j*nocc + k];
            // Send for iakb, if necessary
            ik_index = (i*nocc + k)/nproc;
            ik_proc = (i*nocc + k)%nproc;
            ik_offset = nocc + nocc*nbasis + nbasis*nbasis*ik_index;
            mo_intbuf = (double*) membuf_remote.readonly_on_node(ik_offset,
                                                                 nbasis*nvir-nocc,
                                                                 ik_proc);
            for (a=0; a<nvir; a++) {
              ibja_ptr = &mo_int[nocc + nbasis*(a+nocc + nbasis*ij_index)];
              iajb_ptr = &mo_int[a+nocc + nbasis*(nocc + nbasis*ij_index)];
              /* if (ik_proc == me) iakb_ptr = &mo_int[a + ik_offset];
              else */
                iakb_ptr = &mo_intbuf[a];
              for (b=0; b<nvir; b++) {
                tmpval = 2**iakb_ptr * (*ibja_ptr++ - 2 * *iajb_ptr);
                iakb_ptr += nbasis;
                iajb_ptr += nbasis;
                if (k>=j) *pkj_ptr += tmpval;
                delta_ijab = evals[i_offset+i]+evals[j]-evals[nocc+a]-evals[nocc+b];
                *wjk_ptr += tmpval*delta_ijab;
                } // exit b loop
              }   // exit a loop
            mo_intbuf = 0;
            membuf_remote.release();
            }     // end kloop loop
          ij_index++;
          }       // endif
        }         // exit j loop
      }           // exit i loop
    tim_exit("Pkj and Wkj");

    // debug print
    if (me == 0) {
      fprintf(stdout,"End of Pkj and Wkj\n");
      fflush(outfile);
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
          offset = nocc + nocc*nbasis + nbasis*nbasis*ij_index;
          for (a=0; a<nvir; a++) {
            pab_ptr = &Pab[a*(a+1)/2];
            for (b=0; b<=a; b++) {
              wab_ptr = &Wab[a*nvir + b];
              wba_ptr = &Wab[b*nvir + a];
              ibjc_ptr = &mo_int[offset + b];
              icjb_ptr = &mo_int[offset + b*nbasis];
              iajc_ptr = &mo_int[offset + a];
              icja_ptr = &mo_int[offset + a*nbasis];
              for (c=0; c<nvir; c++) {
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
          ij_index++;
          }     // endif
        }       // exit j loop
      }         // exit i loop
    tim_exit("Pab and Wab");

    // debug print
    if (me == 0) {
    fprintf(stdout,"End of Pab and Wab\n");
    fflush(outfile);
      }
    // end of debug print

    ///////////////////////////////////////
    // Update Waj and Laj with contrib. 
    // from (oo|ov) and (ov|oo) integrals
    ///////////////////////////////////////
    tim_enter("Waj and Laj");

    // (oo|ov) contribution
    index = 0;
    ik_index = 0;
    for (i=0; i<ni; i++) {
      for (k=0; k<nocc; k++) {
        if (index++ % nproc == me) {
          offset = nbasis*nocc + nbasis*nbasis*ik_index;
          for (j=0; j<nocc; j++) {
            for (b=0; b<nvir; b++) {
              ibka_ptr = &mo_int[b+nocc + offset];
              ijkb_ptr = &mo_int[j + nbasis*b + offset];
              waj_ptr = &Waj[j*nvir]; // order as j*nvir+a to make loops more efficient
              laj_ptr = &Laj[j*nvir];
              for (a=0; a<nvir; a++) {
                tmpval = 2**ibka_ptr * *ijkb_ptr;
                ibka_ptr += nbasis;
                *waj_ptr++ += tmpval;
                *laj_ptr++ -= tmpval; // This term had the wrong sign in Frisch's paper
                } // exit a loop
              }   // exit b loop
            }     // exit j loop
          ik_index++;
          }       // endif
        }         // exit k loop
      }           // exit i loop

    // (ov|oo) contribution
    index = 0;
    ik_index = 0;
    for (i=0; i<ni; i++) {
      for (k=0; k<nocc; k++) {
        if (index++ % nproc == me) {
          offset = nocc + nbasis*nbasis*ik_index;
          for (b=0; b<nvir; b++) {
            for (j=0; j<nocc; j++) {
              ibkj_ptr = &mo_int[offset + b + j*nbasis];
              ibka_ptr = &mo_int[offset + b + nocc*nbasis];
              waj_ptr = &Waj[j*nvir];
              laj_ptr = &Laj[j*nvir];
              for (a=0; a<nvir; a++) {
                tmpval = 4 * *ibka_ptr * *ibkj_ptr;
                ibka_ptr += nbasis;
                *waj_ptr++ -= tmpval;
                *laj_ptr++ += tmpval; // This term had the wrong sign in Frisch's paper
                } // exit a loop
              }   // exit j loop
            }     // exit b loop
          ik_index++;
          }       // endif
        }         // exit k loop
      }           // exit i loop

    tim_exit("Waj and Laj");
    /////////////////////////////
    // End of Waj and Laj update
    /////////////////////////////

    // debug print
    if (me == 0) {
      fprintf(stdout,"End of Paj and Waj\n");
      fflush(outfile);
      }
    // end of debug print

    mo_int = 0;
    membuf.release();

    mem->sync(); // Need to synchronize before deleting mo_intbuf

    mo_int = membuf.readwrite_on_node(0, nij*nbasis*nbasis);

    gamma_iajs_tmp = new double[nbasis*nvir];
    if (!gamma_iajs_tmp) {
      fprintf(outfile,"Could not allocate gamma_iajs_tmp\n");
      abort();
      }

    // debug print
    if (me == 0) {
      fprintf(stdout,"Begin 1+2qbt\n");
      fflush(outfile);
      }
    // end of debug print

    /////////////////////////////////////////////////////////
    // Perform first and second quarter back-transformation.
    // Each node produces gamma_iajs, and gamma_iqjs 
    // for a subset of i and j, all a and all s
    /////////////////////////////////////////////////////////

    // Begin first quarter back-transformation
    tim_enter("1. q.b.t.");
    index = 0;
    ij_index = 0;
    for (i=0; i<ni; i++) {
      for (j=0; j<nocc; j++) {
        if (index++ % nproc == me) {
          bzerofast(gamma_iajs_tmp,nbasis*nvir);
          offset = nocc + nocc*nbasis + nbasis*nbasis*ij_index;

          for (a=0; a<nvir; a++) {
            for (s=0; s<nbasis; s++) {
              c_sb = &scf_vector[s][nocc];
              gamma_iajs_ptr = &gamma_iajs_tmp[s*nvir + a];
              ibja_ptr = &mo_int[a*nbasis + offset];
              iajb_ptr = &mo_int[a + offset];

              for (b=0; b<nvir; b++) {
                *gamma_iajs_ptr += 2**c_sb++ * (2**iajb_ptr - *ibja_ptr++);
                iajb_ptr += nbasis;
                } // exit b loop
              }   // exit s loop
            }     // exit a loop
          // Put gamma_iajs_tmp into mo_int for one i,j
          // while overwriting mo_int
          gamma_iajs_ptr = gamma_iajs_tmp;
          for (y=0; y<nbasis; y++) {
            iajy_ptr = &mo_int[nocc + nbasis*(y + nbasis*ij_index)];
            for (a=0; a<nvir; a++) {
              *iajy_ptr++ = *gamma_iajs_ptr++;
              }
            }

          ij_index++;
          }       // endif
        }         // exit j loop
      }           // exit i loop
    // end of first quarter back-transformation
    tim_exit("1. q.b.t.");

    // debug print
    if (me == 0) {
      fprintf(stdout,"End 1+2qbt\n");
      fflush(outfile);
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
      fprintf(outfile,"Could not allocate gamma_iqjs_tmp\n");
      abort();
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
          offset = nbasis*nbasis*ij_index;

          for (s=0; s<nbasis; s++) {
            bzerofast(gamma_iqjs_tmp,nbasis);
            for (q=0; q<nbasis; q++) {
              gamma_iqjs_ptr = &gamma_iqjs_tmp[q];
              gamma_iajs_ptr = &gamma_iajs[nocc + s*nbasis + offset];
              c_qa = &scf_vector[q][nocc];

              for (a=0; a<nvir; a++) {
                *gamma_iqjs_ptr += *c_qa++ * *gamma_iajs_ptr++;
                } // exit a loop
              }   // exit q loop
            // Put gamma_iqjs_tmp into gamma_iajs for one i,j,s
            // while overwriting gamma_iajs
            gamma_iajs_ptr = &gamma_iajs[s*nbasis + offset];
            gamma_iqjs_ptr = gamma_iqjs_tmp;
            for (q=0; q<nbasis; q++) {
              *gamma_iajs_ptr++ = *gamma_iqjs_ptr++;
              }

            }     // exit s loop

          ij_index++;
          }       // endif
        }         // exit j loop
      }           // exit i loop
    tim_exit("2. q.b.t.");
    // end of second quarter back-transformation

    gamma_iajs = 0;
    membuf.release();
    
    mem->sync(); // Keep this here to make sure all nodes have gamma_iqjs
                 // before it is needed below, and that gamma_iajs is not
                 // deleted prematurely

    // The quarter back-transformed elements gamma_iajs have now been
    // overwritten by the half back-transformed elements gamma_iqjs, so rename
    gamma_iqjs = (double*) membuf.readonly_on_node(0, nij*nbasis*nbasis);

    delete[] gamma_iqjs_tmp;

    /////////////////////////////////////////////////
    // End of 1. and 2. quarter back-transformation
    /////////////////////////////////////////////////


    if (nproc > 1) {
      gammabuf = new double[nbasis*nfuncmax];
      if (!gammabuf) {
        fprintf(outfile,"Could not allocate gammabuf\n");
        abort();
        }
      }
 
    // Allocate various arrays
    gamma_iqrs = new double[ni*nbasis*nfuncmax*nfuncmax];
    if (!gamma_iqrs) {
      fprintf(outfile,"Could not allocate gamma_iqrs\n");
      abort();
      }
    gamma_pqrs = new double[nfuncmax*nfuncmax*nfuncmax*nfuncmax];
    if (!gamma_pqrs) {
      fprintf(outfile,"Could not allocate gamma_pqrs\n");
      abort();
      }

    Lpi = new double[nbasis*ni];
    bzerofast(Lpi,nbasis*ni);

    ////////////////////////////////////////////////////////
    // Perform third and fourth quarter back-transformation
    // and compute contrib. to gradient from non-sep 2PDM
    ////////////////////////////////////////////////////////

    index = 0;

    for (S=0; S<nshell; S++) {
      ns = INT_SH_NFUNC((centers),S);
      s_offset = centers->func_num[S];

      for (R=0; R<=S; R++) {
        nr = INT_SH_NFUNC((centers),R);
        r_offset = centers->func_num[R];

        // If both PQRS and PQRS derivative are zero, skip to next R
        if (int_erep_2bound(R,S) < tol && int_erep_2bound_1der(R,S) < tol) continue;

        if (index++ % nproc == me) {

          tim_enter("3. q.b.t.");
          // Begin third quarter back-transformation.

          bzerofast(gamma_iqrs,ni*nbasis*nfuncmax*nfuncmax);

          for (i=0; i<ni; i++) {
            for (jloop=me; jloop<me+nocc; jloop++) {
              // stagger j's to minimize contention
              j = jloop%nocc;

              ij_proc =  (i*nocc + j)%nproc; // ij_proc has this ij pair
              ij_index = (i*nocc + j)/nproc;
              offset = s_offset*nbasis + ij_index*nbasis*nbasis;

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
              // Send for gamma_irjq, if necessary (in array gamma_iqjs)

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

          if (int_erep_2bound(R,S) >= tol) {  // only do this if integral is nonzero

            // Compute contrib to Laj from (ov|vv) integrals
            // (done in AO basis to avoid generating (ov|vv)
            tim_enter("(ov|vv) contrib to Laj");
            for (Q=0; Q<nshell; Q++) {
              nq = INT_SH_NFUNC((centers),Q);
              q_offset = centers->func_num[Q];
              for (P=0; P<=Q; P++) {
                np = INT_SH_NFUNC((centers),P);
                p_offset = centers->func_num[P];
             // if (scf_erep_bound(P,Q,R,S) < tol) {
             //   continue;  // skip ereps less than tol
             //   }
                if (int_erep_4bound(P,Q,R,S) < tol) {
                  continue;  // skip ereps less than tol
                  }
                tim_enter("erep");
                int_erep(INT_EREP|INT_NOBCHK|INT_NOPERM|INT_REDUND,&P,&Q,&R,&S);
                tim_exit("erep");

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

                        if (INT_NONZERO(intbuf[int_index])) {
                          s = s_offset + bf4;

                          if (s < r) {
                            int_index++;
                            continue; // skip to next bf4 value
                            }

                          factor_rs = (r==s ? 0:1);

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

          if (int_erep_2bound_1der(R,S) >= tol) {

            for (Q=0; Q<=S; Q++) {
              nq = INT_SH_NFUNC((centers),Q);

              for (P=0; P<=(Q==S ? R:Q); P++) {
                np = INT_SH_NFUNC((centers),P);

                // If integral derivative is less than threshold skip to next P
                if (int_erep_4bound_1der(P,Q,R,S) < tol) continue;
                aointder_computed++;

                tim_enter("4. q.b.t.");
                bzerofast(gamma_pqrs,nfuncmax*nfuncmax*nfuncmax*nfuncmax);

                offset = nr*ns*nbasis;

                // Begin fourth quarter back-transformation
                gamma_pqrs_ptr = gamma_pqrs;
                for (bf1=0; bf1<np; bf1++) {
                  p = bf1 + centers->func_num[P];
                  for (bf2=0; bf2<nr; bf2++) {
                    for (bf3=0; bf3<nq; bf3++) {
                      q = bf3 + centers->func_num[Q];
                      for (bf4=0; bf4<ns; bf4++) {
                        c_pi = &scf_vector[p][i_offset];
                        c_qi = &scf_vector[q][i_offset];
                        gamma_iqrs_ptr = &gamma_iqrs[bf4 + ns*(q + nbasis*bf2)];
                        gamma_iprs_ptr = &gamma_iqrs[bf4 + ns*(p + nbasis*bf2)];
                        for (i=0; i<ni; i++) {
                          *gamma_pqrs_ptr += *c_pi * *gamma_iqrs_ptr;
                          if (p!=q) *gamma_pqrs_ptr += *c_qi * *gamma_iprs_ptr;
                          c_pi++;
                          c_qi++;
                          gamma_iqrs_ptr += offset;
                          gamma_iprs_ptr += offset;
                          } // exit i loop
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
                int_erep_all1der(flags,&P,&Q,&R,&S,&der_centers);
                tim_exit("erep derivs");

                // Compute contribution to gradient from non-sep 2PDM
                // (i.e., contract derivative integrals with gamma_pqrs)
                int_index = 0;
                tim_enter("non-sep 2PDM contrib.");
                for (derset=0; derset<der_centers.n; derset++) {
                  for (xyz=0; xyz<3; xyz++) {
                    grad_ptr1 = &ginter->d[der_centers.num[derset]][xyz];
                    grad_ptr2 = &ginter->d[der_centers.onum][xyz];
                    for (bf1=0; bf1<np; bf1++) {
                      p = bf1 + centers->func_num[P];
                      for (bf2=0; bf2<nq; bf2++) {
                        q = bf2 + centers->func_num[Q];
                        if (p == q) factor_pq = 0;
                        else factor_pq = 1;
                        qp = q*(q+1)/2 + p;
                        for (bf3=0; bf3<nr; bf3++) {
                          r = bf3 + centers->func_num[R];
                          gamma_pqrs_ptr = &gamma_pqrs[ns*(bf2 + nq*(bf3 + nr*bf1))];
                          for (bf4=0; bf4<ns; bf4++) {
                            s = bf4 + centers->func_num[S];
                            if (r == s) factor_rs = 0;
                            else factor_rs = 1;
                            sr = s*(s+1)/2 + r;
                            if (q == s && p == r) factor_pqrs = 1;
                            else factor_pqrs = 2;
                            tmpval = intbuf[int_index]*factor_pqrs**gamma_pqrs_ptr;
                            gamma_pqrs_ptr++;
                            if (q>=p && s>=r && (P != R || Q != S || sr >= qp)) {
                              *grad_ptr1 += tmpval;
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

    mem->sync(); // Make sure all nodes are done before deleting arrays

    delete[] Lpi;

    delete[] gamma_iqrs;
    delete[] gamma_pqrs;

    gamma_iqjs = 0;
    membuf.release();

    }           // exit loop over i-batches (pass)

  mem = 0;

  // debug print
  if (me == 0) {
    fprintf(stdout,"Exited loop over i-batches\n");
    fflush(outfile);
    }
  // end of debug print

  if (nproc > 1) delete[] iqjs_buf;

  // Accumulate intermediate gradients on node 0
  sum_gradients(ginter, nproc);

  // Add intermediate gradients to the gradient on node 0
  if (me == 0) {
    add_double_matrix(gradient, 0, gradient->n1, 0, gradient->n2,
                      ginter, 0, gradient->n1, 0, gradient->n2);
    }

  // Print out contribution to the gradient from non-sep. 2PDM
  if (me == 0) {
    fprintf(outfile,"Contribution to MP2 gradient from non-separable 2PDM [au]:\n");
    for (i=0; i<natom; i++) {
      fprintf(outfile,"%15.10lf  %15.10lf  %15.10lf\n",
              ginter->d[i][0], ginter->d[i][1], ginter->d[i][2]);
      }
    }

  ///////////////////////////////////////////////////////////////
  // The computation of the MP2 energy is now complete on each
  // node; add the nodes' contributions and print out the energy
  ///////////////////////////////////////////////////////////////
  gsum0(&ecorr_mp2,1,5,mtype_get(),0);
  gsum0(&aoint_computed,1,2,mtype_get(),0);
  gsum0(&aointder_computed,1,2,mtype_get(),0);

  if (me == 0) {
    escf = scf_info->e_elec + scf_info->nuc_rep;
    emp2 = escf + ecorr_mp2;

    // Print out various energies etc.

    fprintf(outfile,"Number of shell quartets for which AO integrals \n"
                    "(or integral derivatives) would have been computed\n"
                    "without bounds checking: %i\n",
                     npass*nshell*nshell*(nshell+1)*(nshell+1)/2);
    fprintf(outfile,"Number of shell quartets for which AO integrals\n"
                    "were computed: %i\n",aoint_computed);
    fprintf(outfile,"Number of shell quartets for which AO integral derivatives\n"
                    "were computed: %i\n",aointder_computed);

    fprintf(outfile,"ROHF energy [au]:                  %13.8lf\n", escf);
    fprintf(outfile,"MP2 correlation energy [au]:       %13.8lf\n", ecorr_mp2);
    fprintf(outfile,"MP2 energy [au]:                   %13.8lf\n", emp2);
    fflush(outfile);
    }


  ////////////////////////////////////////////////////////
  // Add contributions from all nodes to various matrices
  ////////////////////////////////////////////////////////
  tmpint = (nvir > nocc ? nvir:nocc);
  double *tmpmat = new double[tmpint*tmpint];
  gop1(Laj,nvir*nocc,tmpmat,'+',3);
  gop1(Pkj,nocc*(nocc+1)/2,tmpmat,'+',3); // Pkj is now complete
  gop1(Pab,nvir*(nvir+1)/2,tmpmat,'+',3); // Pab is now complete
  gop1(Wab,nvir*nvir,tmpmat,'+',3);
  gop1(Wkj,nocc*nocc,tmpmat,'+',3);
  gop1(Waj,nvir*nocc,tmpmat,'+',3);
  delete[] tmpmat;

  RefSCDimension nocc_dim(new LocalSCDimension(nocc));
  RefSCDimension nvir_dim(new LocalSCDimension(nvir));
  RefSCDimension nbasis_dim(new LocalSCDimension(nbasis));


  // Finish computation of Wab
  tim_enter("Pab and Wab");
  pab_ptr = Pab;
  for (a=0; a<nvir; a++) {
    wba_ptr = &Wab[a];
    wab_ptr = &Wab[a*nvir];
    for (b=0; b<=a; b++) {
      if (a==b) {
        *wab_ptr++ -= evals[nocc+a]**pab_ptr++;
        }
      else {
        *wab_ptr++ -= evals[nocc+a]**pab_ptr;
        *wba_ptr   -= evals[nocc+b]**pab_ptr;
        pab_ptr++;
        }
      wba_ptr += nvir;
      } // exit b loop
    }   // exit a loop
  // Wab is now complete
  tim_exit("Pab and Wab");
  RefSCMatrix Wab_matrix(nvir_dim, nvir_dim);
  Wab_matrix->assign(Wab); // Put elements of Wab into Wab_matrix
  free(Wab);


  // Update Wkj with contribution from Pkj
  tim_enter("Pkj and Wkj");
  pkj_ptr = Pkj;
  for (k=0; k<nocc; k++) {
    wjk_ptr = &Wkj[k];
    wkj_ptr = &Wkj[k*nocc];
    for (j=0; j<=k; j++) {
      if (j==k) {
        *wkj_ptr++ -= evals[k]**pkj_ptr++;
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
  RefSCMatrix Cv(nbasis_dim, nvir_dim); // virtual block of scf_vector
  RefSCMatrix Co(nbasis_dim, nocc_dim); // occupied block of scf_vector
  for (p=0; p<nbasis; p++) {
    c_pq = scf_vector[p];
    for (q=0; q<nbasis; q++) {
      if (q<nocc) Co->set_element(p, q, *c_pq++);
      else Cv->set_element(p, q-nocc, *c_pq++);
      }
    }


  // Compute the density-like matrix Dmat_matrix
  RefSymmSCMatrix Pab_matrix(nvir_dim);
  RefSymmSCMatrix Pkj_matrix(nocc_dim);
  RefSCMatrix Dmat_matrix(nbasis_dim,nbasis_dim);
  Pab_matrix->assign(Pab); // fill in elements of Pab_matrix from Pab
  free(Pab);
  Pkj_matrix->assign(Pkj); // fill in elements of Pkj_matrix from Pkj
  free(Pkj);
  Dmat_matrix = Cv*Pab_matrix*Cv.t() + Co*Pkj_matrix*Co.t();
  // We now have the density-like matrix Dmat_matrix


  // Need to synchronize all nodes here ?

  // Compute the G matrix
  Dmat = new double[nbasis*nbasis];
  Dmat_matrix->convert(Dmat); // convert Dmat_matrix to Dmat (double*)

  RefSymmSCMatrix Gmat(nbasis_dim);
  mbpt_init_gmat(centers, scf_info, intbuf);
  tim_enter("make_gmat for Laj");
  mbpt_make_gmat(scf_info, centers, Gmat, Dmat, outfile);
  tim_exit("make_gmat for Laj");

  // Finish computation of Laj
  RefSCMatrix Laj_matrix(nocc_dim,nvir_dim); // elements are ordered as j*nvir+a
  Laj_matrix->assign(Laj);
  Laj_matrix = Laj_matrix - 2*Co.t()*Gmat*Cv;
  Laj_matrix->convert(Laj);  // Put new Laj_matrix elements into Laj

  tim_exit("Laj");

  //////////////////////////////////////
  // Computation of Laj is now complete
  //////////////////////////////////////

  ////////////////////////////
  // Solve the CPHF equations
  ////////////////////////////
  RefSCMatrix Paj_matrix(nvir_dim, nocc_dim);
  tim_enter("cphf");
  mbpt_cphf(centers, scf_info, outfile, nbasis, nvir, nocc, scf_vector,
            Laj, evals, Paj_matrix);
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
  RefSCMatrix Waj_matrix(nocc_dim, nvir_dim);
  Waj_matrix->assign(Waj); // Put elements of Waj into Waj_matrix
  // NB. Waj_matrix elements are ordered as j*nvir+a
  free(Waj);


  // Finish computation of Wkj
  tim_enter("Pkj and Wkj");
  Dmat_matrix = Co*(Pkj_matrix*Co.t() + Paj_matrix.t()*Cv.t()) +
                Cv*(Paj_matrix*Co.t() + Pab_matrix*Cv.t());
  Dmat_matrix->convert(Dmat); // convert Dmat_matrix to Dmat (double*)
  tim_enter("make_gmat for Wkj");
  mbpt_make_gmat(scf_info, centers, Gmat, Dmat, outfile);
  tim_exit("make_gmat for Wkj");
  mbpt_done_gmat(centers, scf_info);
  RefSCMatrix Wkj_matrix(nocc_dim, nocc_dim);
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
  RefSCMatrix P2AO_matrix(nbasis_dim, nbasis_dim);
  RefSCMatrix P2MO_matrix(nbasis_dim, nbasis_dim);
  RefSCMatrix W2AO_matrix(nbasis_dim, nbasis_dim);
  RefSCMatrix W2MO_matrix(nbasis_dim, nbasis_dim);
  RefSCMatrix SCF_matrix(nbasis_dim, nbasis_dim);
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

  ////////////////////////////////////////////////
  // Compute the contribution to the MP2 gradient 
  // from the separable part of the 2PDM
  ////////////////////////////////////////////////

  fill_double_matrix(ginter, 0, ginter->n1, 0, ginter->n2, 0.0);
  tim_enter("sep 2PDM contrib.");
  s2pdm_contrib(intbuf, centers, PHF, P2AO, ginter, nproc, me);
  tim_exit("sep 2PDM contrib.");
  delete[] PHF;
  delete[] P2AO;

  // The separable 2PDM contribution to the gradient has now been
  // accumulated in ginter on node 0; add it to the total gradients
  if (me == 0) {
    add_double_matrix(gradient, 0, gradient->n1, 0, gradient->n2,
                      ginter, 0, gradient->n1, 0, gradient->n2);
    }
  // Print out the contribution to the gradient from sep. 2PDM
  if (me == 0) {
    fprintf(outfile,"Contribution from separable 2PDM to MP2 gradient [au]:\n");
    for (i=0; i<natom; i++) {
      fprintf(outfile,"%15.10lf  %15.10lf  %15.10lf\n",
              ginter->d[i][0], ginter->d[i][1], ginter->d[i][2]);
      }
    }

  // Done with two-electron routines
  int_done_erep();
  int_done_offsets2(centers,centers,centers,centers);
  int_done_bounds_1der();

  /////////////////////////////////////////////////////////////
  // Compute the one-electron contribution to the MP2 gradient
  /////////////////////////////////////////////////////////////

  fill_double_matrix(ginter, 0, ginter->n1, 0, ginter->n2, 0.0);
  tim_enter("hcore contrib.");
  hcore_grad(centers, PMP2, ginter, nproc, me);
  tim_exit("hcore contrib.");
  delete[] PMP2;
  // The hcore contribution to the gradient has now been accumulated
  // in ginter on node 0; add it to the total gradients
  if (me == 0) {
    add_double_matrix(gradient, 0, gradient->n1, 0, gradient->n2,
                      ginter, 0, gradient->n1, 0, gradient->n2);
    }
  // Print out the contribution to the gradient from hcore
  if (me == 0) {
    fprintf(outfile,"Contribution to MP2 gradient from hcore [au]:\n");
    for (i=0; i<natom; i++) {
      fprintf(outfile,"%15.10lf  %15.10lf  %15.10lf\n",
              ginter->d[i][0], ginter->d[i][1], ginter->d[i][2]);
      }
    }

  fill_double_matrix(ginter, 0, ginter->n1, 0, ginter->n2, 0.0);
  tim_enter("overlap contrib.");
  overlap_grad(centers, WMP2, ginter, nproc, me);
  tim_exit("overlap contrib.");
  delete[] WMP2;
  // The overlap contribution to the gradient has now been accumulated
  // in ginter on node 0; add it to the total gradients
  if (me == 0) {
    add_double_matrix(gradient, 0, gradient->n1, 0, gradient->n2,
                      ginter, 0, gradient->n1, 0, gradient->n2);
    }
  // Print out the overlap contribution to the gradient
  if (me == 0) {
    fprintf(outfile,"Overlap contribution to MP2 gradient [au]:\n");
    for (i=0; i<natom; i++) {
      fprintf(outfile,"%15.10lf  %15.10lf  %15.10lf\n",
              ginter->d[i][0], ginter->d[i][1], ginter->d[i][2]);
      }
    }

  ////////////////////////////////////////////////////////
  // Compute the nuclear contribution to the MP2 gradient
  ////////////////////////////////////////////////////////

  if (me == 0) {
    fill_double_matrix(ginter, 0, ginter->n1, 0, ginter->n2, 0.0);
    allocbn_double_vector(&repder,"n",3);
    for (i=0; i<natom; i++) {
      int_nuclear_repulsion_1der(centers,centers,&repder,centers,i);
      for (xyz=0; xyz<3; xyz++) {
        ginter->d[i][xyz] += repder.d[xyz];
        }
      }
    free_double_vector(&repder);
    add_double_matrix(gradient, 0, gradient->n1, 0, gradient->n2,
                      ginter, 0, gradient->n1, 0, gradient->n2);

    // Print out the nuclear contribution to the gradient
    fprintf(outfile,"Nuclear contribution to MP2 gradient [au]:\n");
    for (i=0; i<natom; i++) {
      fprintf(outfile,"%15.10lf  %15.10lf  %15.10lf\n",
              ginter->d[i][0], ginter->d[i][1], ginter->d[i][2]);
      }
    }


  ////////////////////////////////////////////////////////
  // The computation of the MP2 gradient is now complete;
  // print out the gradient
  ////////////////////////////////////////////////////////
  if (me == 0) {
    fprintf(outfile,"Total MP2 gradient [au]:\n");
    for (i=0; i<natom; i++) {
      fprintf(outfile,"%15.10lf  %15.10lf  %15.10lf\n",
              gradient->d[i][0], gradient->d[i][1], gradient->d[i][2]);
      }
    fflush(outfile);
    }

  if (me == 0) {
      for (i=0; i<natom; i++) {
          for (j=0; j<3; j++) {
              gradientt->d[j][i] = gradient->d[i][j];
            }
        }
      free_double_matrix(gradient);
      delete gradient;
    }

  free_double_matrix(ginter);
  delete ginter;

  for (i=0; i<nbasis; i++) {
      delete[] scf_vector[i];
    }
  delete[] scf_vector;

  return(0);

  }

///////////////////////////////////////////////////////////
// Compute the contribution to the MP2 gradient from hcore
///////////////////////////////////////////////////////////
static void 
hcore_grad(centers_t *centers, double *PMP2, double_matrix_t *ginter,
           int nproc, int me)
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
  int natom;

  double *oneebuf; // 1-electron buffer
  double tmpval1, tmpval2;
  double gxyz[3];

  // Initialize 1e routines
  int_initialize_offsets1(centers,centers);
  oneebuf = int_initialize_1e(0,1,centers,centers);

  nshell = centers->nshell;
  nbasis = centers->nfunc;
  natom  = centers->n;

  ///////////////////////////////////////////////////////////
  // Compute the kinetic energy contribution to the gradient
  ///////////////////////////////////////////////////////////

  for (i=0; i<centers->n; i++) {
    jk_index = 0;

    for (j=0; j<nshell; j++) {
      jsize = INT_SH_NFUNC((centers),j);
      j_offset = centers->func_num[j];

      for (k=0; k<=j; k++) {
        ksize = INT_SH_NFUNC((centers),k);
        k_offset = centers->func_num[k];

        if (jk_index++%nproc == me) {
          int_shell_kinetic_1der(centers,centers,oneebuf,j,k,centers,i);
          int_accum_shell_nuclear_hf_1der(centers,centers,oneebuf,
                                          j,k,centers,i);
          int_accum_shell_nuclear_nonhf_1der(centers,centers,oneebuf,
                                             j,k,centers,i);

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

          for (l=0; l<3; l++) ginter->d[i][l] += gxyz[l];

          } // exit "if"
        }   // exit k loop
      }     // exit j loop
    }       // exit i loop

  /* Accumulate the nodes' intermediate gradients on node 0 */
  sum_gradients(ginter, nproc);

  int_done_1e();
  int_done_offsets1(centers,centers);

}


////////////////////////////////////////////////////
// Compute the overlap contribution to the gradient
////////////////////////////////////////////////////
static void 
overlap_grad(centers_t *centers, double *WMP2, double_matrix_t *ginter,
             int nproc, int me)
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
  int natom;

  double *oneebuf; // 1-electron buffer
  double tmpval1, tmpval2;
  double gxyz[3];

  // Initialize 1e routines
  int_initialize_offsets1(centers,centers);
  oneebuf = int_initialize_1e(0,1,centers,centers);

  nshell = centers->nshell;
  nbasis = centers->nfunc;
  natom  = centers->n;

  for (i=0; i<centers->n; i++) {
    jk_index = 0;

    for (j=0; j<nshell; j++) {
      j_offset = centers->func_num[j];
      jsize = INT_SH_NFUNC((centers),j);

      for (k=0; k<=j; k++) {
        k_offset = centers->func_num[k];
        ksize = INT_SH_NFUNC((centers),k);

        if (jk_index++%nproc == me) {
          int_shell_overlap_1der(centers,centers,oneebuf,j,k,centers,i);

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

          for (l=0; l<3; l++) ginter->d[i][l] += gxyz[l];
          } // exit "if"
        }   // exit k loop
      }     // exit j loop
    }       // exit i loop
  
  /* Accumulate the nodes' intermediate gradients on node 0 */
  sum_gradients(ginter, nproc);

  int_done_1e();
  int_done_offsets1(centers,centers);

}


//////////////////////////////////////////////////////////////
// Compute (in the AO basis) the contribution to the gradient 
// from the separable part of the two particle density matrix
// (NB: intbuf and centers must have been initialized)
//////////////////////////////////////////////////////////////
static void
s2pdm_contrib(double *intbuf, centers_t *centers, double *PHF, double *P2AO, 
               double_matrix_t *ginter, int nproc, int me)
{               

  int P, Q, R, S;
  int QP, SR;
  int p, q, r, s;
  int np, nq, nr, ns;
  int p_offset, q_offset, r_offset, s_offset;
  int bf1, bf2, bf3, bf4;
  int index;
  int xyz;
  int derset;
  int nshell = centers->nshell;
  int nbasis = centers->nfunc;
  int natom  = centers->n;
  int flags;

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
  double *integral_ptr;

  int tol = (int) (-10.0/log10(2.0));  // discard erep derivatives smaller than 10^-10

  gammasym_pqrs = new double[nfuncmax*nfuncmax*nfuncmax*nfuncmax];

  der_centers_t der_centers;

  flags = INT_EREP|INT_REDUND|INT_NOPERM|INT_NOSTRB;

  index = 0;

  tim_enter("PQRS loop");

  for (Q=0; Q<nshell; Q++) {
    nq = INT_SH_NFUNC((centers),Q);
    q_offset = centers->func_num[Q];

    for (S=0; S<=Q; S++) {
      ns = INT_SH_NFUNC((centers),S);
      s_offset = centers->func_num[S];

      for (R=0; R<=S; R++) {
        nr = INT_SH_NFUNC((centers),R);
        r_offset = centers->func_num[R];
        k_SR = (R == S ? 0.5 : 1.0);
        SR = S*(S+1)/2 + R;

        for (P=0; P<=(S==Q ? R:Q); P++) {
          // If integral derivative is 0, skip to next P
          if (int_erep_4bound_1der(P,Q,R,S) < tol) continue;

          index++;

          if (index%nproc == me) {
            np = INT_SH_NFUNC((centers),P);
            p_offset = centers->func_num[P];
            k_QP = (P == Q ? 0.5 : 1.0);
            QP = Q*(Q+1)/2 + P;
            k_QPSR = (QP == SR ? 0.5 : 1.0);
            gamma_factor = k_QP*k_SR*k_QPSR;

            // Evaluate derivative integrals
            int_erep_all1der(flags,&P,&Q,&R,&S,&der_centers);

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
                    s = s_offset + bf4;

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
            integral_ptr = intbuf;
            for (derset=0; derset<der_centers.n; derset++) {

              for (xyz=0; xyz<3; xyz++) {
                grad_ptr1 = &ginter->d[der_centers.num[derset]][xyz];
                grad_ptr2 = &ginter->d[der_centers.onum][xyz];

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
  sum_gradients(ginter, nproc);

}

static void
sum_gradients(double_matrix_t *f, int nproc)
{
  int i;

  if (nproc == 1) return;

  for (i=0; i<f->n1; i++) {
    gsum0(f->d[i],f->n2,5,mtype_get(),0);
    }
}

/////////////////////////////////////
// Compute required (dynamic) memory
// for a given batch size
//
// Only arrays allocated before exiting the loop over
// i-batches are included here  - only these arrays
// affect the batch size.
/////////////////////////////////////
static int
compute_batchsize(int mem_alloc, int mem_static, int nproc, FILE* outfile)
{
  int index;
  int mem1, mem2, mem3;
  int mem_dyn;   // dynamic memory available
  int maxdyn;
  int tmp;
  int dyn_used;  // dynamic memory actually used;
  int i, j;
  int ni;
  int nij_local; // (local to this function; the variable nij is static)

  ///////////////////////////////////////
  // the largest memory requirement will
  // either occur just before the end of
  // the 1. q.b.t. (mem1) or just before
  // the end of the i-batch loop (mem2)
  ///////////////////////////////////////

  mem_dyn = mem_alloc - mem_static;

  // First determine if calculation is possible at all (i.e., if ni=1 possible)

  ni = 1;
  // compute nij as nij on node 0, since nij on node 0 is >= nij on other nodes
  index = 0;
  nij_local = 0;
  for (i=0; i<ni; i++) {
    for (j=0; j<nocc; j++) {
      if (index++ % nproc == 0) nij_local++;
      }
    }
  mem1 = sizeof(double)*(nij_local*nbasis*nbasis + nbasis*nvir);
  mem2 = sizeof(double)*(ni*nbasis*nfuncmax*nfuncmax + nij_local*nbasis*nbasis
                         + ni*nbasis + nbasis*nfuncmax
                         + nfuncmax*nfuncmax*nfuncmax*nfuncmax);
  mem3 = sizeof(double)*(ni*nbasis*nfuncmax*nfuncmax + nij_local*nbasis*nbasis
                         + 2*(2 + nbasis*nfuncmax));
  tmp = (mem2>mem3 ? mem2:mem3);
  maxdyn = (mem1>tmp ? mem1:tmp);

  if (maxdyn > mem_dyn) {
    fprintf(stdout,"At least %i bytes required (for batch size 1)\n"
                   "but only %i bytes allocated; program exits\n",
                   maxdyn+mem_static, mem_alloc);
    abort();
    }

  ni = 2;
  dyn_used = maxdyn;
  while (ni<=nocc) {
    // compute nij
    nij_local = 0;
    index = 0;
    for (i=0; i<ni; i++) {
      for (j=0; j<nocc; j++) {
        if (index++ % nproc == 0) nij_local++;
        }
      }
    // compute max memory required for this ni
    mem1 = sizeof(double)*(nij_local*nbasis*nbasis + nbasis*nvir);
    mem2 = sizeof(double)*(ni*nbasis*nfuncmax*nfuncmax + nij_local*nbasis*nbasis
                           + ni*nbasis + nbasis*nfuncmax
                           + nfuncmax*nfuncmax*nfuncmax*nfuncmax);
    mem3 = sizeof(double)*(ni*nbasis*nfuncmax*nfuncmax + nij_local*nbasis*nbasis
                           + 2*(2 + nbasis*nfuncmax));
    tmp = (mem2>mem3 ? mem2:mem3);
    maxdyn = (mem1>tmp? mem1:tmp);
    if (maxdyn >= mem_dyn) {
      ni--;
      break;
      }
    else dyn_used = maxdyn;
    ni++;
    }
  if (ni > nocc) ni = nocc;

  fprintf(outfile,"Memory available per node:   %i Bytes\n",mem_alloc);
  fprintf(outfile,"Static memory used per node: %i Bytes\n",mem_static);
  fprintf(outfile,"Total memory used per node:  %i Bytes\n",dyn_used+mem_static);
  fprintf(outfile,"Batch size:                  %i\n", ni);

  return ni;
}

void
mbpt_mp2_gradient(scf_struct_t &scf_info,
                  sym_struct_t &sym_info,
                  centers_t &centers,
                  int nfzc, int nfzv,
                  dmt_matrix Scf_Vec, dmt_matrix Fock, dmt_matrix FockO,
                  int mem_alloc,
                  FILE *outfile,
                  RefMessageGrp &grp, double_matrix_t &gradient)
{
  dmt_matrix S = dmt_create("libscfv3 overlap matrix",scf_info.nbfao,SCATTERED);
  dmt_matrix SAHALF;
  dmt_matrix SC;
  dmt_matrix EVECS;
  dmt_matrix SCR;
  double_vector_t occ_num;
  double_vector_t evals;
  dmt_matrix SCR1,SCR2,SCR3;
  // this got free'ed somewhere
  SCR1 = dmt_create("opt2:scr1",scf_info.nbfao,COLUMNS);
  SCR2 = dmt_create("opt2:scr2",scf_info.nbfao,COLUMNS);
  SCR3 = dmt_create("opt2:scr3",scf_info.nbfao,COLUMNS);

  tim_enter("mp2grad");

  mbpt_ffo(S, &scf_info, &sym_info, &centers, Scf_Vec, Fock, FockO);

  /* form 'sahalf' matrix sahalf = u*ei^-0.5*u~ */
  allocbn_double_vector(&evals,"n",scf_info.nbfao);
  SAHALF= dmt_create("libscfv3 scf_core_guess scr4",scf_info.nbfao,COLUMNS);
  EVECS = dmt_create("libscfv3 scf_core_guess scr3",scf_info.nbfao,COLUMNS);
  SCR = dmt_create("libscfv3 scf_core_guess scr5",scf_info.nbfao,COLUMNS);
  SC = dmt_columns("libscfv3 scf_core_guess scr1",S);
  /* diagonalize overlap matrix */
  dmt_diag(SC,EVECS,evals.d);
  /* form SAHALF matrix (s^(-1/2), Sz&Ostl p. 143) */
  for(int i=0; i < scf_info.nbfao; i++) evals.d[i] = 1.0/sqrt(evals.d[i]);
  dmt_fill(SAHALF,0.0);
  dmt_set_diagonal(SAHALF,evals.d);
  /* form the orthogonalization matrix S^(-1/2) (Szabo&Ostlund p. 143)
   * (called SAHALF here) */
  dmt_transpose(EVECS);
  dmt_mult(SAHALF,EVECS,SCR);
  dmt_mult(EVECS,SCR,SAHALF);
  dmt_free(EVECS);
  dmt_free(SC);
  dmt_free(SCR);
  dmt_free(S);
  dmt_copy(Fock,SCR1); /* need a columns distr. matrix */
  dmt_mult(SCR1,SAHALF,SCR2);
  dmt_mult(SAHALF,SCR2,SCR3);
  /* SCR3 is now the Fock matrix in the orthogonalized ao basis */
  dmt_diag(SCR3,Scf_Vec,evals.d);
  dmt_copy(Scf_Vec,SCR1);
  dmt_mult(SAHALF,SCR1,Scf_Vec); /* Sz&Ostl p.146 point 9 */
  dmt_free(SCR1);
  dmt_free(SCR2);
  dmt_free(SCR3);

  sync0();

  tim_enter("mp2grad");
  mp2grad(&centers,&scf_info,Scf_Vec,&evals,nfzc,nfzv,mem_alloc,outfile,
          grp,&gradient);
  tim_exit("mp2grad");
       
  free_double_vector(&evals);

  tim_exit("mp2grad");
}
