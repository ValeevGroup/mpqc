
typedef int dmt_matrix;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <tmpl.h>

extern "C" {
#include <comm/picl/picl.h>
#include <comm/picl/ext/piclext.h>
}

#include <util/class/class.h>
#include <util/state/state.h>
#include <util/keyval/keyval.h>
#include <chemistry/molecule/molecule.h>

extern "C" {
#include <chemistry/qc/dmtqc/libdmtqc.h>
#include <util/misc/libmisc.h>
}

#include "mpqc_int.h"

/////////////////////////////////////////////////////////////////////

int
do_mp2(centers_t& centers, scf_struct_t& scf_info,
        dmt_matrix Scf_Vec, dmt_matrix Fock, FILE* outfile, RefKeyVal keyval)
{
  tim_enter("init");

  int nfzc=0,nfzv=0;
#if defined(PARAGON)
  int bmem=20000000;
#else
  int bmem=10000000;
#endif

  int me=mynode0();
  int nproc=numnodes0();
  int dim=cubedim0();

  if (me==0) {
    if (keyval->exists("base_memory")) {
      bmem = keyval->intvalue("frozen_docc");
      }

    if (keyval->exists("frozen_docc")) {
      nfzc = keyval->intvalue("frozen_docc");
      }

    if (keyval->exists("frozen_uocc")) {
      nfzv = keyval->intvalue("frozen_uocc");
      }
    }

  bcast0(&bmem,sizeof(int),mtype_get(),0);
  bcast0(&nfzc,sizeof(int),mtype_get(),0);
  bcast0(&nfzv,sizeof(int),mtype_get(),0);

  int_initialize_offsets2(&centers,&centers,&centers,&centers);

  int flags = INT_EREP|INT_NOSTRB|INT_NOSTR1|INT_NOSTR2;

  double *intbuf = 
    int_initialize_erep(flags,0,&centers,&centers,&centers,&centers);

  scf_init_bounds(&centers,intbuf);

  int a,b;
  int i,j,k,l;
  int p,q,r,s;
  int P,Q,R,S;
  int bf1,bf2,bf3,bf4;
  int nbasis=centers.nfunc;
  int nfuncmax = int_find_nfuncmax(&centers);
  int nocc=0,nvir;

  if (me==0) {
    if (keyval->exists("docc")) {
      for (i=0; i < keyval->count("docc"); i++)
        nocc += keyval->intvalue("docc",i);
    } else {
      for (i=0; i < centers.n; i++) nocc += (int) centers.center[i].charge;
      nocc = (nocc%2) ? nocc/2 + 1 : nocc/2 ;
    }
  }

  bcast0(&nocc,sizeof(int),mtype_get(),0);
  nvir = nbasis-nocc;

  nvir -= nfzv;
  nocc -= nfzc;

  /*
   * otay, let's get a grip here, and form the MO integrals
   * this is all in core for now till I figure out what the hell
   * I'm doing
   *
   * First thing, let's get the scf vector off the nodes
   */

  DMatrix scf_vector(nbasis,nocc+nvir);  // only doc + uoc
  double *evals = new double[nbasis];

  dmt_get_diagonal(Fock,evals);

  for (i=0; i < nocc+nvir; i++) evals[i] = evals[i+nfzc];

  loop_t *loop = dmt_ngl_create("%mr",Scf_Vec);
  while(dmt_ngl_next(loop)) {
    int iind,isize,jsize;
    double *col;

    dmt_ngl_create_inner(loop,0);
    while(dmt_ngl_next_inner_m(loop,&iind,&isize,&k,&jsize,&col)) {

      if (k >= nfzc && k < nbasis-nfzv) {
        for (i=0; i < nbasis; i++)
          scf_vector(i,k-nfzc) = col[i];
        }
      }
    }

  dmt_ngl_kill(loop);

 /*
  * ok, now we need to figure out how many passes we need to make.
  * for each i we need N(mM+36) m=nocc, M=nvir, N=nbasis
  * plus we need M^2+M all the time, plus NM for the scf vector
  *
  * so  nstatic = N(M+m+1) + M^2 *sizeof(double)
  *       nperi = N(mM/nproc + 2*nfuncmax^2)*sizeof(double)
  *       npass = nperi * m / (mem - nstatic)
  *       ni = m / npass
  */

  // mem use is N + N^2 + V^2 + ni*(2*nmx^2*N + OVN/nproc)

  // let's first allocate the rest of the bits that don't depend on ni
  double *ijab = new double[nvir*nvir];
  int *a2_lengths = new int[nproc];
  int *a2_lens = new int[nproc];

  int nstatic = (nbasis*(nvir+nocc+1)+nvir*nvir)*sizeof(double);
  int mem = bmem-nstatic;

  int nj = nocc*nvir*nbasis/nproc;
  int nperi = (2*nfuncmax*nfuncmax*nbasis+nj);

  int npass = nperi*nocc*sizeof(double)/mem;
  if (nperi*nocc%mem) npass++;

  int ni = nocc/npass;
  if (nocc%npass) ni++;

  int nij = ni*nocc/nproc;
  if (ni*nocc%nproc) nij++;

  if (me==0) {
    fprintf(outfile,"\n  in mp2:\n");
    fprintf(outfile,  "     nstatic = %d\n",nstatic);
    fprintf(outfile,  "         mem = %d\n",mem);
    fprintf(outfile,  "       npass = %d\n",npass);
    fprintf(outfile,  "          ni = %d\n",ni);
    fprintf(outfile,  "         nij = %d\n",nij);
    fprintf(outfile,  "      nshell = %d\n",centers.nshell);
    fprintf(outfile,  "        nfzc = %d\n",nfzc);
    fprintf(outfile,  "        nfzv = %d\n",nfzv);
    fprintf(outfile,  "        nocc = %d\n",nocc);
    fprintf(outfile,  "        nvir = %d\n",nvir);
    }

 // first we must determine how big each nodes chunk will be

  int mylena = nvir/nproc;
  int remaina = nvir%nproc;
  if (me < remaina) mylena++;

  int firsta;
  if (me < remaina) firsta = me*mylena;
  else firsta = remaina*(mylena+1) + (me-remaina)*mylena;

 // fill the a2_lengths vector for gcollect
  for (int ii=0; ii < nproc; ii++)
    a2_lengths[ii] = (nvir/nproc)*sizeof(double);

  for (ii=0; ii < remaina; ii++)
    a2_lengths[ii] += sizeof(double);

  //for (ii=0; ii < nproc; ii++) 
    //a2_lengths[ii] *= ni;


  double *T1 = new double[nfuncmax*nfuncmax*nbasis*ni];
  double *T2 = new double[nfuncmax*nfuncmax*nbasis*ni];
  double *T3 = new double[nij*nbasis*nvir];

  tim_exit("init");

  double Ecorr=0;

  int pass;
  for (pass=0; pass < npass; pass++) {

    int remaini = nocc%npass;

    if (remaini && remaini==pass) ni--;

    int firsti;

    if (pass < remaini) firsti = pass*ni;
    else firsti = remaini*(ni+1) + (pass-remaini)*ni;

    memcpy(T3,0,nij*nbasis*nvir*sizeof(double));

    int tol = (int) (-10.0/log10(2.0));

    tim_enter("RS loop");
    int int_index=0;
    for (R=0; R < centers.nshell; R++) {

if(me==0) printf(" %5d\n",R);
      int nr = INT_SH_NFUNC((&centers),R);

      for (S=0; S <= R; S++) {
        int ns = INT_SH_NFUNC((&centers),S);

        memcpy(T1,0,nr*ns*nbasis*ni*sizeof(double));

        tim_enter("PQ loop");
        for (P=0; P < centers.nshell; P++) {
          int np = INT_SH_NFUNC((&centers),P);

          for (Q=0; Q <= P; Q++) {
            if (scf_erep_bound(P,Q,R,S) < tol){
              continue;
              }

            int_index++;

            if (int_index%nproc != me) continue;

            int nq = INT_SH_NFUNC((&centers),Q);

            tim_enter("erep");
            int_erep(INT_EREP|INT_NOBCHK|INT_NOPERM|INT_REDUND,&P,&Q,&R,&S);
            tim_exit("erep");
        
            int index=0;

            tim_enter("a1 loop");
            for (bf1=0; bf1 < np; bf1++) {
              p = centers.func_num[P] + bf1;

              for (bf2=0; bf2 < nq ; bf2++) {
                q = centers.func_num[Q] + bf2;

                if (q > p) {
                  index += nr*ns;
                  continue;
                  }

                for (bf3=0; bf3 < nr ; bf3++) {
                  for (bf4=0; bf4 < ns ;bf4++,index++) {
                    if (R==S && bf4>bf3) continue;

                    if (INT_NONZERO(intbuf[index])) {
                      double pqrs = intbuf[index];

                      if (p==q) pqrs *= 0.5;

                      double *rsqi = &T1[((bf3*ns+bf4)*nbasis+q)*ni];
                      double *rspi = &T1[((bf3*ns+bf4)*nbasis+p)*ni];
                      double *c_p = &scf_vector[p][firsti];
                      double *c_q = &scf_vector[q][firsti];

                      for (i=0; i < ni; i++) {
                        *rsqi++ += pqrs * *c_p++;
                        *rspi++ += pqrs * *c_q++;
                        }
                      }
                    }
                  }
                }
              }
            tim_exit("a1 loop");
            }
          }
        tim_exit("PQ loop");
        
        tim_enter("a1 sum");
        gop1(T1,nr*ns*nbasis*ni,T2,'+',3);
        tim_exit("a1 sum");

        memcpy(T2,0,ni*nr*ns*nvir*sizeof(double));

        tim_enter("a2 loop");
        for (a=0; a < mylena; a++) {
          for (q=0; q < nbasis ; q++) {
            double c_qa = scf_vector[q][a+firsta+nocc];

            for (bf3=0; bf3 < nr ; bf3++) {
              for (bf4=0; bf4 <= ((R==S) ? bf3 : ns-1) ; bf4++) {
                double *a2ars = &T2[((a*nr+bf3)*ns+bf4)*ni];
                double *a1rsq = &T1[((bf3*ns+bf4)*nbasis+q)*ni];

                for (i=ni; i; i--) *a2ars++ += c_qa * *a1rsq++;
                }
              }
            }
          }
        tim_exit("a2 loop");

        tim_enter("junk");

     // fill the a2_lens vector for gcollect

        for (ii=0; ii < nproc; ii++) 
          a2_lens[ii] = a2_lengths[ii]*nr*ns*ni;

        tim_enter("a2 sum");
        gcollect(T2,a2_lens,T1);
        tim_exit("a2 sum");

        tim_enter("flip");
        double *a2r = T2;
        int rsi=0;
        int nrsi=nr*ns*ni;

        for (i=0; i < ni; i++)
          for (r=0; r < nr; r++)
            for (s=0; s < ns; s++,rsi++) {
              double *a2rsi = &T1[r*ns*ni+s*ni+i];
              for (a=0; a < nvir; a++,a2rsi+=nrsi)
                *a2r++ = *a2rsi;
              }
        tim_exit("flip");
        tim_exit("junk");

        tim_enter("a3 loop");
        for (bf3=0; bf3 < nr ; bf3++) {
          r = centers.func_num[R] + bf3;

          for (bf4=0; bf4 <= ((R==S) ? bf3 : ns-1) ; bf4++) {
            s = centers.func_num[S] + bf4;

            double scal = (r==s) ? 0.5 : 1;

            int ij=0;
            for (i=0; i < ni; i++) {

              double *iars = &T2[i*nr*ns*nvir+bf3*ns*nvir+bf4*nvir];

              for (j=0; j < nocc; j++) {
                if (((i*nocc+j)%nproc) != me) continue;
                double c_rj = scf_vector(r,j)*scal;
                double c_sj = scf_vector(s,j)*scal;

                double *a3ijs = &T3[ij*nbasis*nvir+s*nvir];
                double *a3ijr = &T3[ij*nbasis*nvir+r*nvir];

                for (a=0; a < nvir; a++) {
                  a3ijs[a] += iars[a]*c_rj;
                  a3ijr[a] += iars[a]*c_sj;
                  }

                ij++;
                }
              }
            }
          }

        tim_exit("a3 loop");
        }
      }
    tim_exit("RS loop");

    tim_enter("a4 loop");
    int ij=0;
    for (i=0; i < ni; i++) {
      for (j=0; j < nocc; j++) {
        if (((i*nocc+j)%nproc) != me) continue;

        double eij = -evals[i+firsti]-evals[j];

        memcpy(ijab,0,sizeof(double)*nvir*nvir);

        for (s=0; s < nbasis; s++) {
          for (a=0; a < nvir; a++) {
            double ijsa = T3[ij*nbasis*nvir+s*nvir+a];
            double *c_s = &scf_vector[s][nocc];
            double *ijabp = & ijab[a*nvir];

            for (b=0; b < nvir; b++) {
              *ijabp++ += ijsa * *c_s++;
              }
            }
          }
        
        for (a=0; a < nvir; a++) {
          double ea=evals[a+nocc];
          for (b=0; b < nvir; b++) {
            double eb=evals[b+nocc];
            double t1 = ijab[a*nvir+b]*ijab[a*nvir+b];
            double t2 = ijab[a*nvir+b]*ijab[b*nvir+a];

            Ecorr += (2*t1-t2)/(ea+eb+eij);
            }
          }

        ij++;
        }
      }
    tim_exit("a4 loop");
    }

  gsum0(&Ecorr,1,5,mtype_get(),0);

  if (me==0) {
    printf("\n");
    printf("  E(SCF)  = %20.10f\n",scf_info->e_elec+scf_info->nuc_rep);
    printf("  E(Corr) = %20.10f\n",-Ecorr);
    printf("  E(MP2)  = %20.10f\n\n",scf_info->e_elec+scf_info->nuc_rep-Ecorr);
    }

  scf_done_bounds();
  int_done_erep();
  int_done_offsets2(&centers,&centers,&centers,&centers);

  delete[] evals;
  delete[] ijab;
  delete[] T1;
  delete[] T2;
  delete[] a2_lengths;
  delete[] a2_lens;

  return(0);
  }
