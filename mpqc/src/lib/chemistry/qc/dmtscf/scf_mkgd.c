
#define SHELLS 125
#define HSOS 1

/* Calculates two-electron integrals on the fly and sticks them into the
 * appropriate part of the G matrix
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <tmpl.h>
#include <util/group/picl.h>
#include <math/array/math_lib.h>
#include <util/misc/libmisc.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <chemistry/qc/dmtsym/sym_dmt.h>

#include <chemistry/qc/dmtscf/scf.h>
#include <chemistry/qc/dmtscf/scf_bnd.gbl>
#include <chemistry/qc/dmtscf/scf_gmat.gbl>

#define ioff(i) (((i)*((i)+1))>>1)
#define IOFF(a,b) (((a)>(b))?(ioff(a)+(b)):(ioff(b)+(a)))

extern signed char *scf_bnd_Qvec;

#include <chemistry/qc/dmtscf/scf_mkgd.gbl>
#include <chemistry/qc/dmtscf/scf_mkgd.lcl>

/**************************************************************************
 *
 */

GLOBAL_FUNCTION int
scf_make_g_d(centers,scf_info,sym_info,Gmat,GmatO,DPmat,DPmatO,mgdbuff,outfile)
centers_t *centers;
scf_struct_t *scf_info;
sym_struct_t *sym_info;
dmt_matrix Gmat;
dmt_matrix GmatO;
dmt_matrix DPmat;
dmt_matrix DPmatO;
double *mgdbuff;
FILE *outfile;
{
  int tmax,imax,cpmax,pmaxijk;
  int pmaxik,pmaxjk,pmaxij,Qvecij;
  int i,j,k,l;
  int ij,kl,ijkl;
  int g,gij,gkl,gijkl;
  int n1,n2,n3,n4;
  int e12,e34,e13e24,e_any;
  int bf1,bf2,bf3,bf4;
  int i1,j1,k1,l1;
  int i2,j2,k2,l2;
  int ii,jj,kk,ll;
  int ij1;
  int lij,lkl;
  int index;
  int nijkl,leavel,qijkl=1;
  int int_index,kindex;
  int nproc=numnodes0();
  int me=mynode0();
  int s1,s2,s3,s4;
  int errcod;
  int inttol = (int) ((double) -(scf_info->intcut)/log10(2.0));
  int use_symmetry = (sym_info->g > 1) ? 1 : 0;

  double tnint=0.0;
  double pki_int,value;
  double *gtmp=0, *ptmp=0, *ptmpo=0, *gtmpo=0;

  char *shnfunc;
  signed char *maxp=0;

#if !HSOS
  if (scf_info->iopen) {
    fprintf(stderr,"oops, this version cannot do open shells\n");
    return -1;
  }
#endif

 /* start timing */
  tim_enter("scf_mkgd");

 /* first let's transfer the density matrices to vectors and init gtmp's */
  errcod =
    scf_init_direct_gmat(scf_info,DPmat,DPmatO,&ptmp,&ptmpo,&maxp,&gtmp,&gtmpo);

  if (errcod < 0)  {
    fprintf(stderr,"scf_make_g_d:  scf_init_direct_gmat() failed\n");
    return -1;
  }

 /* form shnfunc array */

  shnfunc = (char *) malloc(centers->nshell);
  if (shnfunc==NULL) {
    fprintf(stderr,"scf_make_g_d: could not malloc shnfunc\n");
    scf_done_direct_gmat(ptmp,ptmpo,maxp,gtmp,gtmpo);
    return -1;
  }

  for (i=0; i < centers->nshell; i++) shnfunc[i]=INT_SH_NFUNC((centers),i);

 /* loop over all shells, calculate a bunch of integrals from each shell
  * quartet, and stick those integrals where they belong
  */

  kindex=int_index=0;
  for (i=0; i<centers->nshell; i++) {
    if (use_symmetry && !sym_info->p1[i])
      continue;

    for (j=0; j<=i; j++) {
      leavel=0;
      ij = ioff(i)+j;

      if (use_symmetry && !sym_info->lamij[ij])
        continue;

      Qvecij=(int)scf_bnd_Qvec[ij];
      if (scf_info->eliminate)
        pmaxij=maxp[ij];

      for (k=0; k<=i; k++,kindex++) {
        if (kindex%nproc!=me)
          continue;

        kl=ioff(k);
        if (scf_info->eliminate) {
          pmaxijk=pmaxij;
          if ((pmaxik=maxp[(ioff(i)+k)]-2)>pmaxijk)
            pmaxijk=pmaxik;
          if ((pmaxjk=maxp[IOFF(j,k)]-2)>pmaxijk)
            pmaxijk=pmaxjk;
        }

        /*tim_enter("l loop");*/
        for (l=0; l<=(k==i?j:k); l++) {
          if (use_symmetry) {
            ijkl=ioff(ij)+kl;
            nijkl=leavel=0;
            for (g=0; g < sym_info->g ; g++) {
              gij=IOFF(sym_info->shell_map[i][g],sym_info->shell_map[j][g]);
              gkl=IOFF(sym_info->shell_map[k][g],sym_info->shell_map[l][g]);
              gijkl = IOFF(gij,gkl);
              if (gijkl > ijkl) {
                leavel=1;
                break;
              }
              if (gijkl == ijkl)
                nijkl++;
            }
            if (leavel) {
              kl++;
              continue;
            }
            qijkl = sym_info->g/nijkl;
          }

          imax = (int) scf_bnd_Qvec[kl]+Qvecij;

          if (scf_info->eliminate) {
            cpmax = (maxp[kl]>pmaxijk) ? maxp[kl] : pmaxijk;
            if ((tmax=maxp[(ioff(i)+l)]-2)>cpmax)
              cpmax=tmax;
            if ((tmax=maxp[IOFF(j,l)]-2)>cpmax)
              cpmax=tmax;

            if (cpmax+imax < inttol) {
              /* If we are trying to save integrals on this node, then
               * int_index must be incremented now. */
              if (scf_info->int_store) int_index++;
              kl++;
              continue;
            }
          }

          /*tim_enter("ints");*/
          s1 = i; s2 = j; s3 = k; s4 = l;

          tim_enter("gd_erep");
          int_erep(INT_EREP|INT_NOBCHK|INT_NOPERM,&s1,&s2,&s3,&s4);
          tim_exit("gd_erep");

          n1 = shnfunc[s1];
          n2 = shnfunc[s2];
          n3 = shnfunc[s3];
          n4 = shnfunc[s4];

         /* Shell equivalency information. */
          e12    = (s2==s1);
          e13e24 = (s3==s1) && (s4==s2);
          e34    = (s4==s3);

          index = 0;

          e_any = (e12||e13e24||e34);
          if (e_any) {
            for (bf1=0; bf1<=INT_MAX1(n1) ; bf1++) {
              i1 = centers->func_num[s1] + bf1;

              for (bf2=0; bf2<=INT_MAX2(e12,bf1,n2) ; bf2++) {
                j1 = centers->func_num[s2] + bf2;
                ij1=ioff(i1)+j1;

                for (bf3=0; bf3<=INT_MAX3(e13e24,bf1,n3) ; bf3++) {
                  k1 = centers->func_num[s3] + bf3;

                  for (bf4=0;bf4<=INT_MAX4(e13e24,e34,bf1,bf2,bf3,n4);bf4++){
                    if (INT_NONZERO(mgdbuff[index])) {
                      l1 = centers->func_num[s4] + bf4;

                      if (ij1 >= ioff(k1)+l1) {
                        ii = i1; jj = j1; kk = k1; ll = l1;
                      } else {
                        ii = k1; jj = l1; kk = i1; ll = j1;
                      }

                      pki_int = mgdbuff[index];
                      if (qijkl > 1)
                        pki_int *= qijkl;

                      if (jj == kk) {
                        /*
                         * if i=j=k or j=k=l, then this integral contributes
                         * to J, K1, and K2 of G(ij), so
                         * value = (ijkl) - 0.25 * ((ikjl)-(ilkj))
                         *       = 0.5 * (ijkl)
                         */
                        if (ii == jj || kk == ll) {
                          lij=ioff(ii)+jj;
                          lkl=ioff(kk)+ll;
                          value = (lij==lkl) ? 0.25*pki_int: 0.5*pki_int;
                          gtmp[lij] += ptmp[lkl]*value;
                          gtmp[lkl] += ptmp[lij]*value;
#if HSOS
                          if (scf_info->hsos) {
                            gtmpo[lij] += ptmpo[lkl]*value;
                            gtmpo[lkl] += ptmpo[lij]*value;
                          }
#endif
                        } else {
                          /*
                           * if j=k, then this integral contributes
                           * to J and K1 of G(ij)
                           *
                           * value = (ijkl) - 0.25 * (ikjl)
                           *       = 0.75 * (ijkl)
                           */
                          lij=ioff(ii)+jj;
                          lkl=ioff(kk)+ll;
                          value = (lij==lkl) ? 0.375*pki_int: 0.75*pki_int;

                          gtmp[lij] += ptmp[lkl]*value;
                          gtmp[lkl] += ptmp[lij]*value;
#if HSOS
                          if (scf_info->hsos) {
                            value *= 0.3333333333333333;
                            gtmpo[lij] += ptmpo[lkl]*value;
                            gtmpo[lkl] += ptmpo[lij]*value;
                          }
#endif

                          /*
                           * this integral also contributes to K1 and K2 of
                           * G(il)
                           *
                           * value = -0.25 * ((ijkl)+(ikjl))
                           *       = -0.5 * (ijkl)
                           */
                          lij=ioff(ii)+ll;
                          lkl=IOFF(kk,jj);
                          value = (lij==lkl)? 0.25*pki_int: 0.5*pki_int;

                          gtmp[lij] -= ptmp[lkl]*value;
                          gtmp[lkl] -= ptmp[lij]*value;
#if HSOS
                          if (scf_info->hsos) {
                            gtmpo[lij] += ptmpo[lkl]*value;
                            gtmpo[lkl] += ptmpo[lij]*value;
                          }
#endif
                        }
                      } else if (ii == kk || jj == ll) {
                        /*
                         * if i=k or j=l, then this integral contributes
                         * to J and K2 of G(ij)
                         *
                         * value = (ijkl) - 0.25 * (ilkj)
                         *       = 0.75 * (ijkl)
                         */
                        lij=ioff(ii)+jj;
                        lkl=ioff(kk)+ll;
                        value = (lij==lkl) ? 0.375*pki_int: 0.75*pki_int;
                        gtmp[lij] += ptmp[lkl]*value;
                        gtmp[lkl] += ptmp[lij]*value;
#if HSOS
                        if (scf_info->hsos) {
                          value *= 0.3333333333333333;
                          gtmpo[lij] += ptmpo[lkl]*value;
                          gtmpo[lkl] += ptmpo[lij]*value;
                        }
#endif

                        /*
                         * this integral also contributes to K1 and K2 of
                         * G(ik)
                         *
                         * value = -0.25 * ((ijkl)+(ilkj))
                         *       = -0.5 * (ijkl)
                         */
                        lij=ioff(ii)+kk;
                        lkl=IOFF(jj,ll);
                        value = (lij==lkl) ? 0.25*pki_int : 0.5*pki_int;

                        gtmp[lij] -= ptmp[lkl]*value;
                        gtmp[lkl] -= ptmp[lij]*value;
#if HSOS
                        if (scf_info->hsos) {
                          gtmpo[lij] += ptmpo[lkl]*value;
                          gtmpo[lkl] += ptmpo[lij]*value;
                        }
#endif
                      } else {
                        /*
                         * This integral contributes to J of G(ij)
                         *
                         * value = (ijkl)
                         */
                        lij=ioff(ii)+jj;
                        lkl=ioff(kk)+ll;
                        value = (lij==lkl)? 0.5*pki_int : pki_int;

                        gtmp[lij] += ptmp[lkl]*value;
                        gtmp[lkl] += ptmp[lij]*value;

                        /*
                         * and to K1 of G(ik)
                         *
                         * value = -0.25 * (ijkl)
                         */
                        lij=ioff(ii)+kk;
                        lkl=IOFF(jj,ll);
                        value = (lij==lkl) ? 0.125*pki_int : 0.25*pki_int;

                        gtmp[lij] -= ptmp[lkl]*value;
                        gtmp[lkl] -= ptmp[lij]*value;
#if HSOS
                        if (scf_info->hsos) {
                          gtmpo[lij] += ptmpo[lkl]*value;
                          gtmpo[lkl] += ptmpo[lij]*value;
                        }
#endif

                        if ((ii != jj) && (kk != ll)) {
                          /*
                           * if i!=j and k!=l, then this integral also
                           * contributes to K2 of G(il)
                           *
                           * value = -0.25 * (ijkl)
                           *
                           * note: if we get here, then ik can't equal jl,
                           * so value wasn't multiplied by 0.5 above.
                           */
                          lij=ioff(ii)+ll;
                          lkl=IOFF(kk,jj);

                          gtmp[lij] -= ptmp[lkl]*value;
                          gtmp[lkl] -= ptmp[lij]*value;
#if HSOS
                          if (scf_info->hsos) {
                            gtmpo[lij] += ptmpo[lkl]*value;
                            gtmpo[lkl] += ptmpo[lij]*value;
                          }
#endif
                        }
                      }
                    }
                    index++;
                  }
                }
              }
            }
          } else {
            for (i1=centers->func_num[s1], bf1=0; bf1<n1 ; bf1++, i1++) {

              for (j1=centers->func_num[s2], bf2=0; bf2<n2 ; bf2++, j1++) {
                ij1=ioff(i1)+j1;

                for (k1=centers->func_num[s3], bf3=0; bf3<n3 ; bf3++, k1++) {

                  for (l1=centers->func_num[s4], bf4=0; bf4<n4 ; bf4++, l1++) {
                    if (INT_NONZERO(mgdbuff[index])) {

                      if (ij1 >= ioff(k1)+l1) {
                        ii = i1; jj = j1; kk = k1; ll = l1;
                      } else {
                        ii = k1; jj = l1; kk = i1; ll = j1;
                      }

                      pki_int = mgdbuff[index];
                      if (qijkl > 1)
                        pki_int *= qijkl;

                      if (jj == kk) {
                        /*
                         * if j=k, then this integral contributes
                         * to J and K1 of G(ij)
                         *
                         * value = (ijkl) - 0.25 * (ikjl)
                         *       = 0.75 * (ijkl)
                         */
                        lij=ioff(ii)+jj;
                        lkl=ioff(kk)+ll;
                        value = 0.75*pki_int;

                        gtmp[lij] += ptmp[lkl]*value;
                        gtmp[lkl] += ptmp[lij]*value;
#if HSOS
                        if (scf_info->hsos) {
                          value *= 0.3333333333333333;
                          gtmpo[lij] += ptmpo[lkl]*value;
                          gtmpo[lkl] += ptmpo[lij]*value;
                        }
#endif

                        /*
                         * this integral also contributes to K1 and K2 of
                         * G(il)
                         *
                         * value = -0.25 * ((ijkl)+(ikjl))
                         *       = -0.5 * (ijkl)
                         */
                        lij=ioff(ii)+ll;
                        lkl=IOFF(kk,jj);
#if HSOS
                        value = 0.5*pki_int;
#else
                        value *= 0.666666666666666;
#endif
                        gtmp[lij] -= ptmp[lkl]*value;
                        gtmp[lkl] -= ptmp[lij]*value;
#if HSOS
                        if (scf_info->hsos) {
                          gtmpo[lij] += ptmpo[lkl]*value;
                          gtmpo[lkl] += ptmpo[lij]*value;
                        }
#endif
                      } else if (ii == kk || jj == ll) {
                        /*
                         * if i=k or j=l, then this integral contributes
                         * to J and K2 of G(ij)
                         *
                         * value = (ijkl) - 0.25 * (ilkj)
                         *       = 0.75 * (ijkl)
                         */
                        lij=ioff(ii)+jj;
                        lkl=ioff(kk)+ll;
                        value = 0.75*pki_int;

                        gtmp[lij] += ptmp[lkl]*value;
                        gtmp[lkl] += ptmp[lij]*value;
#if HSOS
                        if (scf_info->hsos) {
                          value *= 0.3333333333333333;
                          gtmpo[lij] += ptmpo[lkl]*value;
                          gtmpo[lkl] += ptmpo[lij]*value;
                        }
#endif

                        /*
                         * this integral also contributes to K1 and K2 of
                         * G(ik)
                         *
                         * value = -0.25 * ((ijkl)+(ilkj))
                         *       = -0.5 * (ijkl)
                         */
                        lij=ioff(ii)+kk;
                        lkl=IOFF(jj,ll);
#if HSOS
                        value = 0.5*pki_int;
#else
                        value *= 0.666666666666666;
#endif                        

                        gtmp[lij] -= ptmp[lkl]*value;
                        gtmp[lkl] -= ptmp[lij]*value;
#if HSOS
                        if (scf_info->hsos) {
                          gtmpo[lij] += ptmpo[lkl]*value;
                          gtmpo[lkl] += ptmpo[lij]*value;
                        }
#endif
                      } else {
                        /*
                         * This integral contributes to J of G(ij)
                         *
                         * value = (ijkl)
                         */
                        lij=ioff(ii)+jj;
                        lkl=ioff(kk)+ll;
                        value = pki_int;

                        gtmp[lij] += ptmp[lkl]*value;
                        gtmp[lkl] += ptmp[lij]*value;

                        /*
                         * and to K1 of G(ik)
                         *
                         * value = -0.25 * (ijkl)
                         */
                        lij=ioff(ii)+kk;
                        lkl=IOFF(jj,ll);
                        value *= 0.25;

                        gtmp[lij] -= ptmp[lkl]*value;
                        gtmp[lkl] -= ptmp[lij]*value;
#if HSOS
                        if (scf_info->hsos) {
                          gtmpo[lij] += ptmpo[lkl]*value;
                          gtmpo[lkl] += ptmpo[lij]*value;
                        }
#endif

                        /*
                         * and to K2 of G(il)
                         *
                         * value = -0.25 * (ijkl)
                         */
                        lij=ioff(ii)+ll;
                        lkl=IOFF(kk,jj);

                        gtmp[lij] -= ptmp[lkl]*value;
                        gtmp[lkl] -= ptmp[lij]*value;
#if HSOS
                        if (scf_info->hsos) {
                          gtmpo[lij] += ptmpo[lkl]*value;
                          gtmpo[lkl] += ptmpo[lij]*value;
                        }
#endif
                      }
                    }
                    index++;
                  }
                }
              }
            }
          }
          tnint += (double) (n1*n2*n3*n4);
          /*tim_exit("ints");*/

          kl++;
          int_index++;
        }
        /*tim_exit("l loop");*/
      }
    }
  }

  if (scf_info->print_flg & 4) {
    gsum0(&tnint,1,5,mtype_get(),0);
    if (mynode0()==0) {
      fprintf(outfile,"  %8.0f integrals in scf_make_g_d\n",tnint);
      fflush(outfile);
    }
  }

  free(shnfunc);

 /* now sum up contributions to gtmp and gtmpo */
  gop1(gtmp,scf_info->nbatri,ptmp,'+',mtype_get());
  if (scf_info->iopen) gop1(gtmpo,scf_info->nbatri,ptmp,'+',mtype_get());

 /* and stuff gtmp's back into Gmat's */
  scf_lgmat_to_scat(gtmp,Gmat);
  if (scf_info->iopen) scf_lgmat_to_scat(gtmpo,GmatO);

 /* and free up some memory */
  scf_done_direct_gmat(ptmp,ptmpo,maxp,gtmp,gtmpo);

  tim_exit("scf_mkgd");

  int_reduce_storage_threshold();

  return 0;
}
