
#define SHELLS 125
/* Calculates two-electron integrals on the fly and sticks them into the
 * appropriate part of the G matrix
 */

/* $Log$
 * Revision 1.1  1993/12/29 12:53:16  etseidl
 * Initial revision
 *
 * Revision 1.16  1992/06/18  22:48:18  seidl
 * change how work is partitioned
 *
 * Revision 1.15  1992/06/17  21:54:24  jannsen
 * cleaned up for saber-c
 *
 * Revision 1.14  1992/06/16  16:27:36  seidl
 * add open shell stuff
 *
 * Revision 1.13  1992/05/26  20:18:00  jannsen
 * use mtype_get to get message types for global operations
 * check results of memory allocations
 *
 * Revision 1.12  1992/05/19  21:05:12  seidl
 * add cheat stuff
 *
 * Revision 1.11  1992/05/12  10:34:45  seidl
 * use INT_NOPERM flag in int_erep() call
 *
 * Revision 1.10  1992/05/04  11:11:22  seidl
 * inline everything now, remove as much as possible from inner loops
 *
 * Revision 1.9  1992/04/09  17:54:42  seidl
 * use signed char for maxp array
 *
 * Revision 1.8  1992/04/08  20:38:13  seidl
 * use unsigned ints, make tnint double
 *
 * Revision 1.7  1992/04/07  18:08:40  jannsen
 * took out the pmax, findit, and filter tim_enter's
 *
 * Revision 1.6  1992/04/07  18:04:18  jannsen
 *
 * Revision 1.5  1992/04/06  12:37:24  seidl
 * change how number of integrals is calculated
 *
 * Revision 1.4  1992/04/01  01:04:16  seidl
 * fix bounds checking
 *
 * Revision 1.3  1992/03/31  22:26:29  seidl
 * use new bounds checking
 *
 * Revision 1.2  1992/03/21  00:40:25  seidl
 * change sym_libv2.h to chemistry/qc/dmtsym/sym_dmt.h
 * test boolean eliminate
 *
 * Revision 1.1.1.1  1992/03/17  16:26:40  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  16:26:38  seidl
 * Initial revision
 *
 * Revision 1.7  1992/03/09  12:58:15  seidl
 * make max_den more efficient
 *
 * Revision 1.6  1992/03/04  15:58:56  seidl
 * use scf_erep_bound, add some timing stuff
 *
 * Revision 1.5  1992/02/28  19:22:15  seidl
 * remove irrep type stuff (ism, ior)
 *
 * Revision 1.4  1992/02/26  12:20:00  seidl
 * split up work among nodes more efficiently
 *
 * Revision 1.3  1992/02/13  00:39:16  seidl
 * fixed datatype for gsum0. was 5(double), should be 2(int)
 *
 * Revision 1.2  1992/02/10  17:02:34  seidl
 * remove scalar code, remove most timing calls
 *
 * Revision 1.1  1992/02/04  23:48:08  seidl
 * Initial revision
 *
 * Revision 1.10  1992/01/16  19:48:37  seidl
 * use libintv2 now, qijkl is an int, not a double
 *
 * Revision 1.9  1992/01/13  19:16:16  seidl
 * ao density formed in scf_make_gmat now
 *
 * Revision 1.8  1992/01/10  16:18:55  seidl
 * add some timing stuff
 *
 * Revision 1.7  1992/01/09  11:59:56  seidl
 * more efficient parallel version
 *
 * Revision 1.6  1992/01/08  20:13:24  seidl
 * inefficient parallel version
 *
 * Revision 1.5  1992/01/06  11:48:57  seidl
 * attempt to make max_den a little more efficient
 *
 * Revision 1.4  1992/01/02  18:14:28  seidl
 * fix defs of inttol and save_thr, skip loop if imax==0
 *
 * Revision 1.3  1992/01/02  16:21:55  seidl
 * multiply diagonal elements of pmat by 2 before transforming to ao basis
 * set noi and nok = iopen rather than socc[i]
 *
 * Revision 1.2  1991/12/24  19:32:14  seidl
 * fix assignment of ptmpo
 *
 * Revision 1.1  1991/12/20  16:23:14  seidl
 * Initial revision
 *
 * Revision 1.4  1991/11/15  14:18:26  seidl
 * *** empty log message ***
 *
 * Revision 1.3  1991/10/23  22:11:53  seidl
 * move includes of local and global below declaration of struct mgd
 *
 * Revision 1.2  1991/09/23  00:43:55  seidl
 * *** empty log message ***
 *
 * Revision 1.1  1991/06/22  09:03:57  seidl
 * Initial revision
 * */

static char *rcsid = "$Id$";

#include <stdio.h>
#include <stdlib.h>
#include <tmpl.h>
#include <comm/picl/picl.h>
#include <math/array/math_lib.h>
#include <util/misc/libmisc.h>
#include <math.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <chemistry/qc/dmtsym/sym_dmt.h>
#include "scf.h"
#include "scf_bnd.gbl"

#ifdef INT_CE_SH_AM
#undef INT_CE_SH_AM
#endif
#define INT_CE_SH_AM(c,a,s) ((c)->center[(a)].basis.shell[(s)].type[0].am)

#define ioff(i) (((i)*((i)+1))>>1)
#define IOFF(a,b) (((a)>(b))?(ioff(a)+(b)):(ioff(b)+(a)))

extern signed char *Qvec;

#include "scf_mkgd.gbl"
#include "scf_mkgd.lcl"

GLOBAL_FUNCTION int
scf_make_g_d(_centers,_irreps,_scf_info,_sym_info,
                 _gmat,_gmato,_dpmat,_dpmato,maxp,_mgdbuff,iter,_outfile)
centers_t *_centers;
scf_irreps_t *_irreps;
scf_struct_t *_scf_info;
sym_struct_t *_sym_info;
double_vector_t *_gmat;
double_vector_t *_gmato;
double_vector_t *_dpmat;
double_vector_t *_dpmato;
signed char *maxp;
double *_mgdbuff;
int iter;
FILE *_outfile;
{
  register int tmax,imax,cpmax,pmaxijk;
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
  double tnint=0.0;
  double pki_int,value;
  double *gtmp,*ptmp,*ptmpo,*gtmpo;
  int inttol = (int) ((double) -(_scf_info->intcut)/log10(2.0));
  int use_symmetry=(_sym_info->g > 1)?1:0;
  char *shnfunc;

  tim_enter("scf_mkgd");

  if(_scf_info->cheat) {
    switch(iter) {
    case 0:
    case 1:
    case 2:
    case 3:
    /* ~10e-5 */
      if(inttol<-16) inttol = -16;
      break;
    default:
     /* use user input tolerance */
      break;
      }
    }

  gtmp = _gmat->d;
  ptmp = _dpmat->d;
  if(_scf_info->iopen) {
    gtmpo = _gmato->d;
    ptmpo = _dpmato->d;
    }

 /* form shnfunc array */

  shnfunc = (char *) malloc(_centers->nshell);
  if(shnfunc==NULL) {
    fprintf(_outfile,"could not malloc shnfunc\n");
    return(-1);
    }
  for (i=0; i<_centers->nshell; i++) shnfunc[i]=INT_SH_NFUNC((_centers),i);

 /* loop over all shells, calculate a bunch of integrals from each shell
  * quartet, and stick those integrals where they belong
  */

  kindex=int_index=0;
  for (i=0; i<_centers->nshell; i++) {
    if(use_symmetry) if(!_sym_info->p1[i]) continue;

    for (j=0; j<=i; j++) {
      leavel=0;
      ij = ioff(i)+j;
      if(use_symmetry) if(!_sym_info->lamij[ij]) continue;
      Qvecij=(int)Qvec[ij];
      if(_scf_info->eliminate) pmaxij=maxp[ij];

      for (k=0; k<=i; k++,kindex++) {
#if 0
        if((_centers->nshell > SHELLS || use_symmetry) && kindex%nproc!=me) {
          continue;
          }
#else
        if(kindex%nproc!=me) {
          continue;
          }
#endif

        kl=ioff(k);
        if(_scf_info->eliminate) {
          pmaxijk=pmaxij;
          if((pmaxik=maxp[(ioff(i)+k)]-2)>pmaxijk) pmaxijk=pmaxik;
          if((pmaxjk=maxp[IOFF(j,k)]-2)>pmaxijk) pmaxijk=pmaxjk;
          }

        tim_enter("l loop");
        for (l=0; l<=(k==i?j:k); l++) {
          if(use_symmetry) {
            ijkl=ioff(ij)+kl;
            nijkl=leavel=0;
            for(g=0; g < _sym_info->g ; g++) {
              gij=IOFF(_sym_info->shell_map[i][g],_sym_info->shell_map[j][g]);
              gkl=IOFF(_sym_info->shell_map[k][g],_sym_info->shell_map[l][g]);
              gijkl = IOFF(gij,gkl);
              if(gijkl > ijkl) leavel=1;
              if(gijkl == ijkl) nijkl++;
              }
            if(leavel) {
              kl++; continue;
              }
            qijkl = _sym_info->g/nijkl;
            }

          imax = (int) Qvec[kl]+Qvecij;

          if(_scf_info->eliminate) {
            cpmax = (maxp[kl]>pmaxijk) ? maxp[kl] : pmaxijk;
            if((tmax=maxp[(ioff(i)+l)]-2)>cpmax) cpmax=tmax;
            if((tmax=maxp[IOFF(j,l)]-2)>cpmax) cpmax=tmax;

            if(cpmax+imax < inttol) {
              /* If we are trying to save integrals on this node, then
               * int_index must be incremented now. */
              if (_scf_info->int_store) int_index++;
              kl++;
              continue;
              }
            }

#if 0
          if(_centers->nshell > SHELLS || use_symmetry || int_index%nproc==me) {
#endif
            tim_enter("ints");
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
            if(e_any) {
              for (bf1=0; bf1<=INT_MAX1(n1) ; bf1++) {
                i2 = _centers->func_num[s1] + bf1;

                for (bf2=0; bf2<=INT_MAX2(e12,bf1,n2) ; bf2++) {
                  j2 = _centers->func_num[s2] + bf2;
                  if(i2>=j2) { i1=i2; j1=j2; }
                  else { i1=j2; j1=i2; }
                  ij1=ioff(i1)+j1;

                  for (bf3=0; bf3<=INT_MAX3(e13e24,bf1,n3) ; bf3++) {
                    k2 = _centers->func_num[s3] + bf3;

                    for (bf4=0;bf4<=INT_MAX4(e13e24,e34,bf1,bf2,bf3,n4);bf4++){
                      if (INT_NONZERO(_mgdbuff[index])) {
                        l2 = _centers->func_num[s4] + bf4;

                        if(k2>=l2) { k1=k2; l1=l2; }
                        else { k1=l2; l1=k2; }

                        if(ij1 >= ioff(k1)+l1) {
                          ii = i1; jj = j1; kk = k1; ll = l1;
                          }
                        else {
                          ii = k1; jj = l1; kk = i1; ll = j1;
                          }

                        pki_int = (double) qijkl*_mgdbuff[index];

                        if (jj == kk) {
                          if (ii == jj || kk == ll) {
                            lij=ioff(ii)+jj;
                            lkl=ioff(kk)+ll;
                            value=(lij==lkl)? 0.25*pki_int: 0.5*pki_int;
                            gtmp[lij] += ptmp[lkl]*value;
                            gtmp[lkl] += ptmp[lij]*value;
                            if(_scf_info->hsos) {
                              gtmpo[lij] += ptmpo[lkl]*value;
                              gtmpo[lkl] += ptmpo[lij]*value;
                              }
                            }
                          else {
                            lij=ioff(ii)+jj;
                            lkl=ioff(kk)+ll;
                            value=(lij==lkl)? 0.375*pki_int: 0.75*pki_int;
                            gtmp[lij] += ptmp[lkl]*value;
                            gtmp[lkl] += ptmp[lij]*value;
                            if(_scf_info->hsos) {
                              value *= 0.3333333333333333;
                              gtmpo[lij] += ptmpo[lkl]*value;
                              gtmpo[lkl] += ptmpo[lij]*value;
                              }

                            lij=ioff(ii)+ll;
                            lkl=IOFF(kk,jj);
                            value=(lij==lkl)? 0.25*pki_int: 0.5*pki_int;
                            gtmp[lij] -= ptmp[lkl]*value;
                            gtmp[lkl] -= ptmp[lij]*value;
                            if(_scf_info->hsos) {
                              gtmpo[lij] += ptmpo[lkl]*value;
                              gtmpo[lkl] += ptmpo[lij]*value;
                              }
                            }
                          }
                        else if (ii == kk || jj == ll) {
                          lij=ioff(ii)+jj;
                          lkl=ioff(kk)+ll;
                          value=(lij==lkl)? 0.375*pki_int: 0.75*pki_int;
                          gtmp[lij] += ptmp[lkl]*value;
                          gtmp[lkl] += ptmp[lij]*value;
                          if(_scf_info->hsos) {
                            value *= 0.3333333333333333;
                            gtmpo[lij] += ptmpo[lkl]*value;
                            gtmpo[lkl] += ptmpo[lij]*value;
                            }

                          lij=ioff(ii)+kk;
                          lkl=IOFF(jj,ll);
                          value=(lij==lkl)? 0.25*pki_int : 0.5*pki_int;
                          gtmp[lij] -= ptmp[lkl]*value;
                          gtmp[lkl] -= ptmp[lij]*value;
                          if(_scf_info->hsos) {
                            gtmpo[lij] += ptmpo[lkl]*value;
                            gtmpo[lkl] += ptmpo[lij]*value;
                            }
                          }
                        else {
                          lij=ioff(ii)+jj;
                          lkl=ioff(kk)+ll;
                          value=(lij==lkl)? 0.5*pki_int : pki_int;
                          gtmp[lij] += ptmp[lkl]*value;
                          gtmp[lkl] += ptmp[lij]*value;

                          lij=ioff(ii)+kk;
                          lkl=IOFF(jj,ll);
                          value=(lij==lkl)? 0.125*pki_int: 0.25*pki_int;
                          gtmp[lij] -= ptmp[lkl]*value;
                          gtmp[lkl] -= ptmp[lij]*value;
                          if(_scf_info->hsos) {
                            gtmpo[lij] += ptmpo[lkl]*value;
                            gtmpo[lkl] += ptmpo[lij]*value;
                            }

                          if((ii != jj) && (kk != ll)) {
                            lij=ioff(ii)+ll;
                            lkl=IOFF(kk,jj);
                            value=(lij==lkl)? 0.125*pki_int: 0.25*pki_int;
                            gtmp[lij] -= ptmp[lkl]*value;
                            gtmp[lkl] -= ptmp[lij]*value;
                            if(_scf_info->hsos) {
                              gtmpo[lij] += ptmpo[lkl]*value;
                              gtmpo[lkl] += ptmpo[lij]*value;
                              }
                            }
                          }
                        }
                      index++;
                      }
                    }
                  }
                }
              }
            else {
              for (bf1=0; bf1<n1 ; bf1++) {
                i2 = _centers->func_num[s1] + bf1;

                for (bf2=0; bf2<n2 ; bf2++) {
                  j2 = _centers->func_num[s2] + bf2;
                  if(i2>=j2) { i1=i2; j1=j2; }
                  else { i1=j2; j1=i2; }
                  ij1=ioff(i1)+j1;

                  for (bf3=0; bf3<n3 ; bf3++) {
                    k2 = _centers->func_num[s3] + bf3;

                    for (bf4=0; bf4<n4; bf4++) {
                      if (INT_NONZERO(_mgdbuff[index])) {
                        l2 = _centers->func_num[s4] + bf4;

                        if(k2>=l2) { k1=k2; l1=l2; }
                        else { k1=l2; l1=k2; }

                        if(ij1 >= ioff(k1)+l1) {
                          ii = i1; jj = j1; kk = k1; ll = l1;
                          }
                        else {
                          ii = k1; jj = l1; kk = i1; ll = j1;
                          }

                        pki_int = (double) qijkl*_mgdbuff[index];

                        if (jj == kk) {
                          lij=ioff(ii)+jj;
                          lkl=ioff(kk)+ll;
                          value=0.75*pki_int;
                          gtmp[lij] += ptmp[lkl]*value;
                          gtmp[lkl] += ptmp[lij]*value;
                          if(_scf_info->hsos) {
                            value *= 0.3333333333333333;
                            gtmpo[lij] += ptmpo[lkl]*value;
                            gtmpo[lkl] += ptmpo[lij]*value;
                            }

                          lij=ioff(ii)+ll;
                          lkl=IOFF(kk,jj);
                          value=0.5*pki_int;
                          gtmp[lij] -= ptmp[lkl]*value;
                          gtmp[lkl] -= ptmp[lij]*value;
                          if(_scf_info->hsos) {
                            gtmpo[lij] += ptmpo[lkl]*value;
                            gtmpo[lkl] += ptmpo[lij]*value;
                            }
                          }
                        else if (ii == kk || jj == ll) {
                          lij=ioff(ii)+jj;
                          lkl=ioff(kk)+ll;
                          value=0.75*pki_int;
                          gtmp[lij] += ptmp[lkl]*value;
                          gtmp[lkl] += ptmp[lij]*value;
                          if(_scf_info->hsos) {
                            value *= 0.3333333333333333;
                            gtmpo[lij] += ptmpo[lkl]*value;
                            gtmpo[lkl] += ptmpo[lij]*value;
                            }

                          lij=ioff(ii)+kk;
                          lkl=IOFF(jj,ll);
                          value=0.5*pki_int;
                          gtmp[lij] -= ptmp[lkl]*value;
                          gtmp[lkl] -= ptmp[lij]*value;
                          if(_scf_info->hsos) {
                            gtmpo[lij] += ptmpo[lkl]*value;
                            gtmpo[lkl] += ptmpo[lij]*value;
                            }
                          }
                        else {
                          lij=ioff(ii)+jj;
                          lkl=ioff(kk)+ll;
                          value=pki_int;
                          gtmp[lij] += ptmp[lkl]*value;
                          gtmp[lkl] += ptmp[lij]*value;
  
                          lij=ioff(ii)+kk;
                          lkl=IOFF(jj,ll);
                          value*=0.25;
                          gtmp[lij] -= ptmp[lkl]*value;
                          gtmp[lkl] -= ptmp[lij]*value;
                          if(_scf_info->hsos) {
                            gtmpo[lij] += ptmpo[lkl]*value;
                            gtmpo[lkl] += ptmpo[lij]*value;
                            }

                          lij=ioff(ii)+ll;
                          lkl=IOFF(kk,jj);
                          gtmp[lij] -= ptmp[lkl]*value;
                          gtmp[lkl] -= ptmp[lij]*value;
                          if(_scf_info->hsos) {
                            gtmpo[lij] += ptmpo[lkl]*value;
                            gtmpo[lkl] += ptmpo[lij]*value;
                            }
                          }
                        }
                      index++;
                      }
                    }
                  }
                }
              }
            tnint += (double) (n1*n2*n3*n4);
            tim_exit("ints");
#if 0
            }
#endif
          kl++;
          int_index++;
          }
        tim_exit("l loop");
        }
      }
    }

  if(_scf_info->print_flg & 4) {
    gsum0(&tnint,1,5,mtype_get(),0);
    if(mynode0()==0)
      fprintf(_outfile,"  %8.0f integrals in scf_make_g_d\n",tnint);
    }

  free(shnfunc);

  tim_exit("scf_mkgd");

  return 0;
  }
