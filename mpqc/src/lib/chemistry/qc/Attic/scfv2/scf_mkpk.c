
/* $Log$
 * Revision 1.1  1993/12/29 12:53:04  etseidl
 * Initial revision
 *
 * Revision 1.1.1.1  1992/03/17  17:09:09  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  17:09:08  seidl
 * Initial revision
 *
 * Revision 2.1  1992/01/16  19:53:07  seidl
 * use new integral routines
 *
 * Revision 2.0  1992/01/14  19:34:15  seidl
 * totally new way of forming the supermatrix
 * the calculation of a batch is deferred until IJKL > IKJL && IJKL > ILJK
 *
 * Revision 1.8  1992/01/13  19:16:16  seidl
 * ao density formed in scf_make_gmat now
 *
 * Revision 1.7  1992/01/09  12:01:28  seidl
 * first parallel version
 *
 * Revision 1.6  1992/01/06  11:47:04  seidl
 * flush output when done
 *
 * Revision 1.5  1992/01/02  18:13:51  seidl
 * fix definition of save_thr, exit loop if imax==0
 *
 * Revision 1.4  1992/01/02  16:21:55  seidl
 * multiply diagonal elements of pmat by 2 before transforming to ao basis
 * set noi and nok = iopen rather than socc[i]
 *
 * Revision 1.3  1991/12/24  19:32:56  seidl
 * fix assignment of ptmpo
 *
 * Revision 1.2  1991/12/24  11:50:11  seidl
 * *** empty log message ***
 *
 * Revision 1.1  1991/12/20  16:23:16  seidl
 * Initial revision
 *
 * Revision 1.6  1991/11/15  14:18:51  seidl
 * *** empty log message ***
 *
 * Revision 1.5  1991/10/31  11:33:05  seidl
 * remove print statements
 *
 * Revision 1.4  1991/10/23  22:12:23  seidl
 * move includes of local and global below declaration of struct mpk
 *
 * Revision 1.3  1991/09/23  00:43:55  seidl
 * *** empty log message ***
 *
 * Revision 1.2  1991/06/22  09:04:06  seidl
 * test total angular momentum of shell quartet, put in supermatrix if
 * greater than maxam
 *
 * Revision 1.1  1991/06/18  12:51:31  seidl
 * Initial revision
 * */

static char *rcsid = "$Id$";

#include <stdio.h>
#include <tmpl.h>
#include <math/array/math_lib.h>
#include <math.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <chemistry/qc/sym/sym_lib.h>

#include "scf.h"

#ifdef INT_CE_SH_AM
#undef INT_CE_SH_AM
#endif
#define INT_CE_SH_AM(c,a,s) ((c)->center[(a)].basis.shell[(s)].type[0].am)

#define ioff(i) i*(i+1)/2
#define MIN0(a,b) ((a)<(b)) ? (a) : (b)
#define MAX0(a,b) ((a)>(b)) ? (a) : (b)

#ifndef IOFF
#define IOFF(a,b) ((a)>(b))?(a)*((a)+1)/2+(b):(b)*((b)+1)/2+(a)
#endif

static int sijtest,skltest;

struct mpk {
    int intmx;
    int old_nint;
    int iblc;
    int iblo;
    int num_ints_o;
    int num_ints_c;
    int keep;
    int keep2;
    int nint;
    int scr_c;
    int scr_o;
    int *inext;
    unsigned int *lbij;
    unsigned int *lbkl;
    double *pa;
    double *pb;
    double *gtmp;
    double *gtmpo;
    double *ptmp;
    double *ptmpo;
    FILE *outfile;
    } ;

#include "scf_mkpk.gbl"
#include "scf_mkpk.lcl"

GLOBAL_FUNCTION int
scf_make_pk(_centers,_irreps,_scf_info,_sym_info,_mfp,_mast,
       _gmat,_gmato,_dpmat,_dpmato,scr_c,scr_o,spkc,spko,buffer,_outfile)
centers_t *_centers;
scf_irreps_t *_irreps;
scf_struct_t *_scf_info;
sym_struct_t *_sym_info;
scf_mf_ptrs_t *_mfp;
int _mast;
double_vector_t *_gmat;
double_vector_t *_gmato;
double_vector_t *_dpmat;
double_vector_t *_dpmato;
int scr_c;
int scr_o;
scf_pk_buf_c_t *spkc;
scf_pk_buf_o_t *spko;
double *buffer;
FILE *_outfile;
{
  int errcod;
  int iopen;
  int i,j,k,l;
  int ij,kl,ijkl;
  int ik,jl,ikjl;
  int il,jk,iljk;
  struct mpk pkb;

  iopen = _scf_info->iopen;

#if 0
  if(_scf_info->nbatri > 150) {
    pkb.intmx = pkb.old_nint = 25920*4;
    pkb.keep = 1024;
    pkb.keep2 =1023;
    }
  else {
    pkb.intmx = pkb.old_nint = 25920;
    pkb.keep = 128;
    pkb.keep2 = 127;
    }
#else
  pkb.intmx = pkb.old_nint = 25920;
  pkb.keep = 128;
  pkb.keep2 = 127;
#endif
  pkb.iblc = 0;
  pkb.iblo = 0;
  pkb.num_ints_c = 0;
  pkb.num_ints_o = 0;
  pkb.nint=0;
  pkb.scr_c=scr_c;
  pkb.scr_o=scr_o;

  pkb.outfile = _outfile;

  pkb.pa = (double *) init_array(sizeof(double)*pkb.intmx);
  pkb.pb = (double *) init_array(sizeof(double)*pkb.intmx);
  pkb.lbij = (unsigned int *) init_array(sizeof(int)*pkb.intmx);
  pkb.lbkl = (unsigned int *) init_array(sizeof(int)*pkb.intmx);
  pkb.inext = (int *) init_array(sizeof(int)*(pkb.keep+pkb.intmx));

  pkb.ptmp = _dpmat->d;
  pkb.ptmpo = _dpmato->d;

  pkb.gtmp = _gmat->d;
  pkb.gtmpo = _gmato->d;

  for (i=0; i<_centers->nshell; i++) {
    int use_symmetry=(_sym_info->g > 1)?1:0;

    if(use_symmetry) if(!_sym_info->p1[i]) continue;

    for (j=0; j<=i; j++) {
      ij = IOFF(i,j);

      for (k=0; k<=i; k++) {
        ik = IOFF(i,k);
        jk = IOFF(j,k);

        for (l=0; l<=(k==i?j:k); l++) {

#if defined(I860) || defined(NCUBE)
    if(l%numnodes0() != mynode0()) continue;
#endif

          il=IOFF(i,l);
          jl=IOFF(j,l);
          kl=IOFF(k,l);
          ijkl=IOFF(ij,kl);
          ikjl=IOFF(ik,jl);
          iljk=IOFF(il,jk);

          if((ikjl>ijkl) || (iljk>ijkl)) continue;

          pk_int_loop(_centers,_scf_info,_irreps,_sym_info,&pkb,buffer,
               spkc,spko,i,j,k,l);
          }
        }
      }
    }

  packit(&pkb,_irreps,_scf_info,spkc,spko,1);
#if !defined(I860) && !defined(NCUBE)
  fflush(_outfile);
#endif

  free(pkb.pa);
  free(pkb.pb);
  free(pkb.lbij);
  free(pkb.lbkl);
  free(pkb.inext);
  return(0);
  }


LOCAL_FUNCTION VOID
filter_int(ism,jsm,ksm,lsm,ior,jor,kor,lor,ii,jj,kk,ll,pki_int,
           pkb,_scf_info,_irreps)
int ism;
int jsm;
int ksm;
int lsm;
int ior;
int jor;
int kor;
int lor;
int ii;
int jj;
int kk;
int ll;
double pki_int;
struct mpk *pkb;
scf_struct_t *_scf_info;
scf_irreps_t *_irreps;
{
  double fabs();
  int p1,p2;
                                                         
#if 0
  if (ism == jsm && ksm == lsm && ism == ksm) {
#endif
    if ((jor==kor) && ((ior == jor) || (kor == lor))) {
      findit(_scf_info,_irreps,pkb,ii,jj,kk,ll,ism,ksm,pki_int,5);
      }
    else if (ior == kor || jor == lor) {
      findit(_scf_info,_irreps,pkb,ii,jj,kk,ll,ism,ksm,pki_int,3);

      findit(_scf_info,_irreps,pkb,ii,kk,jj,ll,ism,ksm,pki_int,4);
      }
    else if (jor == kor) {
      findit(_scf_info,_irreps,pkb,ii,jj,kk,ll,ism,ksm,pki_int,3);

      findit(_scf_info,_irreps,pkb,ii,ll,jj,kk,ism,ksm,pki_int,4);
      }
    else if (ior == jor || kor == lor) {
      findit(_scf_info,_irreps,pkb,ii,jj,kk,ll,ism,ksm,pki_int,1);

      findit(_scf_info,_irreps,pkb,ii,kk,jj,ll,ism,ksm,pki_int,2);
      }
    else {
      findit(_scf_info,_irreps,pkb,ii,jj,kk,ll,ism,ksm,pki_int,1);

      findit(_scf_info,_irreps,pkb,ii,kk,jj,ll,ism,ksm,pki_int,2);

      findit(_scf_info,_irreps,pkb,ii,ll,jj,kk,ism,ksm,pki_int,2);
      }
#if 0
    }
  else if (ism == jsm) {
    findit(_scf_info,_irreps,pkb,ii,jj,kk,ll,ism,ksm,pki_int,1);
    }
  else if (ism == ksm) {
    if (ior == kor || jor == lor) {
      p1 = MAX0(jj,ll);
      p2 = MIN0(jj,ll);

      findit(_scf_info,_irreps,pkb,ii,kk,p1,p2,ism,jsm,pki_int,4);
      }
    else {
      p1 = MAX0(jj,ll);
      p2 = MIN0(jj,ll);

      findit(_scf_info,_irreps,pkb,ii,kk,p1,p2,ism,jsm,pki_int,2);
      }
    }
#endif
  }

LOCAL_FUNCTION VOID
findit(_scf_info,_irreps,pkb,ii,jj,kk,ll,ism,ksm,value,iab)
scf_struct_t *_scf_info;
scf_irreps_t *_irreps;
struct mpk *pkb;
int ii;
int jj;
int kk;
int ll;
int ism;
int ksm;
double value;
int iab;
{
  register int i,j;
  unsigned int *ijtmp, *kltmp;
  int *nxtmp;
  int p1,p2,p3;
  int noi,nok;
  int next,start;
  int lij,lkl;
  double *patmp, *pbtmp;
  scf_irrep_t *irr_i,*irr_k;

  irr_i = _irreps[ism].ir;
  irr_k = _irreps[ksm].ir;

#if 0
  noi = irr_i->nopen;
  nok = irr_k->nopen;
#else
  noi = _scf_info->iopen;
  nok = _scf_info->iopen;
#endif

  lij = IOFF(ii,jj);
  lkl = IOFF(kk,ll);

  if(!pkb->nint) {
    bzero(pkb->inext,sizeof(int) * pkb->old_nint);
    bzero(&pkb->inext[pkb->intmx],sizeof(int) * pkb->keep);
    }

  start = 2*lij + lkl;
  start = (start & pkb->keep2) + pkb->intmx;

L1:
  next=pkb->inext[start];
  if(next) {
    if (pkb->lbij[next-1] == lij && pkb->lbkl[next-1] == lkl) i=next-1;
    else {
      start = next;
      goto L1;
      }
    }
  else {
    i=pkb->nint;
    if(pkb->nint >= pkb->intmx) {
      fprintf(pkb->outfile,"\n  increasing size of buffers in findit\n");
      fprintf(pkb->outfile,"  intmx was %d, is %d\n",pkb->intmx,pkb->intmx*2);
      fflush(pkb->outfile);
      pkb->intmx*=2;

 /* i don't use realloc because strange things were happening */

      nxtmp = (int *) init_array(sizeof(int)*(pkb->keep+pkb->intmx));
      bcopy(pkb->inext,nxtmp,(int)sizeof(int)*(pkb->intmx/2));
      for(j=0; j < pkb->keep ; j++) 
        nxtmp[j+pkb->intmx]=pkb->inext[j+pkb->intmx/2];
      free(pkb->inext);
      pkb->inext=nxtmp;

      ijtmp = (unsigned int *) init_array(sizeof(int)*pkb->intmx);
      bcopy(pkb->lbij,ijtmp,sizeof(int)*(pkb->intmx/2));
      free(pkb->lbij);
      pkb->lbij = ijtmp;

      kltmp = (unsigned int *) init_array(sizeof(int)*pkb->intmx);
      bcopy(pkb->lbkl,kltmp,sizeof(int)*(pkb->intmx/2));
      free(pkb->lbkl);
      pkb->lbkl = kltmp;

      patmp = (double *) init_array(sizeof(double)*pkb->intmx);
      bcopy(pkb->pa,patmp,sizeof(double)*(pkb->intmx/2));
      free(pkb->pa);
      pkb->pa = patmp;

      pbtmp = (double *) init_array(sizeof(double)*pkb->intmx);
      bcopy(pkb->pb,pbtmp,sizeof(double)*(pkb->intmx/2));
      free(pkb->pb);
      pkb->pb = pbtmp;

      if(pkb->inext==NULL || pkb->lbij==NULL || pkb->lbkl==NULL) {
        fprintf(pkb->outfile,
         "\n pathological problems with realloc in findit\n");
        fprintf(pkb->outfile," try upping intmx to %d\n",pkb->intmx);
        fflush(pkb->outfile);
        exit(1);
        }

      start = 2*lij + lkl;
      start = (start & pkb->keep2) + pkb->intmx;
      }
    pkb->inext[start] = ++pkb->nint;
    pkb->lbij[i] = lij;
    pkb->lbkl[i] = lkl;
    pkb->pa[i] = pkb->pb[i] = 0.0;
    }

  value = (lij == lkl) ? value*0.5 : value;

  switch(iab) {
  case 1:
    pkb->pa[i] += value;
    if(noi && nok) {
      if(_scf_info->special) {
        p3 = IOFF(irr_i->os_num,irr_k->os_num);
        pkb->pb[i] += _scf_info->alpha[p3]*value;
        }
      }
    break;
  case 2:
    pkb->pa[i] -= 0.25*value;
    if(noi && nok) {
      if(_scf_info->hsos) pkb->pb[i] += 0.25*value;
      else if(_scf_info->singlet) {
        if (ism != ksm) pkb->pb[i] -= 0.75*value;
        else pkb->pb[i] += 0.25*value;
        }
      else if(_scf_info->twocon) {
        if (ism != ksm) pkb->pb[i] += 0.25*value;
        }
      else {
        p3 = IOFF(irr_i->os_num,irr_k->os_num);
        pkb->pb[i] += 0.25*_scf_info->beta[p3]*value;
        }
      }
    break;
  case 3:
    pkb->pa[i] += 0.75*value;
    if(noi && nok) {
      if(_scf_info->hsos) pkb->pb[i] += 0.25*value;
      else if(_scf_info->singlet) {
        if (ism != ksm) pkb->pb[i] -= 0.75*value;
        else pkb->pb[i] += 0.25*value;
        }
      else if(_scf_info->twocon) {
        if (ism != ksm) pkb->pb[i] += 0.25*value;
        }
      else {
        p3 = IOFF(irr_i->os_num,irr_k->os_num);
        pkb->pb[i] += (_scf_info->alpha[p3] + 0.25*_scf_info->beta[p3])*value;
        }
      }
    break;
  case 4:
    pkb->pa[i] -= 0.5*value;
    if(noi && nok) {
      if(_scf_info->hsos) pkb->pb[i] = 0.5*value;
      else if(_scf_info->singlet) {
        if (ism != ksm) pkb->pb[i] -= 1.5*value;
        else pkb->pb[i] += 0.5*value;
        }
      else if(_scf_info->twocon) {
        if (ism != ksm) pkb->pb[i] += 0.5*value;
        }
      else {
        p3 = IOFF(irr_i->os_num,irr_k->os_num);
        pkb->pb[i] += 0.5*_scf_info->beta[p3]*value;
        }
      }
    break;
  case 5:
    pkb->pa[i] = 0.5*value;
    if(noi && nok) {
      if(_scf_info->hsos) pkb->pb[i] = 0.5*value;
      else if(_scf_info->singlet) {
        if (ism != ksm) pkb->pb[i] = -1.5*value;
        else pkb->pb[i] = 0.5*value;
        }
      else if(_scf_info->twocon) {
        if (ism != ksm) pkb->pb[i] += 0.5*value;
        }
      else {
        p3 = IOFF(irr_i->os_num,irr_k->os_num);
        pkb->pb[i] = (_scf_info->alpha[p3] + 0.5*_scf_info->beta[p3])*value;
        }
      }
    }
  pkb->old_nint=pkb->nint;
  }


LOCAL_FUNCTION VOID
packit(pkb,_irreps,_scf_info,spkc,spko,endflg)
struct mpk *pkb;
scf_irreps_t *_irreps;
scf_struct_t *_scf_info;
scf_pk_buf_c_t *spkc;
scf_pk_buf_o_t *spko;
int endflg;
{
  int i,j,k,ij,kl;
  double pval, kval;
  double tol = 10e-14;


  if(!endflg) {
    for(i=0; i < pkb->nint ; i++) {
      pval=pkb->pa[i];
      kval=pkb->pb[i];
      ij = pkb->lbij[i];
      kl = pkb->lbkl[i];
      if (fabs(pval) >= tol || fabs(kval) >= tol) {
        if(fabs(kval) >= tol) {
          spko->pk[pkb->iblo].ij = ij;
          spko->pk[pkb->iblo].kl = kl;
          spko->pk[pkb->iblo].pval = pval;
          spko->pk[pkb->iblo].kval = kval;

          if(!_scf_info->twocon) {
            pkb->gtmp[ij] += pkb->ptmp[kl]*pval;
            pkb->gtmp[kl] += pkb->ptmp[ij]*pval;
            pkb->gtmpo[ij] += pkb->ptmpo[kl]*kval;
            pkb->gtmpo[kl] += pkb->ptmpo[ij]*kval;
            }

          pkb->iblo++;

          if (pkb->iblo >= _scf_info->maxbufo) {
            _scf_info->readflgo=1;
            bio_write(pkb->scr_o,spko->pk,
                sizeof(scf_pkint_o_t)*_scf_info->maxbufo);
            pkb->num_ints_o += pkb->iblo;
            _scf_info->num_bufs_o++;
            pkb->iblo=0;
            }
          }
        else {
          spkc->pk[pkb->iblc].ij = ij;
          spkc->pk[pkb->iblc].kl = kl;
          spkc->pk[pkb->iblc].pval = pval;

          if(!_scf_info->twocon) {
            pkb->gtmp[ij] += pkb->ptmp[kl]*pval;
            pkb->gtmp[kl] += pkb->ptmp[ij]*pval;
            }

          pkb->iblc++;

          if (pkb->iblc >= _scf_info->maxbufc) {
            _scf_info->readflgc=1;
            bio_write(pkb->scr_c,spkc->pk,
              sizeof(scf_pkint_c_t)*_scf_info->maxbufc);
            pkb->num_ints_c += pkb->iblc;
            _scf_info->num_bufs_c++;
            pkb->iblc=0;
            }
          }
        }
      }
    pkb->nint=0;
    }
  else {
    pkb->num_ints_o += pkb->iblo;
    pkb->num_ints_c += pkb->iblc;
    _scf_info->num_bufs_c++;
    _scf_info->num_bufs_o++;
#if !defined(I860) && !defined(NCUBE)
    fprintf(pkb->outfile,"%10d closed-shell integrals written to disk",
                                     pkb->num_ints_c);
    fprintf(pkb->outfile," in %3d buffers\n",_scf_info->num_bufs_c);
    fprintf(pkb->outfile,"%10d open-shell integrals written to disk",
                                     pkb->num_ints_o);
    fprintf(pkb->outfile," in %3d buffers\n",_scf_info->num_bufs_o);
#endif
                                     
    _scf_info->lasto = pkb->iblo;
    _scf_info->lastc = pkb->iblc;
    if(_scf_info->readflgo)
      bio_write(pkb->scr_o,spko->pk,sizeof(scf_pkint_o_t)*pkb->iblo);
    if(_scf_info->readflgc)
      bio_write(pkb->scr_c,spkc->pk,sizeof(scf_pkint_c_t)*pkb->iblc);
    }
  }

LOCAL_FUNCTION VOID_PTR
init_array(n)
int n;
{
  void *arr;

  arr = (void *) malloc(n);
  if(arr==NULL) {
    fprintf(stderr,"scf_make_pk: init_array:\n");
    fprintf(stderr,"malloc trouble\n");
    exit(99);
    }
  bzero(arr,n);
  return arr;
  }


LOCAL_FUNCTION VOID
pk_int_loop(_centers,_scf_info,_irreps,_sym_info,
   pkb,buffer,spkc,spko,i,j,k,l)
centers_t *_centers;
scf_struct_t *_scf_info;
scf_irreps_t *_irreps;
sym_struct_t *_sym_info;
struct mpk *pkb;
double *buffer;
scf_pk_buf_c_t *spkc;
scf_pk_buf_o_t *spko;
int i;
int j;
int k;
int l;
{
  int g,m,leavel,nijkl,qijkl=1;
  int num,tam,tl;
  int s1,s2,s3,s4;
  int n1,n2,n3,n4;
  int e12,e34,e13e24;
  int bf1,bf2,bf3,bf4;
  int i1,j1,k1,l1;
  int i2,j2,k2,l2;
  int ii,jj,kk,ll;
  int ism,jsm,ksm,lsm;
  int ior,jor,kor,lor;
  int ij,kl,ijkl;
  int gi,gj,gk,gl,gij,gkl,gijkl;
  int index;
  int test_size;
  int use_symmetry=(_sym_info->g >1);
  double pki_int;
  double save_thr = pow(10.0,(double) -(_scf_info->save_thr));
  double imax;

  ism=jsm=ksm=lsm=0;

  tam = INT_CE_SH_AM((_centers),(_centers->center_num[i]),
                                (_centers->shell_num[i]));
  tam += INT_CE_SH_AM((_centers),(_centers->center_num[j]),
                                (_centers->shell_num[j]));
  tam += INT_CE_SH_AM((_centers),(_centers->center_num[k]),
                                (_centers->shell_num[k]));
  tam += INT_CE_SH_AM((_centers),(_centers->center_num[l]),
                                (_centers->shell_num[l]));

  tl=6;
  if(tam==0) tl=1;
  else if(tam==1) tl=2;
  else if(tam==2) tl=3;

  tl *= _centers->center[_centers->center_num[i]].basis.shell[
                                       _centers->shell_num[i]].nprim;
  tl *= _centers->center[_centers->center_num[j]].basis.shell[
                                       _centers->shell_num[j]].nprim;
  tl *= _centers->center[_centers->center_num[k]].basis.shell[
                                       _centers->shell_num[k]].nprim;
  tl *= _centers->center[_centers->center_num[l]].basis.shell[
                                       _centers->shell_num[l]].nprim;

  test_size = (tl < _scf_info->expense);

  if((j==k) && ((i==j) || (k==l))) num=1;
  else if((k==l)||(j==k)||(j==l)||(i==k)||(i==j)) num=2;
  else num=3;

  for(m=0; m < num ; m++) {
    switch(m) {
    case 0:
      s1=i; s2=j; s3=k; s4=l;
      break;
    case 1:
      if(j==k) {
        s1=i; s2=l; s3=k; s4=j;
        }
      else {
        s1=i; s2=k; s3=j; s4=l;
        }
      break;
    case 2:
      s1=i; s2=l; s3=k; s4=j;
      }

    if(use_symmetry) {
      ij=IOFF(s1,s2);
      if(!_sym_info->lamij[ij]) continue;
      }

    imax = fabs(int_erep_bound(0,s1,s2,s3,s4));

    if(test_size && imax < save_thr) continue;

    if(use_symmetry) {
      kl=IOFF(s3,s4);
      ijkl=IOFF(ij,kl);

      nijkl=leavel=0;
      for(g=0; g < _sym_info->g ; g++) {
        gi = _sym_info->shell_map[s1][g];
        gj = _sym_info->shell_map[s2][g];
        gk = _sym_info->shell_map[s3][g];
        gl = _sym_info->shell_map[s4][g];
        gij = IOFF(gi,gj);
        gkl = IOFF(gk,gl);
        gijkl = IOFF(gij,gkl);
        if(gijkl > ijkl) leavel=1;
        if(gijkl == ijkl) nijkl++;
        }
      if(leavel) continue;
      qijkl = _sym_info->g/nijkl;
      }

    int_erep(INT_EREP|INT_NOBCHK,&s1,&s2,&s3,&s4);

    n1 = INT_SH_NFUNC((_centers),s1);
    n2 = INT_SH_NFUNC((_centers),s2);
    n3 = INT_SH_NFUNC((_centers),s3);
    n4 = INT_SH_NFUNC((_centers),s4);

 /* Shell equivalency information. */
    e12    = (s2==s1);
    e13e24 = (s3==s1) && (s4==s2);
    e34    = (s4==s3);

    index = 0;
    for (bf1=0; bf1<=INT_MAX1(n1); bf1++) {
      for (bf2=0; bf2<=INT_MAX2(e12,bf1,n2); bf2++) {
        for (bf3=0; bf3<=INT_MAX3(e13e24,bf1,n3); bf3++) {
          for (bf4=0; bf4<=INT_MAX4(e13e24,e34,bf1,bf2,bf3,n4); bf4++) {
            if (INT_NONZERO(buffer[index])) {
              i2 = _centers->func_num[s1] + bf1;
              j2 = _centers->func_num[s2] + bf2;
              k2 = _centers->func_num[s3] + bf3;
              l2 = _centers->func_num[s4] + bf4;

              i1=(i2 >= j2) ? i2 : j2;
              j1=(i2 >= j2) ? j2 : i2;
              k1=(k2 >= l2) ? k2 : l2;
              l1=(k2 >= l2) ? l2 : k2;

              if((i1*(i1+1)/2+j1) >= (k1*(k1+1)/2+l1)) {
                ior = ii = i1;
                jor = jj = j1;
                kor = kk = k1;
                lor = ll = l1;
                }
              else {
                ior = ii = k1;
                jor = jj = l1;
                kor = kk = i1;
                lor = ll = j1;
                }
  
              pki_int = (double) qijkl*buffer[index];

              filter_int(ism,jsm,ksm,lsm,ior,jor,kor,lor,
                         ii,jj,kk,ll,pki_int,pkb,_scf_info,_irreps);
              }
            index++;
            }
          }
        }
      }
    }
  packit(pkb,_irreps,_scf_info,spkc,spko,0);
  }
