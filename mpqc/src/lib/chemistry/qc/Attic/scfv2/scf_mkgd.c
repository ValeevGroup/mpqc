
/* Calculates two-electron integrals on the fly and sticks them into the
 * appropriate part of the G matrix
 */

/* $Log$
 * Revision 1.1  1993/12/29 12:53:04  etseidl
 * Initial revision
 *
 * Revision 1.2  1992/06/17  22:14:15  jannsen
 * clean up for saber-c
 *
 * Revision 1.1.1.1  1992/03/17  17:09:06  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  17:09:05  seidl
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
#include <tmpl.h>
#include <math/array/math_lib.h>
#include <math.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <chemistry/qc/sym/sym_lib.h>
#include <comm/picl/picl.h>
#include "scf.h"

#ifdef INT_CE_SH_AM
#undef INT_CE_SH_AM
#endif
#define INT_CE_SH_AM(c,a,s) ((c)->center[(a)].basis.shell[(s)].type[0].am)

#define ioff(i) (i)*((i)+1)/2
#define MIN0(a,b) ((a)<(b)) ? (a) : (b)
#define MAX0(a,b) ((a)>(b)) ? (a) : (b)

#ifndef IOFF
#define IOFF(a,b) ((a)>(b))?(a)*((a)+1)/2+(b):(b)*((b)+1)/2+(a)
#endif

#define NEED_SHELL 10001
#define SHELL_INDEX 10002

struct mgd {
    double *gtmp;
    double *gtmpo;
    double *ptmp;
    double *ptmpo;
    FILE *outfile;
    } ;

#include "scf_mkgd.gbl"
#include "scf_mkgd.lcl"

GLOBAL_FUNCTION int
scf_make_g_d(_centers,_irreps,_scf_info,_sym_info,_mfp,_mast,
                 _gmat,_gmato,_dpmat,_dpmato,_mgdbuff,iter,_outfile)
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
double *_mgdbuff;
int iter;
FILE *_outfile;
{
  int index;
  int errcod;
  int iopen;
  int s1,s2,s3,s4;
  int n1,n2,n3,n4;
  int am1,am2,am3,am4,tam;
  int nc1,nc2,nc3,nc4,tl;
  int i,j,k,l;
  int ij,kl,ijkl;
  int si,sj,sij;
  int g,gi,gj,gij,gk,gl,gkl,gijkl;
  int nijkl,leavel,qijkl=1;
  int joff,nn;
  int tnint=0;
  int tshell=0;
  int indices[4];
  double pki_int;
  scf_irrep_t *s;
  struct mgd gdb;
  double pmax,imax;
  double inttol = pow(10.0,(double) -(_scf_info->intcut));
  double save_thr = pow(10.0,(double) -(_scf_info->save_thr));
  int use_symmetry=(_sym_info->g > 1)?1:0;

  iopen = _scf_info->iopen;

  gdb.outfile = _outfile;

  gdb.gtmp = _gmat->d;
  gdb.ptmp = _dpmat->d;
  if(iopen) {
    gdb.gtmpo = _gmato->d;
    gdb.ptmpo = _dpmato->d;
    }

 /* loop over all shells, calculate a bunch of integrals from each shell
  * quartet, and stick those integrals where they belong
  */

#if defined(I860) || defined(NCUBE)

 /* for parallel architectures, if only a few nodes are used, then split
  * work up in a very simple way */

  if(numnodes0()==1 || numnodes0()==2) {
    for (i=0; i<_centers->nshell; i++) {
      if(use_symmetry) if(!_sym_info->p1[i]) continue;

      am1 = INT_CE_SH_AM((_centers),(_centers->center_num[i]),
                                        (_centers->shell_num[i]));
      nc1 = _centers->center[_centers->center_num[i]].basis.shell[
                                               _centers->shell_num[i]].nprim;

      for (j=0; j<=i; j++) {
        leavel=0;
        if(use_symmetry) {
          ij = IOFF(i,j);
          if(!_sym_info->lamij[ij]) continue;
          }

        am2 = INT_CE_SH_AM((_centers),(_centers->center_num[j]),
                                        (_centers->shell_num[j]));
        nc2 = _centers->center[_centers->center_num[j]].basis.shell[
                                               _centers->shell_num[j]].nprim;

        for (k=0; k<=i; k++) {
          am3 = INT_CE_SH_AM((_centers),(_centers->center_num[k]),
                                        (_centers->shell_num[k]));
          nc3 = _centers->center[_centers->center_num[k]].basis.shell[
                                               _centers->shell_num[k]].nprim;
          for (l=0; l<=(k==i?j:k); l++) {

            if(l%numnodes0() != mynode0()) continue;

            if(use_symmetry) {
              kl=IOFF(k,l);
              ijkl=IOFF(ij,kl);
              nijkl=leavel=0;
              for(g=0; g < _sym_info->g ; g++) {
                gi = _sym_info->shell_map[i][g];
                gj = _sym_info->shell_map[j][g];
                gk = _sym_info->shell_map[k][g];
                gl = _sym_info->shell_map[l][g];
                gij = IOFF(gi,gj);
                gkl = IOFF(gk,gl);
                gijkl = IOFF(gij,gkl);
                if(gijkl > ijkl) leavel=1;
                if(gijkl == ijkl) nijkl++;
                }
              if(leavel) continue;
              qijkl = _sym_info->g/nijkl;
              }

            s1 = i;
            s2 = j;
            s3 = k;
            s4 = l;

            am4 = INT_CE_SH_AM((_centers),(_centers->center_num[s4]),
                                        (_centers->shell_num[s4]));
            nc4 = _centers->center[_centers->center_num[s4]].basis.shell[
                                               _centers->shell_num[s4]].nprim;

            imax = fabs(int_erep_bound(0,s1,s2,s3,s4));
            if(imax==0.0) continue;

            tam=am1+am2+am3+am4;
            tl=6;
            if(tam==0) tl=1;
            else if(tam==1) tl=2;
            else if(tam==2) tl=3;

            if((tl*nc1*nc2*nc3*nc4) >= _scf_info->expense ||
               imax >= save_thr) continue;

            pmax = fabs(max_den(_centers,gdb.ptmp,s1,s2,s3,s4));
            if(pmax*imax < inttol) continue;

            int_erep(INT_EREP|INT_NOBCHK,&s1,&s2,&s3,&s4);

            mgd_int_loop(_centers,_scf_info,_irreps,&gdb,_mgdbuff,
              s1,s2,s3,s4,qijkl,&tnint);
            }
          }
        }
      }
    }
  else if(mynode0()==0) {

 /* for a large number of nodes, have node 0 pass work out to all the 
  * other nodes
  */
    int iln=sizeof(int);
    int node;
    int asdf=0;

    for (i=_centers->nshell-1; i>-1; i--) {
      if(use_symmetry) if(!_sym_info->p1[i]) continue;
      indices[0]=i;

      for (j=i; j>-1; j--) {
        leavel=0;
        if(use_symmetry) {
          ij = IOFF(i,j);
          if(!_sym_info->lamij[ij]) continue;
          }
        indices[1]=j;

        recv0(&node,iln,NEED_SHELL);
        send0(indices,iln*2,SHELL_INDEX,node);
        }
      }

    indices[0]=indices[1]= -1;
    for (k=1; k < numnodes0() ; k++) {
      recv0(&node,iln,NEED_SHELL);
      send0(indices,iln*2,SHELL_INDEX,node);
      }
    }
  else {
    int iln=sizeof(int);
    int asdf=0;
    do {
      i=mynode0();
      send0(&i,iln,NEED_SHELL,0);
      recv0(indices,iln*2,SHELL_INDEX);
      i=indices[0];
      j=indices[1];

      if(i==-1) continue;

      tshell++;

      am1 = INT_CE_SH_AM((_centers),(_centers->center_num[i]),
                                        (_centers->shell_num[i]));
      nc1 = _centers->center[_centers->center_num[i]].basis.shell[
                                               _centers->shell_num[i]].nprim;

      am2 = INT_CE_SH_AM((_centers),(_centers->center_num[j]),
                                      (_centers->shell_num[j]));
      nc2 = _centers->center[_centers->center_num[j]].basis.shell[
                                             _centers->shell_num[j]].nprim;
      ij = IOFF(i,j);

      for (k=0; k<=i; k++) {
        am3 = INT_CE_SH_AM((_centers),(_centers->center_num[k]),
                                      (_centers->shell_num[k]));
        nc3 = _centers->center[_centers->center_num[k]].basis.shell[
                                             _centers->shell_num[k]].nprim;
        for (l=0; l<=(k==i?j:k); l++) {
          if(use_symmetry) {
            kl=IOFF(k,l);
            ijkl=IOFF(ij,kl);
            nijkl=leavel=0;
            for(g=0; g < _sym_info->g ; g++) {
              gi = _sym_info->shell_map[i][g];
              gj = _sym_info->shell_map[j][g];
              gk = _sym_info->shell_map[k][g];
              gl = _sym_info->shell_map[l][g];
              gij = IOFF(gi,gj);
              gkl = IOFF(gk,gl);
              gijkl = IOFF(gij,gkl);
              if(gijkl > ijkl) leavel=1;
              if(gijkl == ijkl) nijkl++;
              }
            if(leavel) continue;
            qijkl = _sym_info->g/nijkl;
            }

          s1 = i;
          s2 = j;
          s3 = k;
          s4 = l;

          am4 = INT_CE_SH_AM((_centers),(_centers->center_num[s4]),
                                      (_centers->shell_num[s4]));
          nc4 = _centers->center[_centers->center_num[s4]].basis.shell[
                                             _centers->shell_num[s4]].nprim;

          imax = fabs(int_erep_bound(0,s1,s2,s3,s4));
          if(imax==0.0) continue;

          tam=am1+am2+am3+am4;
          tl=6;
          if(tam==0) tl=1;
          else if(tam==1) tl=2;
          else if(tam==2) tl=3;

          if((tl*nc1*nc2*nc3*nc4) >= _scf_info->expense ||
             imax >= save_thr) continue;

          pmax = fabs(max_den(_centers,gdb.ptmp,s1,s2,s3,s4));
          if(pmax*imax < inttol) continue;

          int_erep(INT_EREP|INT_NOBCHK,&s1,&s2,&s3,&s4);

          mgd_int_loop(_centers,_scf_info,_irreps,&gdb,_mgdbuff,
            s1,s2,s3,s4,qijkl,&tnint);
          }
        }
      } while(i != -1);
    }

#else /* !I860 && !NCUBE */
  for (i=0; i<_centers->nshell; i++) {
    int use_symmetry=(_sym_info->g > 1)?1:0;
    if(use_symmetry) if(!_sym_info->p1[i]) continue;

    am1 = INT_CE_SH_AM((_centers),(_centers->center_num[i]),
                                        (_centers->shell_num[i]));
    nc1 = _centers->center[_centers->center_num[i]].basis.shell[
                                               _centers->shell_num[i]].nprim;
    for (j=0; j<=i; j++) {
      leavel=0;
      if(use_symmetry) {
        ij = IOFF(i,j);
        if(!_sym_info->lamij[ij]) continue;
        }

      am2 = INT_CE_SH_AM((_centers),(_centers->center_num[j]),
                                        (_centers->shell_num[j]));
      nc2 = _centers->center[_centers->center_num[j]].basis.shell[
                                               _centers->shell_num[j]].nprim;
      for (k=0; k<=i; k++) {
        am3 = INT_CE_SH_AM((_centers),(_centers->center_num[k]),
                                        (_centers->shell_num[k]));
        nc3 = _centers->center[_centers->center_num[k]].basis.shell[
                                               _centers->shell_num[k]].nprim;
        for (l=0; l<=(k==i?j:k); l++) {
          if(use_symmetry) {
            kl=IOFF(k,l);
            ijkl=IOFF(ij,kl);
            nijkl=leavel=0;
            for(g=0; g < _sym_info->g ; g++) {
              gi = _sym_info->shell_map[i][g];
              gj = _sym_info->shell_map[j][g];
              gk = _sym_info->shell_map[k][g];
              gl = _sym_info->shell_map[l][g];
              gij = IOFF(gi,gj);
              gkl = IOFF(gk,gl);
              gijkl = IOFF(gij,gkl);
              if(gijkl > ijkl) leavel=1;
              if(gijkl == ijkl) nijkl++;
              }
            if(leavel) continue;
            qijkl = _sym_info->g/nijkl;
            }

          s1 = i;
          s2 = j;
          s3 = k;
          s4 = l;

          am4 = INT_CE_SH_AM((_centers),(_centers->center_num[s4]),
                                        (_centers->shell_num[s4]));
          nc4 = _centers->center[_centers->center_num[s4]].basis.shell[
                                               _centers->shell_num[s4]].nprim;

          imax = fabs(int_erep_bound(0,s1,s2,s3,s4));
          if(imax==0.0) continue;

          tam=am1+am2+am3+am4;
          tl=6;
          if(tam==0) tl=1;
          else if(tam==1) tl=2;
          else if(tam==2) tl=3;

          if((tl*nc1*nc2*nc3*nc4) >= _scf_info->expense ||
             imax >= save_thr) continue;

          pmax = fabs(max_den(_centers,gdb.ptmp,s1,s2,s3,s4));
          if(pmax*imax < inttol) continue;

          int_erep(INT_EREP|INT_NOBCHK,&s1,&s2,&s3,&s4);

          mgd_int_loop(_centers,_scf_info,_irreps,&gdb,_mgdbuff,
            s1,s2,s3,s4,qijkl,&tnint);
          }
        }
      }
    }

  fprintf(_outfile,"%10d integrals in scf_make_g_d\n",tnint);
  fflush(_outfile);
#endif /* I860 || NCUBE */

  return 0;
  }


LOCAL_FUNCTION VOID
g_d_filter_int(ism,jsm,ksm,lsm,ior,jor,kor,lor,ii,jj,kk,ll,pki_int,
           gdb,_scf_info,_irreps)
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
struct mgd *gdb;
scf_struct_t *_scf_info;
scf_irreps_t *_irreps;
{
  int p1,p2;
                                                         
  if (ism == jsm && ksm == lsm && ism == ksm) {
    if (ior == jor && ior == kor || jor == kor && jor == lor) {
      g_d_findit(_scf_info,_irreps,gdb,ii,jj,kk,ll,ism,ksm,pki_int,5);
      }
    else if (ior == kor || jor == lor) {
      g_d_findit(_scf_info,_irreps,gdb,ii,jj,kk,ll,ism,ksm,pki_int,3);

      p1 = MAX0(jj,ll);
      p2 = MIN0(jj,ll);

      g_d_findit(_scf_info,_irreps,gdb,ii,kk,p1,p2,ism,ksm,pki_int,4);
      }
    else if (jor == kor) {
      g_d_findit(_scf_info,_irreps,gdb,ii,jj,kk,ll,ism,ksm,pki_int,3);

      g_d_findit(_scf_info,_irreps,gdb,ii,ll,jj,kk,ism,ksm,pki_int,4);
      }
    else if (ior == jor || kor == lor) {
      g_d_findit(_scf_info,_irreps,gdb,ii,jj,kk,ll,ism,ksm,pki_int,1);

      p1 = MAX0(jj,ll);
      p2 = MIN0(jj,ll);

      g_d_findit(_scf_info,_irreps,gdb,ii,kk,p1,p2,ism,ksm,pki_int,2);
      }
    else {
      g_d_findit(_scf_info,_irreps,gdb,ii,jj,kk,ll,ism,ksm,pki_int,1);

      p1 = MAX0(jj,ll);
      p2 = MIN0(jj,ll);

      g_d_findit(_scf_info,_irreps,gdb,ii,kk,p1,p2,ism,ksm,pki_int,2);

      p1 = MAX0(jj,kk);
      p2 = MIN0(jj,kk);

      g_d_findit(_scf_info,_irreps,gdb,ii,ll,p1,p2,ism,ksm,pki_int,2);
      }
    }
  else if (ism == jsm) {
    g_d_findit(_scf_info,_irreps,gdb,ii,jj,kk,ll,ism,ksm,pki_int,1);
    }
  else if (ism == ksm) {
    if (ior == kor || jor == lor) {
      p1 = MAX0(jj,ll);
      p2 = MIN0(jj,ll);

      g_d_findit(_scf_info,_irreps,gdb,ii,kk,p1,p2,ism,jsm,pki_int,4);
      }
    else {
      p1 = MAX0(jj,ll);
      p2 = MIN0(jj,ll);

      g_d_findit(_scf_info,_irreps,gdb,ii,kk,p1,p2,ism,jsm,pki_int,2);
      }
    }
  }

LOCAL_FUNCTION VOID
g_d_findit(_scf_info,_irreps,gdb,ii,jj,kk,ll,ism,ksm,value,iab)
scf_struct_t *_scf_info;
scf_irreps_t *_irreps;
struct mgd *gdb;
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
  int p1,p2,p3;
  int noi,nok;
  int lij,lkl;
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

  lij = ioff(ii)+jj;
  lkl = ioff(kk)+ll;

  value = (lij == lkl) ? value*0.5 : value;

  switch(iab) {
  case 1:
    gdb->gtmp[lij] += gdb->ptmp[lkl]*value;
    gdb->gtmp[lkl] += gdb->ptmp[lij]*value;
    if(noi && nok) {
      if(_scf_info->special) {
        p1 = MAX0(irr_i->os_num,irr_k->os_num);
        p2 = MIN0(irr_i->os_num,irr_k->os_num);
        p3 = ioff(p1)+p2;
        gdb->gtmpo[lij] += _scf_info->alpha[p3]*value*gdb->ptmpo[lkl];
        gdb->gtmpo[lkl] += _scf_info->alpha[p3]*value*gdb->ptmpo[lij];
        }
      }
    break;
  case 2:
    gdb->gtmp[lij] -= 0.25*gdb->ptmp[lkl]*value;
    gdb->gtmp[lkl] -= 0.25*gdb->ptmp[lij]*value;
    if(noi && nok) {
      if(_scf_info->hsos) {
        gdb->gtmpo[lij] += 0.25*gdb->ptmpo[lkl]*value;
        gdb->gtmpo[lkl] += 0.25*gdb->ptmpo[lij]*value;
        }
      else if(_scf_info->singlet) {
        if (ism != ksm) {
          gdb->gtmpo[lij] -= 0.75*gdb->ptmpo[lkl]*value;
          gdb->gtmpo[lkl] -= 0.75*gdb->ptmpo[lij]*value;
          }
        else {
          gdb->gtmpo[lij] += 0.25*gdb->ptmpo[lkl]*value;
          gdb->gtmpo[lkl] += 0.25*gdb->ptmpo[lij]*value;
          }
        }
      else if(_scf_info->twocon) {
        if (ism != ksm) {
          gdb->gtmpo[lij] += 0.25*gdb->ptmpo[lkl]*value;
          gdb->gtmpo[lkl] += 0.25*gdb->ptmpo[lij]*value;
          }
        }
      else {
        p1 = MAX0(irr_i->os_num,irr_k->os_num);
        p2 = MIN0(irr_i->os_num,irr_k->os_num);
        p3 = ioff(p1)+p2;
        gdb->gtmpo[lij] += 0.25*gdb->ptmpo[lkl]*value*_scf_info->beta[p3];
        gdb->gtmpo[lkl] += 0.25*gdb->ptmpo[lij]*value*_scf_info->beta[p3];
        }
      }
    break;
  case 3:
    gdb->gtmp[lij] += 0.75*gdb->ptmp[lkl]*value;
    gdb->gtmp[lkl] += 0.75*gdb->ptmp[lij]*value;
    if(noi && nok) {
      if(_scf_info->hsos) {
        gdb->gtmpo[lij] += 0.25*gdb->ptmpo[lkl]*value;
        gdb->gtmpo[lkl] += 0.25*gdb->ptmpo[lij]*value;
        }
      else if(_scf_info->singlet) {
        if (ism != ksm) {
          gdb->gtmpo[lij] -= 0.75*gdb->ptmpo[lkl]*value;
          gdb->gtmpo[lkl] -= 0.75*gdb->ptmpo[lij]*value;
          }
        else {
          gdb->gtmpo[lij] += 0.25*gdb->ptmpo[lkl]*value;
          gdb->gtmpo[lkl] += 0.25*gdb->ptmpo[lij]*value;
          }
        }
      else if(_scf_info->twocon) {
        if (ism != ksm) {
          gdb->gtmpo[lij] += 0.25*gdb->ptmpo[lkl]*value;
          gdb->gtmpo[lkl] += 0.25*gdb->ptmpo[lij]*value;
          }
        }
      else {
        p1 = MAX0(irr_i->os_num,irr_k->os_num);
        p2 = MIN0(irr_i->os_num,irr_k->os_num);
        p3 = ioff(p1)+p2;
        gdb->gtmpo[lij] += (_scf_info->alpha[p3] + 0.25*_scf_info->beta[p3])*
                                   value*gdb->ptmpo[lkl];
        gdb->gtmpo[lkl] += (_scf_info->alpha[p3] + 0.25*_scf_info->beta[p3])*
                                   value*gdb->ptmpo[lij];
        }
      }
    break;
  case 4:
    gdb->gtmp[lij] -= 0.5*gdb->ptmp[lkl]*value;
    gdb->gtmp[lkl] -= 0.5*gdb->ptmp[lij]*value;
    if(noi && nok) {
      if(_scf_info->hsos) {
        gdb->gtmpo[lij] += 0.5*gdb->ptmpo[lkl]*value;
        gdb->gtmpo[lkl] += 0.5*gdb->ptmpo[lij]*value;
        }
      else if(_scf_info->singlet) {
        if (ism != ksm) {
          gdb->gtmpo[lij] -= 1.5*gdb->ptmpo[lkl]*value;
          gdb->gtmpo[lkl] -= 1.5*gdb->ptmpo[lij]*value;
          }
        else {
          gdb->gtmpo[lij] += 0.5*gdb->ptmpo[lkl]*value;
          gdb->gtmpo[lkl] += 0.5*gdb->ptmpo[lij]*value;
          }
        }
      else if(_scf_info->twocon) {
        if (ism != ksm) {
          gdb->gtmpo[lij] += 0.5*gdb->ptmpo[lkl]*value;
          gdb->gtmpo[lkl] += 0.5*gdb->ptmpo[lij]*value;
          }
        }
      else {
        p1 = MAX0(irr_i->os_num,irr_k->os_num);
        p2 = MIN0(irr_i->os_num,irr_k->os_num);
        p3 = ioff(p1)+p2;
        gdb->gtmpo[lij] += 0.5*gdb->ptmpo[lkl]*value*_scf_info->beta[p3];
        gdb->gtmpo[lkl] += 0.5*gdb->ptmpo[lij]*value*_scf_info->beta[p3];
        }
      }
    break;
  case 5:
    gdb->gtmp[lij] += 0.5*gdb->ptmp[lkl]*value;
    gdb->gtmp[lkl] += 0.5*gdb->ptmp[lij]*value;
    if(noi && nok) {
      if(_scf_info->hsos) {
        gdb->gtmpo[lij] += 0.5*gdb->ptmpo[lkl]*value;
        gdb->gtmpo[lkl] += 0.5*gdb->ptmpo[lij]*value;
        }
      else if(_scf_info->singlet) {
        if (ism != ksm) {
          gdb->gtmpo[lij] -= 1.5*gdb->ptmpo[lkl]*value;
          gdb->gtmpo[lkl] -= 1.5*gdb->ptmpo[lij]*value;
          }
        else {
          gdb->gtmpo[lij] += 0.5*gdb->ptmpo[lkl]*value;
          gdb->gtmpo[lkl] += 0.5*gdb->ptmpo[lij]*value;
          }
        }
      else if(_scf_info->twocon) {
        if (ism != ksm) {
          gdb->gtmpo[lij] += 0.5*gdb->ptmpo[lkl]*value;
          gdb->gtmpo[lkl] += 0.5*gdb->ptmpo[lij]*value;
          }
        }
      else {
        p1 = MAX0(irr_i->os_num,irr_k->os_num);
        p2 = MIN0(irr_i->os_num,irr_k->os_num);
        p3 = ioff(p1)+p2;
        gdb->gtmpo[lij] += (_scf_info->alpha[p3] + 0.5*_scf_info->beta[p3])*
                                   value*gdb->ptmpo[lkl];
        gdb->gtmpo[lkl] += (_scf_info->alpha[p3] + 0.5*_scf_info->beta[p3])*
                                   value*gdb->ptmpo[lij];
        }
      }
    }
  }

LOCAL_FUNCTION double
max_den(_centers,ptmp,s1,s2,s3,s4)
centers_t *_centers;
double *ptmp;
int s1;
int s2;
int s3;
int s4;
{
  double max=0.0;
  int i,j;
  int p1,p2,p3;
  int s1i =_centers->func_num[s1];
  int s1u = s1i+_centers->center[_centers->center_num[s1]].basis.shell[
                               _centers->shell_num[s1]].nfunc - 1;
  int s2i =_centers->func_num[s2];
  int s2u = s2i+_centers->center[_centers->center_num[s2]].basis.shell[
                                          _centers->shell_num[s2]].nfunc - 1;
  int s3i =_centers->func_num[s3];
  int s3u = s3i+_centers->center[_centers->center_num[s3]].basis.shell[
                                          _centers->shell_num[s3]].nfunc - 1;
  int s4i =_centers->func_num[s4];
  int s4u = s4i+_centers->center[_centers->center_num[s4]].basis.shell[
                                          _centers->shell_num[s4]].nfunc - 1;

  for(i=s1i; i <= s1u ; i++) {
    for(j=s2i; j <= ((s1i==s2i)?i:s2u) ; j++) {
      p1=MAX0(i,j);
      p2=MIN0(i,j);
      p3=ioff(p1)+p2;
      max = MAX0(max,(fabs(ptmp[p3])));
      }
    if(s3!=s2) {
      for(j=s3i; j <= ((s1i==s3i)?i:s3u) ; j++) {
        p1=MAX0(i,j);
        p2=MIN0(i,j);
        p3=ioff(p1)+p2;
        max = MAX0(max,(0.25*fabs((ptmp[p3]))));
        }
      }
    if(s4!=s2) {
      for(j=s4i; j <= ((s1i==s4i)?i:s4u) ; j++) {
        p1=MAX0(i,j);
        p2=MIN0(i,j);
        p3=ioff(p1)+p2;
        max = MAX0(max,(0.25*fabs((ptmp[p3]))));
        }
      }
    }

  if(s3!=s1) {
    for(i=s3i; i <= s3u ; i++) {
      for(j=s4i; j <= ((s3i==s4i)?i:s4u) ; j++) {
        p1=MAX0(i,j);
        p2=MIN0(i,j);
        p3=ioff(p1)+p2;
        max = MAX0(max,(fabs(ptmp[p3])));
        }
      }
    }

  if(s2!=s1) {
    for(i=s2i; i <= s2u ; i++) {
      for(j=s3i; j <= ((s2i==s3i)?i:s3u) ; j++) {
        p1=MAX0(i,j);
        p2=MIN0(i,j);
        p3=ioff(p1)+p2;
        max = MAX0(max,(0.25*fabs(ptmp[p3])));
        }
      if(s4!=s3) {
        for(j=s4i; j <= ((s2i==s4i)?i:s4u) ; j++) {
          p1=MAX0(i,j);
          p2=MIN0(i,j);
          p3=ioff(p1)+p2;
          max = MAX0(max,(0.25*fabs(ptmp[p3])));
          }
        }
      }
    }

  return(max);
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
mgd_int_loop(_centers,_scf_info,_irreps,gdb,_mgdbuff,s1,s2,s3,s4,qijkl,tnint)
centers_t *_centers;
scf_struct_t *_scf_info;
scf_irreps_t *_irreps;
struct mgd *gdb;
double *_mgdbuff;
int s1;
int s2;
int s3;
int s4;
int qijkl;
int *tnint;
{
  int n1,n2,n3,n4;
  int e12,e34,e13e24;
  int bf1,bf2,bf3,bf4;
  int i1,j1,k1,l1;
  int i2,j2,k2,l2;
  int ii,jj,kk,ll;
  int ism,jsm,ksm,lsm;
  int ior,jor,kor,lor;
  int index;
  double pki_int;

  ism=jsm=ksm=lsm=0;

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
          if (INT_NONZERO(_mgdbuff[index])) {
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
        
            pki_int = (double) qijkl*_mgdbuff[index];

            g_d_filter_int(ism,jsm,ksm,lsm,ior,jor,kor,lor,
                     ii,jj,kk,ll,pki_int,gdb,_scf_info,_irreps);
            (*tnint)++;
            }
          index++;
          }
        }
      }
    }   
  }
