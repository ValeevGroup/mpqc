
/* $Log$
 * Revision 1.1  1993/12/29 12:53:15  etseidl
 * Initial revision
 *
 * Revision 1.12  1992/06/23  20:04:30  seidl
 * change dmt matrices to uppercase,
 * get rid of unnecessary matrice
 *
 * Revision 1.11  1992/06/17  21:53:13  jannsen
 * cleaned up for saber-c and changed to ngl loops
 *
 * Revision 1.10  1992/06/16  16:25:08  seidl
 * add tim_print(), remove call to scf_make_g_d_o
 *
 * Revision 1.9  1992/05/26  20:17:45  jannsen
 * use mtype_get to get message types for global operations
 * check results of memory allocations
 *
 * Revision 1.8  1992/05/19  20:56:36  seidl
 * use message types 7000-7999
 *
 * Revision 1.7  1992/05/04  11:03:29  seidl
 * use gdcomb for getting local matrix off loop, remove pk integral stuff
 *
 * Revision 1.6  1992/04/22  15:54:05  seidl
 * add timing call
 *
 * Revision 1.5  1992/04/09  17:55:08  seidl
 * use signed char for maxp array
 *
 * Revision 1.4  1992/04/06  12:35:46  seidl
 * include math.h
 *
 * Revision 1.3  1992/04/01  01:02:47  seidl
 * fixed new bounds
 *
 * Revision 1.2  1992/03/21  00:38:38  seidl
 * change sym_libv2.h to chemistry/qc/dmtsym/sym_dmt.h
 *
 * Revision 1.1.1.1  1992/03/17  16:26:15  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  16:26:14  seidl
 * Initial revision
 *
 * Revision 1.8  1992/03/09  13:00:16  seidl
 * for local_P make maxp array
 *
 * Revision 1.7  1992/03/04  15:55:47  seidl
 * add scf bounds checking
 *
 * Revision 1.6  1992/02/28  19:00:40  seidl
 * add call to scf_make_g_l
 *
 * Revision 1.5  1992/02/21  15:10:19  seidl
 * start adding support for non-local density matrices
 *
 * Revision 1.4  1992/02/10  16:59:33  seidl
 * don't call scf_make_g if no integrals written to disk
 *
 * Revision 1.3  1992/02/07  12:58:44  seidl
 * remove timing stuff
 *
 * Revision 1.2  1992/02/05  14:11:35  seidl
 * use util/misc/libmisc timing routines
 *
 * Revision 1.1  1992/02/04  23:48:08  seidl
 * Initial revision
 *
 * Revision 1.3  1992/01/16  19:53:07  seidl
 * use new integral routines
 *
 * Revision 1.2  1992/01/14  19:33:40  seidl
 * change a few print statements
 *
 * Revision 1.1  1992/01/13  19:13:15  seidl
 * Initial revision
 * */

static char rcsid[] = "$Id$";

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <tmpl.h>
#include <comm/picl/picl.h>
#include <comm/picl/ext/piclext.h>
#include <math/dmt/libdmt.h>
#include <util/misc/libmisc.h>
#include <math/array/math_lib.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <chemistry/qc/dmtsym/sym_dmt.h>
#include "scf.h"

#ifndef IOFF
#define IOFF(a,b) ((a)>(b))?(a)*((a)+1)/2+(b):(b)*((b)+1)/2+(a)
#endif

#include "scf_mkgd.gbl"
#include "scf_loopg.gbl"
#include "scf_bnd.gbl"

#include "scf_gmat.gbl"
#include "scf_gmat.lcl"

static double *intbuf=NULL;

GLOBAL_FUNCTION int
scf_make_gmat(_scf_info,_sym_info,_irreps,_centers,
       GMAT,GMATO,DPMAT,DPMATO,SSCR1,SSCR2,_init,_iter,_outfile)
scf_struct_t *_scf_info;
sym_struct_t *_sym_info;
scf_irreps_t *_irreps;
centers_t *_centers;
dmt_matrix GMAT;
dmt_matrix GMATO;
dmt_matrix DPMAT;
dmt_matrix DPMATO;
dmt_matrix SSCR1;
dmt_matrix SSCR2;
int _init;
int _iter;
FILE *_outfile;
{
  int iopen=_scf_info->iopen;
  int errcod;
  int ntri=_scf_info->nbatri;
  int flags=0;
  signed char *maxp;

  double_vector_t ptmp,ptmpo,gtmp,gtmpo;

  if((_scf_info->print_flg & 16) && _init) tim_print(0);

  tim_enter("gmat");
  tim_enter("init");

/* if using local density matrices, then allocate memory for them */

  if(_scf_info->local_p) {
    errcod = allocbn_double_vector(&gtmp,"n",ntri);
    if(errcod != 0) {
      fprintf(_outfile,"scf_make_gmat: malloc trouble 1\n");
      return(-1);
      }
    zero_double_vector(&gtmp);

    errcod = allocbn_double_vector(&ptmp,"n",ntri);
    if(errcod != 0) {
      fprintf(_outfile,"scf_make_gmat: malloc trouble 2\n");
      return(-1);
      }
    zero_double_vector(&ptmp);

    if(_scf_info->eliminate) {
      int nshell=dmt_nblocks(DPMAT);
      int nsht=nshell*(nshell+1)/2;
      maxp = (signed char *) malloc(nsht);
      if(maxp==NULL) {
        fprintf(_outfile,"trouble mallocing maxp\n");
        return(-1);
        }
      }
    scat_to_dv(&ptmp,DPMAT,maxp,_scf_info->eliminate);

    if(iopen) {
      errcod = allocbn_double_vector(&gtmpo,"n",ntri);
      if(errcod != 0) {
        fprintf(_outfile,"scf_make_gmat: malloc trouble 3\n");
        return(-1);
        }
      zero_double_vector(&gtmpo);

      errcod = allocbn_double_vector(&ptmpo,"n",ntri);
      if(errcod != 0) {
        fprintf(_outfile,"scf_make_gmat: malloc trouble 4\n");
        return(-1);
        }
      zero_double_vector(&ptmpo);
      scat_to_dv(&ptmpo,DPMATO,NULL,0);
      }
    }

 /* if this is the first time through, initialize the buffer which will
  * contain the integrals
  */
  if(_init) {
    int_initialize_offsets2(_centers,_centers,_centers,_centers);

    flags = INT_EREP|INT_NOSTRB|INT_NOSTR1|INT_NOSTR2;
#if 0
    if(!_scf_info->local_p) flags |= (INT_NOSTR1|INT_NOSTR2);
#endif

    intbuf = 
      int_initialize_erep(flags,0,_centers,_centers,_centers,_centers);

    int_storage(_scf_info->int_store);
    if(_scf_info->eliminate || _scf_info->local_p) 
      scf_init_bounds(_centers,intbuf);

    }
  tim_exit("init");
  if((_scf_info->print_flg & 16) && _init) tim_print(0);

  if(_scf_info->local_p) {

  /* calculate integrals directly and stuff into the g matrix */
    errcod = scf_make_g_d(_centers,_irreps,_scf_info,_sym_info,
               &gtmp,&gtmpo,&ptmp,&ptmpo,maxp,intbuf,_iter,_outfile);
    if(errcod != 0) {
      fprintf(_outfile,"scf_iter: trouble forming gmat 2\n");
      return(-1);
      }
    }

/* if not using local density, then form the g matrix loopwise */
  else {
    errcod = scf_make_g_l(_centers,_scf_info,_sym_info,
               GMAT,GMATO,DPMAT,DPMATO,SSCR1,SSCR2,intbuf,_iter,_outfile);
    if(errcod != 0) {
      fprintf(_outfile,"scf_iter: trouble forming gmat 3\n");
      return(-1);
      }
    }

  tim_enter("clean");
  int_reduce_storage_threshold();

  if(_scf_info->local_p) {
    dv_to_scat(&gtmp,GMAT,&ptmp);
    if(iopen) dv_to_scat(&gtmpo,GMATO,&ptmpo);
  
    free_double_vector(&ptmp);
    free_double_vector(&gtmp);
    if(iopen) {
      free_double_vector(&ptmpo);
      free_double_vector(&gtmpo);
      }
    if(_scf_info->eliminate) free(maxp);
    }

  tim_exit("clean");
  tim_exit("gmat");

  return 0;
  }

LOCAL_FUNCTION VOID
scat_to_dv(dv,sm,maxp,elim)
double_vector_t *dv;
dmt_matrix sm;
signed char *maxp;
int elim;
{
  int i,j,ij,ib,jb,isz,jsz,ist,jst,lij;
  double *blk;
  double linv=1.0/log(2.0);
  double tol=pow(2.0,-126.0);
  double tmp,ftmp;
  loop_t *loop;

  loop = dmt_ngl_create("%mr",sm);

  while(dmt_ngl_next(loop)) {
    dmt_ngl_create_inner(loop,0);
    while(dmt_ngl_next_inner_m(loop,&ib,&isz,&jb,&jsz,&blk)) {
      dmt_describe_block(sm,ib,&ist,&isz);
      dmt_describe_block(sm,jb,&jst,&jsz);

      if(ib!=jb) {
        for(i=0; i < isz ; i++) {
          lij=IOFF((i+ist),jst);
          for(j=0; j < jsz ; j++,lij++) {
            dv->d[lij] = blk[i*jsz+j];
            }
          }
        }
      else {
        for(i=0; i < isz ; i++) {
          lij=IOFF((i+ist),jst);
          for(j=0; j <= i ; j++,lij++) {
            dv->d[lij] = blk[i*jsz+j];
            }
          }
        }

      if(elim) {
        tmp=0.0;
        for(i=0; i < isz*jsz; i++) 
          if ((ftmp=fabs(blk[i])) > tmp) tmp=ftmp;

        tmp=(tmp>tol) ? tmp:tol;

        ij=ib*(ib+1)/2+jb;
        maxp[ij]=(signed char) (log(tmp)*linv);
        }

      }
    }

  dmt_ngl_kill(loop);
  }


LOCAL_FUNCTION VOID
dv_to_scat(dv,sm,scr)
double_vector_t *dv;
dmt_matrix sm;
double_vector_t *scr;
{
  int i,j,ib,jb,isz,jsz,ist,jst,lij;
  int nlocal,nl;
  int dim = cubedim0();
  double *lblk;

  nlocal = dmt_nlocal(sm);

 /* sum up locally held matrices */
  gop1(dv->d,dv->n,scr->d,'+',mtype_get());

 /* and transfer to locally held blocks of distributed matrix */
  for(nl=0; nl < nlocal ; nl++) {
    dmt_get_block(sm,nl,&ib,&jb,&lblk);
    dmt_describe_block(sm,ib,&ist,&isz);
    dmt_describe_block(sm,jb,&jst,&jsz);

    if(ib!=jb) {
      for(i=0; i < isz ; i++) {
        lij=IOFF((i+ist),jst);
        for(j=0; j < jsz ; j++,lij++) {
          lblk[i*jsz+j] += dv->d[lij];
          }
        }
      }
    else {
      for(i=0; i < isz ; i++) {
        lij=IOFF((i+ist),jst);
        for(j=0; j <= i ; j++,lij++) {
          lblk[i*jsz+j] += dv->d[lij];
          if(i!=j) lblk[j*jsz+i] += dv->d[lij];
          }
        }
      }
    }
  }
