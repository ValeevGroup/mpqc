
/* This function calculates the closed and open-shell density matrices
 * P(i,j) = sum(k=1,ndocc) 2.0*C(i,k)*C(j,k)
 * Po(i,j) = sum(k=ndocc+1,nsocc) C(i,k)*C(j,k)
 * 
 * Because we are only storing the lower triangle, the off-diagonal elements
 * are multiplied by 2
 */

/* $Log$
 * Revision 1.1  1993/12/29 12:53:16  etseidl
 * Initial revision
 *
 * Revision 1.4  1992/06/23  20:04:41  seidl
 * change dmt matrices to uppercase,
 * get rid of unnecessary matrice
 *
 * Revision 1.3  1992/06/17  21:54:22  jannsen
 * cleaned up for saber-c
 *
 * Revision 1.2  1992/03/21  00:39:48  seidl
 * change sym_libv2.h to chemistry/qc/dmtsym/sym_dmt.h
 *
 * Revision 1.1.1.1  1992/03/17  16:26:33  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  16:26:32  seidl
 * Initial revision
 *
 * Revision 1.2  1992/02/04  23:48:08  seidl
 * use math/dmt/libdmt-like density matrix formation
 *
 * Revision 1.1  1992/01/30  19:30:45  seidl
 * Initial revision
 *
 * Revision 1.4  1992/01/16  19:53:07  seidl
 * use new integral routines
 *
 * Revision 1.3  1992/01/13  19:15:08  seidl
 * add some comments
 *
 * Revision 1.2  1992/01/09  11:45:36  seidl
 * no longer include util/ipv2/ip_libv2
 *
 * Revision 1.1  1991/12/20  16:23:10  seidl
 * Initial revision
 * */

static char *rcsid = "$Id$";

#include <stdio.h>
#include <tmpl.h>
#include <math/dmt/libdmt.h>
#include <math/array/math_lib.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <chemistry/qc/dmtsym/sym_dmt.h>

#include "scf.h"

#include "scf_mkden.gbl"
#include "scf_mkden.lcl"

GLOBAL_FUNCTION int
scf_make_density(_scf_info,_irreps,
       SCF_VEC,PMAT,DPMAT,PMATO,DPMATO,_occnum,_outfile)
scf_struct_t *_scf_info;
scf_irreps_t *_irreps;
dmt_matrix SCF_VEC;
dmt_matrix PMAT;
dmt_matrix DPMAT;
dmt_matrix PMATO;
dmt_matrix DPMATO;
double_vector_t *_occnum;
FILE *_outfile;
{
  int nlb;
  int nlocalb;
  int ib,jb,i,j,ist,jst,isz,jsz;
  int k;
  int ndocc,nsocc;
  double *col,*occ,*pblk;
  scf_irrep_t *s;
  loop_t *loop;

  s = &_irreps->ir[0];
  occ = _occnum->d;
  ndocc = s->nclosed;
  nsocc = s->nopen;

  nlocalb = dmt_nlocal(PMAT);

/* save old density in DPMAT, and scale by -1.0 */

  dmt_copy(PMAT,DPMAT);
  dmt_scale(DPMAT,-1.0);
  dmt_fill(PMAT,0.0);

  if(_scf_info->iopen) {
    dmt_copy(PMATO,DPMATO);
    dmt_scale(DPMATO,-1.0);
    dmt_fill(PMATO,0.0);
    }

 /* now put scfvec on the loop */

  loop = dmt_ngl_create("%mr",SCF_VEC);
  while(dmt_ngl_next(loop)) {
    int iind,isize,jsize;
    dmt_ngl_create_inner(loop,0);
    while(dmt_ngl_next_inner_m(loop,&iind,&isize,&k,&jsize,&col)) {

      for(nlb=0; nlb < nlocalb ; nlb++) {
        if(k < ndocc) {
          dmt_get_block(PMAT,nlb,&ib,&jb,&pblk);
          dmt_describe_block(PMAT,ib,&ist,&isz);
          dmt_describe_block(PMAT,jb,&jst,&jsz);

          for(i=0; i < isz ; i++) {
            for(j=0; j < jsz ; j++) {
              pblk[i*jsz+j] += col[ist+i]*col[jst+j];
              }
            }
          }
        else if(k < ndocc+nsocc) {
          dmt_get_block(PMATO,nlb,&ib,&jb,&pblk);
          dmt_describe_block(PMATO,ib,&ist,&isz);
          dmt_describe_block(PMATO,jb,&jst,&jsz);

          for(i=0; i < isz ; i++) {
            for(j=0; j < jsz ; j++) {
              pblk[i*jsz+j] += occ[k]*col[ist+i]*col[jst+j];
              }
            }
          }
        }
      }
    }

 /* now kill loop */

  dmt_ngl_kill(loop);
  
 /* now scale PMAT */

  dmt_scale(PMAT,4.0);
  dmt_scale_diagonal(PMAT,0.5);

  if(_scf_info->iopen) {
    dmt_scale(PMATO,2.0);
    dmt_scale_diagonal(PMATO,0.5);
    dmt_sum(PMATO,DPMATO);
    dmt_sum(PMATO,PMAT);
    }

  dmt_sum(PMAT,DPMAT);
  return 0;
  }
