
/* This function calculates the closed and open-shell density matrices
 * P(i,j) = sum(k=1,ndocc) 2.0*C(i,k)*C(j,k)
 * Po(i,j) = sum(k=ndocc+1,nsocc) C(i,k)*C(j,k)
 * 
 * Because we are only storing the lower triangle, the off-diagonal elements
 * are multiplied by 2
 */

#include <stdio.h>
#include <assert.h>
#include <tmpl.h>
#include <util/misc/libmisc.h>
#include <math/dmt/libdmt.h>
#include <math/array/math_lib.h>
#include <chemistry/qc/intv2/int_libv2.h>

#include <chemistry/qc/dmtscf/scf.h>

#include <chemistry/qc/dmtscf/scf_mkden.gbl>
#include <chemistry/qc/dmtscf/scf_mkden.lcl>

GLOBAL_FUNCTION int
scf_make_density(scf_info,Scf_Vec,Pmat,DPmat,PmatO,DPmatO,occnum)
scf_struct_t *scf_info;
dmt_matrix Scf_Vec;
dmt_matrix Pmat;
dmt_matrix DPmat;
dmt_matrix PmatO;
dmt_matrix DPmatO;
double *occnum;
{
  int nlb;
  int nlocalb;
  int ib,jb,i,j,ist,jst,isz,jsz;
  int k;
  int ndocc,nsocc;
  double *col,*pblk;
  loop_t *loop;

  assert(dmt_distribution(Scf_Vec) == COLUMNS);
  assert(dmt_distribution(Pmat) == SCATTERED);
  assert(dmt_distribution(DPmat) == SCATTERED);
  if (scf_info->iopen) {
    assert(dmt_distribution(PmatO) == SCATTERED);
    assert(dmt_distribution(DPmatO) == SCATTERED);
  }

  ndocc = scf_info->nclosed;
  nsocc = scf_info->nopen;

  nlocalb = dmt_nlocal(Pmat);

/* save old density in DPmat, and scale by -1.0 */

  dmt_copy(Pmat,DPmat);
  dmt_scale(DPmat,-1.0);
  dmt_fill(Pmat,0.0);

  if (scf_info->iopen) {
    dmt_copy(PmatO,DPmatO);
    dmt_scale(DPmatO,-1.0);
    dmt_fill(PmatO,0.0);
  }

 /* now put scfvec on the loop */

  loop = dmt_ngl_create("%mr",Scf_Vec);
  while (dmt_ngl_next(loop)) {
    int iind,isize,jsize;
    dmt_ngl_create_inner(loop,0);
    while (dmt_ngl_next_inner_m(loop,&iind,&isize,&k,&jsize,&col)) {

      for (nlb=0; nlb < nlocalb ; nlb++) {
        if (k < ndocc) {
          dmt_get_block(Pmat,nlb,&ib,&jb,&pblk);
          dmt_describe_block(Pmat,ib,&ist,&isz);
          dmt_describe_block(Pmat,jb,&jst,&jsz);

          for (i=0; i < isz ; i++) {
            for (j=0; j < jsz ; j++) {
              pblk[i*jsz+j] += col[ist+i]*col[jst+j];
            }
          }
        } else if (k < ndocc+nsocc) {
          dmt_get_block(PmatO,nlb,&ib,&jb,&pblk);
          dmt_describe_block(PmatO,ib,&ist,&isz);
          dmt_describe_block(PmatO,jb,&jst,&jsz);

          for (i=0; i < isz ; i++) {
            for (j=0; j < jsz ; j++) {
              pblk[i*jsz+j] += occnum[k]*col[ist+i]*col[jst+j];
            }
          }
        }
      }
    }
  }

 /* now kill loop */

  dmt_ngl_kill(loop);
  
 /* now scale Pmat */

  dmt_scale(Pmat,4.0);
  dmt_scale_diagonal(Pmat,0.5);

  if (scf_info->iopen) {
    dmt_scale(PmatO,2.0);
    dmt_scale_diagonal(PmatO,0.5);
    dmt_sum(PmatO,DPmatO);
    dmt_sum(PmatO,Pmat);
  }

  dmt_sum(Pmat,DPmat);
  return 0;
}
