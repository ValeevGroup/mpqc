
/* This function calculates the closed and open-shell density matrices
 * P(i,j) = sum(k=1,ndocc) 2.0*C(i,k)*C(j,k)
 * Po(i,j) = sum(k=ndocc+1,nsocc) C(i,k)*C(j,k)
 * 
 * Because we are only storing the lower triangle, the off-diagonal elements
 * are multiplied by 2
 */

/* $Log$
 * Revision 1.1  1993/12/29 12:53:04  etseidl
 * Initial revision
 *
 * Revision 1.1.1.1  1992/03/17  17:08:57  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  17:08:56  seidl
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

#define HOCC 0

#include <stdio.h>
#include <tmpl.h>
#include <math/array/math_lib.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <chemistry/qc/sym/sym_lib.h>

#include "scf.h"

#include "scf_mkden.gbl"
#include "scf_mkden.lcl"

GLOBAL_FUNCTION int
scf_make_density(_scf_info,_irreps,
       _scfvec,_pmat,_dpmat,_pmato,_dpmato,_occnum,_outfile)
scf_struct_t *_scf_info;
scf_irreps_t *_irreps;
sym_d_matrix_t *_scfvec;
sym_d_vector_t *_pmat;
sym_d_vector_t *_dpmat;
sym_d_vector_t *_pmato;
sym_d_vector_t *_dpmato;
sym_d_vector_t *_occnum;
FILE *_outfile;
{
  int i,j,k,l,ij;
  int n;
  int ndocc,nsocc,nhocc;
  double ptempc,ptempo,ctmp;
  double **sv,*p,*dp,*po,*dpo,*occ;
  scf_irrep_t *s;

  for(l=0; l < _irreps->nirrep ; l++) {
    s = &_irreps->ir[l];
    if(n=s->num_so) {
      sv=_scfvec->ir[l].d;
      p=_pmat->ir[l].d;
      dp=_dpmat->ir[l].d;
      occ=_occnum->ir[l].d;
      if(_scf_info->iopen) {
        po=_pmato->ir[l].d;
        dpo=_dpmato->ir[l].d;
        }
      ndocc = s->nclosed;
      nsocc = s->nopen;

      for(i=ij=0; i < n ; i++ ) {
        for(j=0; j < i; j++,ij++) {
          ptempc=ptempo=0.0;
          for (k=0; k < ndocc ; k++)
            ptempc += 4.0*sv[i][k]*sv[j][k];

          for (k=ndocc; k < ndocc+nsocc ; k++)
            ptempo += 
               2.0*occ[k]*sv[i][k]*sv[j][k];

#if HOCC
          for (k=ndocc+nsocc; k < ndocc+nsocc+nhocc ; k++)
            ptempo += 
               2.0*occ[k]*sv[i][k]*sv[j][k];
#endif

          if(_scf_info->iopen) {
            dpo[ij] = ptempo - po[ij];
            po[ij] = ptempo;
            }
          dp[ij] = ptempc+ptempo - p[ij];
          p[ij] = ptempc+ptempo;
          }
        ptempc=ptempo=0.0;
        for (k=0; k < ndocc ; k++) {
          ctmp=sv[i][k];
          ptempc += 2.0*ctmp*ctmp;
          }
        for (k=ndocc; k < ndocc+nsocc ; k++) {
          ctmp=sv[i][k];
          ptempo += ctmp*ctmp*occ[k];
          }
#if HOCC
        for (k=ndocc+nsocc; k < ndocc+nsocc+nhocc ; k++) {
          ctmp=sv[i][k];
          ptempo += ctmp*ctmp*occ[k];
          }
#endif
        if(_scf_info->iopen) {
          dpo[ij] = ptempo - po[ij];
          po[ij] = ptempo;
          }
        dp[ij] = ptempc+ptempo - p[ij];
        p[ij] = ptempc+ptempo;
        ij++;
        }
      }
    }
  return 0;
  }
