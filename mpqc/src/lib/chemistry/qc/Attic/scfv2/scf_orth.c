/* $Log$
 * Revision 1.1  1993/12/29 12:53:04  etseidl
 * Initial revision
 *
 * Revision 1.1.1.1  1992/03/17  17:09:19  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  17:09:17  seidl
 * Initial revision
 *
 * Revision 1.2  1992/01/16  19:53:07  seidl
 * use new integral routines
 *
 * Revision 1.1  1992/01/03  12:44:49  seidl
 * Initial revision
 * */

static char rcsid[] = "$Id$";

#include <stdio.h>
#include <math.h>
#include <tmpl.h>
#include <math/array/math_lib.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <chemistry/qc/sym/sym_lib.h>
#include "scf.h"

#include "scf_orth.gbl"
#include "scf_orth.lcl"

GLOBAL_FUNCTION int
scf_schmit(_scf_info,_irreps,_mfp,_mast,_scf_vec,
           _scrsv1,_scr1,_scrv1,all,_outfile)
scf_struct_t *_scf_info;
scf_irreps_t *_irreps;
scf_mf_ptrs_t *_mfp;
int _mast;
sym_d_matrix_t *_scf_vec;
sym_d_vector_t *_scrsv1;
double_matrix_t *_scr1;
double_vector_t *_scrv1;
int all;
FILE *_outfile;
{
  int i,j,ij,nn;
  int n,m,ncol;
  int *mptr=&_mfp->cur_ptr;
  int errcod;
  double *v=_scrv1->d;
  double **ctmp=_scr1->d;
  double *smat;
  double vtmp;
  scf_irrep_t *s;

/* read overlap integrals from master file */
  free_sym_d_vector(_scrsv1);
  init_sym_d_vector(_scrsv1);
  *mptr=_mfp->overlap;
  errcod = bread_sym_d_vector(_mast,_scrsv1,mptr);
  if(errcod < 0) {
    fprintf(_outfile,"scf_schmit:");
    fprintf(_outfile," trouble reading overlap integrals\n");
    return(-1);
    }

  for(n=0; n < _irreps->nirrep ; n++) {
    s=&_irreps->ir[n];
    smat=_scrsv1->ir[n].d;
    if(nn=s->num_so) {
      for (i=0; i < nn ; i++)
        for (j=0; j < nn ; j++)
          ctmp[j][i] = _scf_vec->ir[n].d[i][j];

      ncol = s->nclosed + s->nopen;
      if(!ncol) ncol++;
      if(all) ncol = nn;
      for(m=0; m < ncol ; m++) {
        v[0]=ctmp[m][0]*smat[0];
        for(i=ij=1; i < nn ; i++) {
          for(j=0,vtmp=0.0; j < i ; j++,ij++) {
            vtmp += ctmp[m][j]*smat[ij];
            v[j] += ctmp[m][i]*smat[ij];
            }
          v[i] = vtmp+ctmp[m][i]*smat[ij];
          ij++;
          }
        for(i=0,vtmp=0.0; i < nn ; i++) vtmp += v[i]*ctmp[m][i];
        if(!vtmp) return(-1);
        if(vtmp < 1.0e-20) vtmp = 1.0e-20;
        vtmp = 1.0/sqrt(vtmp);

        for(i=0; i < nn ; i++) {
          v[i] *= vtmp;
          ctmp[m][i] *= vtmp;
          }

        if(m < ncol-1) {
          for(i=m+1,vtmp=0.0; i < ncol ; i++) {
            for(j=0,vtmp=0.0; j<nn ;j++) vtmp += v[j]*ctmp[i][j];
            for(j=0; j < nn ; j++) ctmp[i][j] -= vtmp*ctmp[m][j];
            }
          }
        }

      for (i=0; i < nn ; i++)
        for (j=0; j < nn ; j++)
          _scf_vec->ir[n].d[i][j] = ctmp[j][i];
      }
    }
  return(0);
  }
