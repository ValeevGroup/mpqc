
#include <stdio.h>
#include <math.h>
#include <tmpl.h>
#include <math/array/math_lib.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <math/dmt/libdmt.h>

#include <chemistry/qc/dmtsym/symm.h>

#include <chemistry/qc/dmtsym/skeleton.gbl>
#include <chemistry/qc/dmtsym/skeleton.lcl>

/************************************************************************
 * 
 * this is really just a debugging function.  given a symmetric dmt matrix
 * and the sym struct, create a skeleton matrix a la Dupuis and King.
 * mat is overwritten with the skeleton matrix
 *
 */

GLOBAL_FUNCTION void
sym_skel_matrix(sym_info,mat)
sym_struct_t *sym_info;
dmt_matrix mat;
{
  int i,j;
  int si,sj,sij;
  int isiz,jsiz;
  int nlocal;
  double lij;
  double *tmp;

  if (dmt_distribution(mat) != SCATTERED) {
    fprintf(stderr,"sym_skel_matrix: want scattered matrices only\n");
    exit(-1);
  }

  nlocal = dmt_nlocal(mat);

  for (i=0; i < nlocal ; i++) {
    dmt_get_block_dsc(mat,i,&si,&isiz,&sj,&jsiz,&tmp);
    sij = (si > sj) ? si*(si+1)/2 + sj : sj*(sj+1)/2 +si;
    lij = (double) sym_info->lamij[sij];
    for (j=0; j < isiz*jsiz ; j++) tmp[j] *= lij;
  }
}
