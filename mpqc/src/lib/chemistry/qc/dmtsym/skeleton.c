
/* $Log$
 * Revision 1.1  1993/12/29 12:52:58  etseidl
 * Initial revision
 *
 * Revision 1.1.1.1  1992/03/17  16:27:33  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  16:27:32  seidl
 * Initial revision
 *
 * Revision 1.3  1992/01/27  16:49:52  seidl
 * libint needs libmath
 *
 * Revision 1.2  1992/01/27  16:41:51  seidl
 * remove some unnecessary includes
 *
 * Revision 1.1  1992/01/27  12:53:04  seidl
 * Initial revision
 *
 * Revision 1.2  1992/01/21  11:32:49  seidl
 * use libintv2
 *
 * Revision 1.1  1991/12/20  15:45:23  seidl
 * Initial revision
 *
 * Revision 1.2  1991/12/03  16:42:53  etseidl
 * global and local changed to gbl and lcl
 *
 * Revision 1.1  1991/11/22  18:28:21  seidl
 * Initial revision
 * */

static char rcsid[] = "$Id$";

#include <stdio.h>
#include <math.h>
#include <tmpl.h>
#include <math/array/math_lib.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <math/dmt/libdmt.h>

#include "symm.h"
#include "symerr.gbl"

#include "skeleton.gbl"
#include "skeleton.lcl"

GLOBAL_FUNCTION VOID
sym_skel_matrix(_sym_info,mat,_outfile)
sym_struct_t *_sym_info;
dmt_matrix mat;
FILE *_outfile;
{
  int i,j;
  int si,sj,sij;
  int isiz,jsiz;
  int nlocal;
  double lij;
  double *tmp;

  if(dmt_distribution(mat)!=SCATTERED) {
    fprintf(_outfile,"sym_skel_matrix: want scattered matrices only\n");
    exit(-1);
    }

  nlocal = dmt_nlocal(mat);

  for(i=0; i < nlocal ; i++) {
    dmt_get_block_dsc(mat,i,&si,&isiz,&sj,&jsiz,&tmp);
    sij = (si > sj) ? si*(si+1)/2 + sj : sj*(sj+1)/2 +si;
    lij = (double) _sym_info->lamij[sij];
    for(j=0; j < isiz*jsiz ; j++) tmp[j] *= lij;
    }

  }
