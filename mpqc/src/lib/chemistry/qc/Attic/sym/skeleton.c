
/* $Log$
 * Revision 1.1  1993/12/29 12:53:13  etseidl
 * Initial revision
 *
 * Revision 1.1.1.1  1992/03/17  17:11:00  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  17:10:59  seidl
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
#include <util/ipv2/ip_libv2.h>
#include <chemistry/qc/intv2/int_libv2.h>

#include "symm.h"
#include "symerr.gbl"

#include "skeleton.gbl"
#include "skeleton.lcl"

GLOBAL_FUNCTION VOID
sym_skel_matrix(_centers,_sym_info,mat,_outfile)
centers_t *_centers;
sym_struct_t *_sym_info;
double_matrix_t *mat;
FILE *_outfile;
{
  int i,j,ij;
  int si,sj,sij;
  int fi,fj;
  double lij;
  double **tmp=mat->d;

  for(si=sij=0; si < _centers->nshell ; si++) {
    fi = _centers->func_num[si];
    for(sj=0; sj <= si ; sj++,sij++) {
      fj = _centers->func_num[sj];
      lij = _sym_info->lamij[sij];
      for(i=fi; i < fi + INT_SH_NFUNC(_centers,si); i++) {
        for(j=fj; j < fj + INT_SH_NFUNC(_centers,sj); j++) {
          tmp[i][j] *= lij;
          if(si!=sj) tmp[j][i] = tmp[i][j];
          }
        }
      }
    }
  }

GLOBAL_FUNCTION VOID
sym_skel_vector(_centers,_sym_info,vec,_outfile)
centers_t *_centers;
sym_struct_t *_sym_info;
double_vector_t *vec;
FILE *_outfile;
{
  int i,j,ij;
  int si,sj,sij;
  int fi,fj;
  double lij;
  double *tmp=vec->d;

  ij=0;
  for(si=sij=0; si < _centers->nshell ; si++) {
    fi = _centers->func_num[si];
    for(sj=0; sj <= si ; sj++,sij++) {
      fj = _centers->func_num[sj];
      lij = _sym_info->lamij[sij];
      for(i=fi; i < fi + INT_SH_NFUNC(_centers,si); i++) {
        for(j=fj; j <((si==sj)? i+1:fj+INT_SH_NFUNC(_centers,sj)); j++,ij++) {
          tmp[ij] *= lij;
          }
        }
      }
    }
  }
