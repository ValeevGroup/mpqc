
/* $Log$
 * Revision 1.1  1993/12/29 12:53:13  etseidl
 * Initial revision
 *
 * Revision 1.1.1.1  1992/03/17  17:10:41  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  17:10:40  seidl
 * Initial revision
 *
 * Revision 1.1  1991/12/20  15:45:23  seidl
 * Initial revision
 *
 * Revision 1.3  1991/12/03  16:40:34  etseidl
 * global and local become gbl and lcl
 *
 * Revision 1.2  1991/12/02  19:53:17  seidl
 * *** empty log message ***
 *
 * Revision 1.1  1991/11/22  18:27:53  seidl
 * Initial revision
 * */

static char rcsid[] = "$Id$";

#include <stdio.h>
#include <tmpl.h>
#include <math/array/math_lib.h>

#include "symm.h"
#include "symminit.h"

#include "alloc.gbl"
#include "alloc.lcl"

#define NAMES_LENGTH 256

GLOBAL_FUNCTION int
alloc_s_d_vector(_vector, _irreps)
sym_d_vector_t *_vector;
scf_irreps_t *_irreps;
{
  int i;
  int errcod;

  init_sym_d_vector(_vector);

  _vector->nirrep = _irreps->nirrep;
  _vector->ir = 
      (double_vector_t *) malloc(sizeof(double_vector_t)*_irreps->nirrep);
  if(_vector->ir == NULL) return(-1);

  for(i=0; i < _irreps->nirrep ; i++) {
    errcod = allocbn_double_vector(&_vector->ir[i],"n",
                                               _irreps->ir[i].num_so_tri);
    if(errcod != 0) return(-1);
    }
  return(0);
  }

GLOBAL_FUNCTION int
alloc_s_d_matrix(_matrix, _irreps)
sym_d_matrix_t *_matrix;
scf_irreps_t *_irreps;
{
  int i;
  int errcod;

  init_sym_d_matrix(_matrix);

  _matrix->nirrep = _irreps->nirrep;
  _matrix->ir = 
        (double_matrix_t *) malloc(sizeof(double_matrix_t)*_irreps->nirrep);
  if(_matrix->ir == NULL) return(-1);

  for(i=0; i < _irreps->nirrep ; i++) {
    errcod = allocbn_double_matrix(&_matrix->ir[i],"n1 n2",
                               _irreps->ir[i].num_so,_irreps->ir[i].num_so);
    if(errcod != 0) return(-1);
    }
  return(0);
  }
