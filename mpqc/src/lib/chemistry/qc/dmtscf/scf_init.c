
/* $Log$
 * Revision 1.1  1993/12/29 12:53:15  etseidl
 * Initial revision
 *
 * Revision 1.5  1992/06/17  21:54:13  jannsen
 * cleaned up for saber-c
 *
 * Revision 1.4  1992/05/12  10:33:04  seidl
 * no longer calculate nuclear rep. energy here
 *
 * Revision 1.3  1992/05/04  10:58:53  seidl
 * allow coordinate input in angstroms
 *
 * Revision 1.2  1992/03/21  00:38:54  seidl
 * change sym_libv2.h to chemistry/qc/dmtsym/sym_dmt.h
 *
 * Revision 1.1.1.1  1992/03/17  16:26:18  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  16:26:17  seidl
 * Initial revision
 *
 * Revision 1.1  1992/02/04  23:48:08  seidl
 * Initial revision
 *
 * Revision 1.4  1992/01/16  19:53:07  seidl
 * use new integral routines
 *
 * Revision 1.3  1992/01/13  19:13:26  seidl
 * remove some include statements, add some comments
 *
 * Revision 1.2  1991/12/20  16:29:24  seidl
 * many changes to allow use of Dupuis&King symmetry, but without
 * transforming to so basis
 *
 * Revision 1.1  1991/12/17  21:43:01  seidl
 * Initial revision
 *
 * Revision 1.2  1991/12/02  19:58:14  seidl
 * *** empty log message ***
 *
 * Revision 1.1  1991/11/26  19:08:23  seidl
 * Initial revision
 * */

static char rcsid[] = "$Id$";

#include <stdio.h>
#include <tmpl.h>
#include <math/array/math_lib.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <chemistry/qc/dmtsym/sym_dmt.h>
#include "scf.h"

#include "scf_init.gbl"
#include "scf_init.lcl"

#define CLOSED 0
#define OPENSH 1
#define TWOCON 2

GLOBAL_FUNCTION int
scf_initialize(_outfile,_centers,_scf_info,_sym_info,_irreps,pg,angs)
FILE *_outfile;
centers_t *_centers;
scf_struct_t *_scf_info;
sym_struct_t *_sym_info;
scf_irreps_t *_irreps;
char *pg;
int angs;
{
  int i;
  int errcod;
  int nbasis,ntri;

/* initialize the centers structure */

/* first read in information about atoms, basis set, etc, and place
 * in the centers struct.  note, sym_init_centers calls
 * int_normalize_centers()
 */

  errcod = sym_init_centers(_centers,_sym_info,pg,_outfile);
  if(errcod != 0) {
    fprintf(_outfile,"scf_initialize: trouble symmetrizing centers\n");
    return(-1);
    }

/* if the input coordinates are in angstroms, then convert to a.u. */

  if(angs) {
    double atoau=1.0/0.52917706;

    for(i=0; i < _centers->n; i++) {
      _centers->center[i].r[0] *= atoau;
      _centers->center[i].r[1] *= atoau;
      _centers->center[i].r[2] *= atoau;
      }
    }

  int_initialize_1e(0,0,_centers,_centers);
  int_initialize_offsets1(_centers,_centers);

  nbasis = _centers->nfunc;
  ntri = nbasis*(nbasis+1)/2;

 /* now form the _irreps struct */
 /* if not using symmetry, then set up _irreps in C1 symmetry */
 /* for dmt version, don't use so basis yet */

  errcod = allocbn_scf_irreps(_irreps,"nirrep",1);
  if(errcod != 0) {
    fprintf(_outfile,"scf_initialize: could not allocate _irreps\n");
    return(-1);
    }

  errcod = allocbn_scf_irrep(&_irreps->ir[0],"degeneracy irrep_label",1,"A");
  if(errcod != 0) {
    fprintf(_outfile,
      "scf_initialize: could not allocate _irreps->ir[0]\n");
    return(-1);
    }

  _irreps->ir[0].num_so=nbasis;
  _irreps->ir[0].num_so_tri = ntri;

/* fill in some redundant info */

  if(_scf_info->hsos || _scf_info->singlet || _scf_info->twocon ||
     _scf_info->special) _scf_info->iopen=1;

  _scf_info->nbfao = nbasis;
  _scf_info->nbatri = ntri;

 /* for now, use cartesian functions, so nbfso=nbfao */
  _scf_info->nbfso = nbasis;
  ntri = nbasis*(nbasis+1)/2;
  _scf_info->nbstri = ntri;

 /* nsomax is the size of the largest symmetry block
  * mxcoef is the sum of the sizes of the symmetry blocks squared
  * mxcoef2 is the ioff version of mxcoef
  * ideg is vestigal
  */
  for(i=0; i < _irreps->nirrep ; i++) {
    _scf_info->nsomax = (_irreps->ir[i].num_so > _scf_info->nsomax) ?
                          _irreps->ir[i].num_so : _scf_info->nsomax;
    _scf_info->mxcoef += _irreps->ir[i].num_so*_irreps->ir[i].num_so;
    _scf_info->mxcoef2 += _irreps->ir[i].num_so_tri;
    if(i) _irreps->ir[i].ideg = _irreps->ir[i-1].num_so+_irreps->ir[i-1].ideg;
    }

/* deallocate integral stuff */
  int_done_offsets1(_centers,_centers);
  int_done_1e();

  return(0);
  }
