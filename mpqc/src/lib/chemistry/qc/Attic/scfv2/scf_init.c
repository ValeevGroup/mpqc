
/* $Log$
 * Revision 1.1  1993/12/29 12:53:04  etseidl
 * Initial revision
 *
 * Revision 1.1.1.1  1992/03/17  17:08:46  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  17:08:45  seidl
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
#include <chemistry/qc/sym/sym_lib.h>
#include "scf.h"

#include "scf_init.gbl"
#include "scf_init.lcl"

#define CLOSED 0
#define OPENSH 1
#define TWOCON 2

GLOBAL_FUNCTION int
scf_initialize(_outfile,_centers,_scf_info,_sym_info,_irreps,pg)
FILE *_outfile;
centers_t *_centers;
scf_struct_t *_scf_info;
sym_struct_t *_sym_info;
scf_irreps_t *_irreps;
char *pg;
{
  int i,j,k,l;
  int errcod;
  int nbasis,ntri;
  double nso,nf,rj;
  double_matrix_t scr;
  double_vector_t trace;
  center_t *c;

/* initialize the centers structure */

/* first read in information about atoms, basis set, etc, and place
 * in the centers struct */

  errcod = sym_init_centers(_centers,_sym_info,pg,_outfile);
  if(errcod != 0) {
    fprintf(_outfile,"scf_initialize: trouble symmetrizing centers\n");
    return(-1);
    }

  nbasis = _centers->nfunc;
  ntri = nbasis*(nbasis+1)/2;

 /* now form the _irreps struct */

  if(_scf_info->use_symmetry) {
    errcod = allocbn_scf_irreps(_irreps,"nirrep",_sym_info->nirrep);
    if(errcod != 0) {
      fprintf(_outfile,"scf_initialize: could not allocate _irreps\n");
      return(-1);
      }

    for(i=0; i < _sym_info->nirrep ; i++) {
      errcod = allocbn_scf_irrep(&_irreps->ir[i],
       "degeneracy irrep_label",
         _sym_info->ct.gamma[i].degen,_sym_info->ct.gamma[i].label);
      if(errcod != 0) {
        fprintf(_outfile,
          "scf_initialize: could not allocate _irreps->ir[%d]\n",i);
        return(-1);
        }
      }

  /* now figure out how many so's of each symmetry there will be */

    errcod = allocbn_double_matrix(&scr,"n1 n2",nbasis,nbasis);
    if(errcod != 0) {
      fprintf(_outfile, "scf_initialize: could not allocate scr");
      return(-1);
      }
    errcod = allocbn_double_vector(&trace,"n",_sym_info->g);
    if(errcod != 0) {
      fprintf(_outfile, "scf_initialize: could not allocate trace");
      return(-1);
      }

    for(i=0; i < _sym_info->g ; i++) {
      sym_create_r(_centers,_sym_info,&scr,i,_outfile);
  
      trace.d[i]=0.0;
      for(j=0; j < nbasis ;j++) trace.d[i] += scr.d[j][j];
      }

    for(i=0; i < _sym_info->nirrep ; i++) {
      nso=0.0;
      for(j=0; j < _sym_info->g ; j++) 
        nso += trace.d[j]*_sym_info->ct.gamma[i].rep[j];
      nso /= (double) _sym_info->g;

      if(_sym_info->pg == _PG_CN || _sym_info->pg == _PG_SN || 
         _sym_info->pg == _PG_CNH) nso /= (double) _irreps->ir[i].degeneracy;
#if 1
      nso *= _irreps->ir[i].degeneracy;
#endif

      nf = modf(nso,&rj);
      j = (int) rj;
      if(nf >= 0.5) j++;
      _irreps->ir[i].num_so=j;
      _irreps->ir[i].num_so_tri = j*(j+1)/2;
      }

    free_double_matrix(&scr);
    free_double_vector(&trace);
    }

 /* if not using symmetry, then set up _irreps in C1 symmetry */

  else {
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
    }

/* fill in some redundant info */

  if(_scf_info->hsos || _scf_info->singlet || _scf_info->twocon ||
     _scf_info->special) _scf_info->iopen=1;

  _scf_info->nbfao = nbasis;
  _scf_info->nbatri = ntri;

  nbasis=0;
  for(i=0; i < _centers->n ; i++) {
    c = &_centers->center[i];
    for(j=0; j < c->basis.n ; j++) {
      for(k=0; k < c->basis.shell[j].ncon ; k++) {
        int am = c->basis.shell[j].type[k].am;
        nbasis += 2*am+1;
        }
      }
    }

  _scf_info->nbfso = nbasis;
  ntri = nbasis*(nbasis+1)/2;
  _scf_info->nbstri = ntri;

  for(i=0; i < _irreps->nirrep ; i++) {
    _scf_info->nsomax = (_irreps->ir[i].num_so > _scf_info->nsomax) ?
                          _irreps->ir[i].num_so : _scf_info->nsomax;
    _scf_info->mxcoef += _irreps->ir[i].num_so*_irreps->ir[i].num_so;
    _scf_info->mxcoef2 += _irreps->ir[i].num_so_tri;
    if(i) _irreps->ir[i].ideg = _irreps->ir[i-1].num_so+_irreps->ir[i-1].ideg;
    }

/* calculate nuclear repulsion energy */

  _scf_info->nuc_rep = (double) int_nuclear_repulsion(_centers,_centers);

  int_done_offsets1(_centers,_centers);

  return(0);
  }
