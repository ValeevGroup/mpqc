
/* $Log$
 * Revision 1.1  1993/12/29 12:53:04  etseidl
 * Initial revision
 *
 * Revision 1.1.1.1  1992/03/17  17:09:16  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  17:09:14  seidl
 * Initial revision
 *
 * Revision 1.5  1992/01/16  19:53:07  seidl
 * use new integral routines
 *
 * Revision 1.4  1992/01/09  11:48:38  seidl
 * add parallel code, master file pointers initialized in scf_open_files
 *
 * Revision 1.3  1991/12/20  16:31:11  seidl
 * fix bug in sq_to_tri
 *
 * Revision 1.2  1991/12/17  21:44:48  seidl
 * change to new util/sgen/sgen format
 *
 * Revision 1.1  1991/12/02  19:58:51  seidl
 * Initial revision
 * */

static char rcsid[] = "$Id$";

#include <stdio.h>
#include <tmpl.h>
#include <math/array/math_lib.h>
#include <util/bio/libbio.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <chemistry/qc/sym/sym_lib.h>
#include "scf.h"
#include "scfallc.h"
#include "scfbrd.h"
#include "scfbwr.h"
#include "scffree.h"
#include "scfinit.h"
#include "scfprnt.h"

#include "scf_oeis.gbl"
#include "scf_oeis.lcl"

GLOBAL_FUNCTION int
scf_oeis(_scf_info,_sym_info,_irreps,_centers,_mfp,_mast,_outfile)
scf_struct_t *_scf_info;
sym_struct_t *_sym_info;
scf_irreps_t *_irreps;
centers_t *_centers;
scf_mf_ptrs_t *_mfp;
int _mast;
FILE *_outfile;
{
  int i,j,k,l;
  int *mptr;
  int errcod;
  int nbasis=_scf_info->nbfao;
  int ntri=_scf_info->nbatri;
  double *junk;

  double_matrix_t uaoso,usoao;
  double_matrix_t scr,scr2;
  sym_d_vector_t sscr,sscr2;

  if(_scf_info->use_symmetry) {
    errcod = allocbn_double_matrix(&uaoso,"n1 n2",nbasis,nbasis);
    if(errcod != 0) {
      fprintf(_outfile,"scf_oeis:\n");
      fprintf(_outfile,"could not allocate memory for uaoso matrix\n");
      return(-1);
      }
    errcod = allocbn_double_matrix(&usoao,"n1 n2",nbasis,nbasis);
    if(errcod != 0) {
      fprintf(_outfile,"scf_oeis:\n");
      fprintf(_outfile,"could not allocate memory for usoao matrix\n");
      return(-1);
      }
    errcod = allocbn_double_matrix(&scr2,"n1 n2",nbasis,nbasis);
    if(errcod != 0) {
      fprintf(_outfile,"scf_oeis:\n");
      fprintf(_outfile,"could not allocate memory for scr2 matrix\n");
      return(-1);
      }
    }

  errcod = allocbn_double_matrix(&scr,"n1 n2",nbasis,nbasis);
  if(errcod != 0) {
    fprintf(_outfile,"scf_oeis:\n");
    fprintf(_outfile,"could not allocate memory for scr matrix\n");
    return(-1);
    }
  errcod = alloc_s_d_vector(&sscr,_irreps);
  if(errcod != 0) {
    fprintf(_outfile,"scf_oeis:\n");
    fprintf(_outfile,"could not allocate memory for sscr vector\n");
    return(-1);
    }
  errcod = alloc_s_d_vector(&sscr2,_irreps);
  if(errcod != 0) {
    fprintf(_outfile,"scf_oeis:\n");
    fprintf(_outfile,"could not allocate memory for sscr2 vector\n");
    return(-1);
    }

/* initialize _centers struct */

#if defined(I860)
  junk = int_initialize_1e(0,0,_centers,_centers);
#endif
  int_initialize_offsets1(_centers,_centers);

/* get pointer to beginning of one-electron part of master file */

  mptr = &_mfp->cur_ptr;

/* get ao to so transformation matrices */

  if(_scf_info->use_symmetry) {
    errcod = sym_aotoso(_centers,_sym_info,_irreps,&uaoso,&usoao,_outfile);
    if(errcod != 0) {
      fprintf(_outfile,"scf_oeis:\n");
      fprintf(_outfile,"trouble forming so to ao matrices\n");
      return(-1);
      }

    *mptr = _mfp->uaoso;
    errcod = bwrite_double_matrix(_mast,&uaoso,mptr);
    if(errcod < 0) {
      fprintf(_outfile,"scf_oeis:\n");
      fprintf(_outfile,"trouble writing uaoso to master file\n");
      return(errcod);
      }
    *mptr = _mfp->usoao;
    errcod = bwrite_double_matrix(_mast,&usoao,mptr);
    if(errcod < 0) {
      fprintf(_outfile,"scf_oeis:\n");
      fprintf(_outfile,"trouble writing usoao to master file\n");
      return(errcod);
      }
    }

/* calculate one-electron integrals and write to master file */

  *mptr = _mfp->overlap;
  int_overlap(_centers,_centers,&scr);

  if(_scf_info->use_symmetry) {
    math_dmxdm_dm(&uaoso,1,&scr,0,&scr2,0,nbasis,nbasis,nbasis,0);
    math_dmxdm_dm(&scr2,0,&uaoso,0,&scr,0,nbasis,nbasis,nbasis,0);
    sym_dm_t_svec(&scr,&sscr);
    }
  else {
    sq_to_tri(&scr,&sscr.ir[0]);
    }
  errcod = bwrite_sym_d_vector(_mast,&sscr,mptr);
  if(errcod < 0) {
    fprintf(_outfile,"scf_oeis:\n");
    fprintf(_outfile,"trouble writing overlap ints to master file\n");
    return(errcod);
    }

  *mptr = _mfp->kinetic;
  int_kinetic(_centers,_centers,&scr);

  if(_scf_info->use_symmetry) {
    math_dmxdm_dm(&uaoso,1,&scr,0,&scr2,0,nbasis,nbasis,nbasis,0);
    math_dmxdm_dm(&scr2,0,&uaoso,0,&scr,0,nbasis,nbasis,nbasis,0);
    sym_dm_t_svec(&scr,&sscr);
    }
  else {
    sq_to_tri(&scr,&sscr.ir[0]);
    }
  errcod = bwrite_sym_d_vector(_mast,&sscr,mptr);
  if(errcod < 0) {
    fprintf(_outfile,"scf_oeis:\n");
    fprintf(_outfile,"trouble writing kinetic energy ints to master file\n");
    return(errcod);
    }

  *mptr = _mfp->nuclear;
  int_nuclear(_centers,_centers,&scr);

  if(_scf_info->use_symmetry) {
    math_dmxdm_dm(&uaoso,1,&scr,0,&scr2,0,nbasis,nbasis,nbasis,0);
    math_dmxdm_dm(&scr2,0,&uaoso,0,&scr,0,nbasis,nbasis,nbasis,0);
    sym_dm_t_svec(&scr,&sscr);
    }
  else {
    sq_to_tri(&scr,&sscr.ir[0]);
    }
  errcod = bwrite_sym_d_vector(_mast,&sscr,mptr);
  if(errcod < 0) {
    fprintf(_outfile,"scf_oeis:\n");
    fprintf(_outfile,"trouble writing potential energy ints to master file\n");
    return(errcod);
    }

  *mptr = _mfp->kinetic;

  errcod = bread_sym_d_vector(_mast,&sscr2,mptr);
  if(errcod < 0) {
    fprintf(_outfile,"scf_oeis:\n");
    fprintf(_outfile,"trouble reading kinetic energy ints from master file\n");
    return(errcod);
    }

  for(i=0; i < _irreps->nirrep ; i++)
#if 0
    add_double_vector(&sscr.ir[i],&sscr2.ir[i],&sscr.ir[i]);
#else
    for(j=0; j < sscr.ir[i].n ; j++)
      sscr.ir[i].d[j] += sscr2.ir[i].d[j];
#endif

  *mptr = _mfp->hcore;
  errcod = bwrite_sym_d_vector(_mast,&sscr,mptr);
  if(errcod < 0) {
    fprintf(_outfile,"scf_oeis:\n");
    fprintf(_outfile,"trouble writing core hamiltonian to master file\n");
    return(errcod);
    }

  i=sizeof(int);
  errcod = bwrite_scf_mf_ptrs(_mast,_mfp,&i);
  if(errcod < 0) {
    fprintf(_outfile,"scf_oeis:\n");
    fprintf(_outfile,"trouble writing matrix pointers to master file\n");
    return(errcod);
    }

  int_done_offsets1(_centers,_centers);
  int_done_1e();

  free_double_matrix(&scr);
  free_sym_d_vector(&sscr);
  free_sym_d_vector(&sscr2);

  if(_scf_info->use_symmetry) {
    free_double_matrix(&usoao);
    free_double_matrix(&uaoso);
    free_double_matrix(&scr2);
    }

  return(0);
  }

LOCAL_FUNCTION VOID
sq_to_tri(mat,vec)
double_matrix_t *mat;
double_vector_t *vec;
{
  int i,j,ij;

  if(mat->n1 != mat->n2) {
    fprintf(stderr,"scf_oeis: sq_to_tri:\n");
    fprintf(stderr,"n1 != n2. this no good\n");
    exit(1);
    }

  for(i=ij=0; i < mat->n1 ; i++)
    for(j=0; j <= i ; j++,ij++)
      vec->d[ij] = mat->d[i][j];
  }
