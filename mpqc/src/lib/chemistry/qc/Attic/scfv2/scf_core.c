
/* $Log$
 * Revision 1.1  1993/12/29 12:53:04  etseidl
 * Initial revision
 *
 * Revision 1.1.1.1  1992/03/17  17:08:25  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  17:08:24  seidl
 * Initial revision
 *
 * Revision 1.6  1992/01/16  19:53:07  seidl
 * use new integral routines
 *
 * Revision 1.5  1992/01/13  19:10:44  seidl
 * remove debug printing, add some comments
 *
 * Revision 1.4  1992/01/09  11:35:59  seidl
 * no longer include util/ipv2/ip_libv2.h
 *
 * Revision 1.3  1991/12/20  16:25:47  seidl
 * add debug print statements, and fix bug in forming scf_vec
 *
 * Revision 1.2  1991/12/17  21:42:21  seidl
 * modified for new util/sgen/sgen
 *
 * Revision 1.1  1991/12/17  14:09:14  seidl
 * Initial revision
 *
 * Revision 1.1  1991/12/02  19:58:51  seidl
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

#include "scf_core.gbl"
#include "scf_core.lcl"


GLOBAL_FUNCTION int
scf_core_guess(_scf_vec,_irreps,_scf_info,_mfp,_mast,_outfile)
sym_d_matrix_t *_scf_vec;
scf_irreps_t *_irreps;
scf_struct_t *_scf_info;
scf_mf_ptrs_t *_mfp;
int _mast;
FILE *_outfile;
{
  int i;
  int ii,jj,kk;
  int *mptr=&_mfp->cur_ptr;
  int errcod;
  int nn,num_irrep,nsomax;
  double tol=1.0e-14;

  double_vector_t evals;
  double_matrix_t ctrans,sahalf;
  sym_d_vector_t overlap;
  sym_d_vector_t hcore;
  scf_irrep_t *ir;

  num_irrep=_irreps->nirrep;
  nsomax = _scf_info->nsomax;
  
  errcod = allocbn_double_vector(&evals,"n",nsomax);
  if(errcod != 0) {
    fprintf(_outfile,"scf_core_guess:\n");
    fprintf(_outfile,"could not allocate memory for evals vector\n");
    return(-1);
    }
  errcod = allocbn_double_matrix(&ctrans,"n1 n2",nsomax,nsomax);
  if(errcod != 0) {
    fprintf(_outfile,"scf_core_guess:\n");
    fprintf(_outfile,"could not allocate memory for ctrans matrix\n");
    return(-1);
    }
  errcod = allocbn_double_matrix(&sahalf,"n1 n2",nsomax,nsomax);
  if(errcod != 0) {
    fprintf(_outfile,"scf_core_guess:\n");
    fprintf(_outfile,"could not allocate memory for sahalf matrix\n");
    return(-1);
    }

/* read in overlap integrals and core hamiltonian from master file
 * bread allocates memory for overlap and hcore, so just init them */

  *mptr = _mfp->overlap;
  init_sym_d_vector(&overlap);
  errcod = bread_sym_d_vector(_mast,&overlap,mptr);
  
  *mptr = _mfp->hcore;
  init_sym_d_vector(&hcore);
  errcod = bread_sym_d_vector(_mast,&hcore,mptr);
  
  for(i=0; i < num_irrep ; i++) {
    ir = &_irreps->ir[i];
    if(nn=ir->num_so) {

  /* diagonalize overlap matrix */
 
      math_diag_dv(&overlap.ir[i],&evals,&ctrans,1,tol,1);

  /* form 'sahalf' matrix sahalf = u*ei^-0.5*u~  */

       for(ii=0; ii < nn ; ii++) {
         for(jj=0; jj < nn ; jj++) {
           sahalf.d[ii][jj]=0.0;
           for(kk=0; kk < nn ; kk++) sahalf.d[ii][jj] +=
             ctrans.d[ii][kk]*ctrans.d[jj][kk]/sqrt(evals.d[kk]);
           }
         }

  /* sahalf is formed in the so basis, note that to obtain S-1/2 in
   * the ao basis, the transformation Ua*S-1/2*Ua~ must be performed
   * where Ua is the aotoso transformation matrix (not sotoao!)
   */

  /* now diagonalize core hamiltonian, and multiply vector by sahalf
   * this way the scf vector will also transform to the S-1/2 basis
   *  c' = S-1/2*c
   *  M(mo) = ~c'*M*c' = ~c*S-1/2*M*S-1/2*c = ~c*M(s-1/2)*c
   */

      math_diag_dv(&hcore.ir[i],&evals,&ctrans,1,tol,1);
      math_dmxdm_dm(&sahalf,0,&ctrans,0,&_scf_vec->ir[i],0,nn,nn,nn,0);
      }
    }

  *mptr = _mfp->scf_vector;
  errcod = bwrite_sym_d_matrix(_mast,_scf_vec,mptr);
  if(errcod < 0) {
    fprintf(_outfile,"scf_core_guess:\n");
    fprintf(_outfile,"could not write scf vector to master file\n");
    return(-1);
    }

  if(!_scf_info->restart) {
    i=sizeof(int);
    errcod = bwrite_scf_mf_ptrs(_mast,_mfp,&i);
    }

  free_double_vector(&evals);
  free_double_matrix(&ctrans);
  free_double_matrix(&sahalf);
  free_sym_d_vector(&overlap);
  free_sym_d_vector(&hcore);

  return 0;
  }

