
/* $Log$
 * Revision 1.1  1993/12/29 12:53:04  etseidl
 * Initial revision
 *
 * Revision 1.1.1.1  1992/03/17  17:08:38  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  17:08:37  seidl
 * Initial revision
 *
 * Revision 1.3  1992/01/16  19:53:07  seidl
 * use new integral routines
 *
 * Revision 1.2  1992/01/09  11:44:36  seidl
 * no longer include util/ipv2/ip_libv2.h
 *
 * Revision 1.1  1992/01/02  16:19:07  seidl
 * Initial revision
 * */

static char rcsid[]="$Id$";

#include <stdio.h>
#include <tmpl.h>
#include <math/array/math_lib.h>
#include <math.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <chemistry/qc/sym/sym_lib.h>
#include "scf.h"

#include "scf_fock.gbl"
#include "scf_fock.lcl"


GLOBAL_FUNCTION int
scf_make_fock(_irreps,_scf_info,_fock,_focko,_scf_vec,_occ_num,_tfock,_tfocko,
              iter,_outfile)
scf_irreps_t *_irreps;
scf_struct_t *_scf_info;
sym_d_vector_t *_fock;
sym_d_vector_t *_focko;
sym_d_matrix_t *_scf_vec;
sym_d_vector_t *_occ_num;
double_vector_t *_tfock;
double_vector_t *_tfocko;
int iter;
FILE *_outfile;
{
  int i,j,ij,m;
  int nn;
  int errcod;
  double occi, occj, occ0;
  double_matrix_t scr1;
  scf_irrep_t *s;

  errcod = 
    allocbn_double_matrix(&scr1,"n1 n2",_scf_info->nsomax,_scf_info->nsomax);
  if(errcod != 0) {
    fprintf(_outfile,"scf_make_fock: ");
    fprintf(_outfile,"could not allocate memory for scr1\n");
    return(-1);
    }

  for (m=0; m < _irreps->nirrep ; m++) {
    s = &_irreps->ir[m];
    if (nn=s->num_so) {

   /* transform fock to mo basis */
      math_dmxdv_dm(&_scf_vec->ir[m],1,&_fock->ir[m],0,&scr1,0,nn,nn,nn,0);
      math_dmxdm_dv(&scr1,0,&_scf_vec->ir[m],0,_tfock,0,nn,nn,nn,0);

   /* transform fock_open to mo basis */
      math_dmxdv_dm(&_scf_vec->ir[m],1,&_focko->ir[m],0,&scr1,0,nn,nn,nn,0);
      math_dmxdm_dv(&scr1,0,&_scf_vec->ir[m],0,_tfocko,0,nn,nn,nn,0);

   /* form effective fock matrix in mo basis */

      zero_double_matrix(&scr1);

      occ0 = _occ_num->ir[m].d[0];
      for (i=ij=0; i < nn; i++ ) {
        for (j=0; j <= i; j++,ij++) {
          occi = _occ_num->ir[m].d[i];
          occj = _occ_num->ir[m].d[j];

     /* default: Guest & Saunders general form */
          if(iter < _scf_info->maxiter-1 && 
             !_scf_info->converged && !_scf_info->fock_typ) {
            if(occi == occj) 
              _fock->ir[m].d[ij] = _tfock->d[ij];
            else if(occi)
              _fock->ir[m].d[ij] = 2.0*_tfock->d[ij]-_tfocko->d[ij];
            else {
              if(occj==2.0)
                _fock->ir[m].d[ij] = _tfock->d[ij];
              else
                _fock->ir[m].d[ij] = _tfocko->d[ij];
              }
            }

     /* Guest & Saunders' form for high spin */
          else if(iter < _scf_info->maxiter-1 && !_scf_info->converged 
                                            && _scf_info->fock_typ == 1) {
            if (occi == occj || occi)
              _fock->ir[m].d[ij] = 2.0*_tfock->d[ij]-_tfocko->d[ij];
            else if (occj == 2.0) _fock->ir[m].d[ij] = _tfock->d[ij];
            else _fock->ir[m].d[ij] = _tfocko->d[ij];
            }

     /* test form (fo fo fo) */
          else if(iter < _scf_info->maxiter-1 && !_scf_info->converged &&
                                                _scf_info->fock_typ == 2) {
            if (occi == occj) _fock->ir[m].d[ij] = _tfocko->d[ij];
            else if(occi)
              _fock->ir[m].d[ij] = 2.0*_tfock->d[ij]-_tfocko->d[ij];
            else if(occj == 2.0)
              _fock->ir[m].d[ij] = _tfock->d[ij];
            else _fock->ir[m].d[ij] = _tfocko->d[ij];
            }

     /* form for converged wavefunction */
          else {
            if (occi == 2.0) _fock->ir[m].d[ij]=_tfock->d[ij];
            else if(occj != 2.0)
              _fock->ir[m].d[ij]=_tfocko->d[ij];
            else {
              if(occi)
                _fock->ir[m].d[ij] = 2.0*_tfock->d[ij]-_tfocko->d[ij];
              else _fock->ir[m].d[ij]=_tfock->d[ij];
              }
            }

          if(j==i) {
            if (occi == occ0 && occi)
               _fock->ir[m].d[ij] -= _scf_info->lvl_shift;
            else 
               if (occi) _fock->ir[m].d[ij] -= 0.5*_scf_info->lvl_shift;
            }
          }
        }
      }
    }

  free_double_matrix(&scr1);
  return(0);
  }
