
/* $Log$
 * Revision 1.1  1993/12/29 12:53:14  etseidl
 * Initial revision
 *
 * Revision 1.1.1.1  1992/03/17  17:11:17  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  17:11:16  seidl
 * Initial revision
 *
 * Revision 1.2  1992/01/21  11:33:45  seidl
 * use libintv2
 *
 * Revision 1.1  1991/12/20  15:45:23  seidl
 * Initial revision
 *
 * Revision 1.1  1991/12/03  16:50:37  etseidl
 * Initial revision
 *
 * Revision 1.1  1991/12/02  19:55:10  seidl
 * Initial revision
 * */

static char rcsid[] = "$Id$";

#define IOFF(a) (a)*((a)+1)/2

#include <stdio.h>
#include <math.h>
#include <tmpl.h>
#include <math/array/math_lib.h>
#include <util/ipv2/ip_libv2.h>
#include <chemistry/qc/intv2/int_libv2.h>

#include "symm.h"
#include "symerr.gbl"

#include "sympack.gbl"
#include "sympack.lcl"

GLOBAL_FUNCTION VOID
sym_dm_t_svec(_mat,_svec)
double_matrix_t *_mat;
sym_d_vector_t *_svec;
{
  int i,j,ij;
  int irrep;
  int nso,nsot;
  int ioff=0;
  double_vector_t *v;

  for(irrep=0; irrep < _svec->nirrep ; irrep++) {
    v = &_svec->ir[irrep];
    if(nsot=v->n) {
      nso=(int)sqrt(1.0+(double)nsot*8)/2;
      for(i=ij=0; i < nso ; i++) {
        for(j=0; j <= i ; j++,ij++) {
          v->d[ij]=_mat->d[ioff+i][ioff+j];
          }
        }
      ioff += nso;
      }
    }
  }

GLOBAL_FUNCTION VOID
sym_svec_t_dm(_mat,_svec)
double_matrix_t *_mat;
sym_d_vector_t *_svec;
{
  int i,j,ij;
  int irrep;
  int nso,nsot;
  int ioff=0;
  double_vector_t *v;

  for(irrep=0; irrep < _svec->nirrep ; irrep++) {
    v = &_svec->ir[irrep];
    if(nsot=v->n) {
      nso=(int)sqrt(1.0+(double)nsot*8)/2;
      for(i=ij=0; i < nso ; i++) {
        for(j=0; j <= i ; j++,ij++) {
          _mat->d[ioff+i][ioff+j]=v->d[ij];
          }
        }
      ioff += nso;
      }
    }
  }

GLOBAL_FUNCTION VOID
sym_dm_t_smat(_mat,_smat)
double_matrix_t *_mat;
sym_d_matrix_t *_smat;
{
  int i,j;
  int irrep;
  int nso;
  int ioff=0;
  int joff=0;
  double_matrix_t *m;

  for(irrep=0; irrep < _smat->nirrep ; irrep++) {
    m = &_smat->ir[irrep];
    if(nso=m->n1) {
      for(i=0; i < nso ; i++) {
        for(j=0; j < m->n2 ; j++) {
          m->d[i][j]=_mat->d[ioff+i][joff+j];
          }
        }
      ioff += nso;
      joff += m->n2;
      }
    }
  }

GLOBAL_FUNCTION VOID
sym_smat_t_dm(_mat,_smat)
double_matrix_t *_mat;
sym_d_matrix_t *_smat;
{
  int i,j;
  int irrep;
  int nso;
  int ioff=0;
  int joff=0;
  double_matrix_t *m;

  for(irrep=0; irrep < _smat->nirrep ; irrep++) {
    m = &_smat->ir[irrep];
    if(nso=m->n1) {
      for(i=0; i < nso ; i++) {
        for(j=0; j < m->n2 ; j++) {
          _mat->d[ioff+i][joff+j]=m->d[i][j];
          }
        }
      ioff += nso;
      joff += m->n2;
      }
    }
  }

GLOBAL_FUNCTION VOID
sym_dv_t_svec(_vec,_svec)
double_vector_t *_vec;
sym_d_vector_t *_svec;
{
  int i,j,ij;
  int irrep;
  int nso,nsot;
  int ioff=0;
  double_vector_t *v;

  for(irrep=0; irrep < _svec->nirrep ; irrep++) {
    v = &_svec->ir[irrep];
    if(nsot=v->n) {
      nso=(int)sqrt(1.0+(double)nsot*8)/2;
      for(i=ij=0; i < nso ; i++) {
        for(j=0; j <= i ; j++,ij++) {
          v->d[ij]=_vec->d[IOFF((ioff+i))+ioff+j];
          }
        }
      ioff += nso;
      }
    }
  }

GLOBAL_FUNCTION VOID
sym_svec_t_dv(_vec,_svec)
double_vector_t *_vec;
sym_d_vector_t *_svec;
{
  int i,j,ij;
  int irrep;
  int nso,nsot;
  int ioff=0;
  double_vector_t *v;

  for(irrep=0; irrep < _svec->nirrep ; irrep++) {
    v = &_svec->ir[irrep];
    if(nsot=v->n) {
      nso=(int)sqrt(1.0+(double)nsot*8)/2;
      for(i=ij=0; i < nso ; i++) {
        for(j=0; j <= i ; j++,ij++) {
          _vec->d[IOFF((ioff+i))+ioff+j]=v->d[ij];
          }
        }
      ioff += nso;
      }
    }
  }
