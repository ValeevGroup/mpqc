
/* $Log$
 * Revision 1.1  1993/12/29 12:53:04  etseidl
 * Initial revision
 *
 * Revision 1.1.1.1  1992/03/17  17:09:03  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  17:09:02  seidl
 * Initial revision
 *
 * Revision 1.5  1992/01/16  19:53:07  seidl
 * use new integral routines
 *
 * Revision 1.4  1992/01/13  19:15:41  seidl
 * ao density formed in scf_make_gmat now
 *
 * Revision 1.3  1992/01/09  11:54:02  seidl
 * no longer include util/ipv2/ip_libv2.h
 *
 * Revision 1.2  1992/01/02  16:20:37  seidl
 * multiply diagonal elements of pmat by 2 before transforming to ao basis
 *
 * Revision 1.1  1991/12/20  16:23:13  seidl
 * Initial revision
 *
 * Revision 1.1  1991/06/22  09:01:31  seidl
 * Initial revision
 * */

static char rcsid[] = "$Id$";

#include <stdio.h>
#include <tmpl.h>
#include <math/array/math_lib.h>
#include <math.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <chemistry/qc/sym/sym_lib.h>
#include "scf.h"

/* calculates G matrix, where G(ij)=sum(k,l) P(kl)*T(ijkl)
 * P is the density difference matrix, so G is really DELTA G
 */

#include "scf_mkg.gbl"
#include "scf_mkg.lcl"

#define ioff(i) i*(i+1)/2
#define MIN0(a,b) ((a)<(b)) ? (a) : (b)
#define MAX0(a,b) ((a)>(b)) ? (a) : (b)

int wherec=0;
int whereo=0;
int *int_nums_c;
int *int_nums_o;

GLOBAL_FUNCTION int 
scf_make_g(_centers,_irreps,_scf_info,_sym_info,_mfp,_mast,
            _gmat,_gmato,_dpmat,_dpmato,scr_c,scr_o,spkc,spko,_outfile)
centers_t *_centers;
scf_irreps_t *_irreps;
scf_struct_t *_scf_info;
sym_struct_t *_sym_info;
scf_mf_ptrs_t *_mfp;
int _mast;
double_vector_t *_gmat;
double_vector_t *_gmato;
double_vector_t *_dpmat;
double_vector_t *_dpmato;
int scr_c;
int scr_o;
scf_pk_buf_c_t *spkc;
scf_pk_buf_o_t *spko;
FILE *_outfile;
{
  register int i,j,k,joff,nn;
  register int ij,kl;
  int num;
  int iopen=_scf_info->iopen;
  int errcod;
  double dotest,tmpval;
  scf_pkint_o_t *o_temp;
  scf_pkint_c_t *c_temp;
  scf_irrep_t *s;

  if(_scf_info->lastc + _scf_info->lasto == 0) return 0;

  if(!wherec) {
    int_nums_c = (int *) init_array(sizeof(int)*(_scf_info->num_bufs_c+1));
    for(i=1; i < _scf_info->num_bufs_c ; i++) int_nums_c[i]=_scf_info->maxbufc;
    int_nums_c[_scf_info->num_bufs_c]=_scf_info->lastc;
    wherec=_scf_info->num_bufs_c;
    if(iopen) {
      int_nums_o = (int *) init_array(sizeof(int)*(_scf_info->num_bufs_o+1));
      for(i=1; i < _scf_info->num_bufs_o ; i++) 
        int_nums_o[i]=_scf_info->maxbufo;
      int_nums_o[_scf_info->num_bufs_o]=_scf_info->lasto;
      whereo=_scf_info->num_bufs_o;
      }
    }

  if(iopen) {
    num=int_nums_o[whereo];
    for (j=0; j < _scf_info->num_bufs_o ; j++) {
      o_temp = &spko->pk[0];
  
      for (i=num; i ; i--,o_temp++) {
        ij = (*o_temp).ij;
        kl = (*o_temp).kl;
        tmpval = (*o_temp).pval;
        dotest = (*o_temp).kval;

        _gmat->d[ij] += _dpmat->d[kl]*tmpval;
        _gmat->d[kl] += _dpmat->d[ij]*tmpval;
        _gmato->d[ij] += _dpmato->d[kl]*dotest;
        _gmato->d[kl] += _dpmato->d[ij]*dotest;
        }
   
      if (_scf_info->readflgo && j < _scf_info->num_bufs_o-1) {
        if(whereo==_scf_info->num_bufs_o) {
          bio_rewind(scr_o);
          whereo=0;
          }
        whereo++;
        num=int_nums_o[whereo];
        bio_read(scr_o,spko->pk,sizeof(scf_pkint_o_t)*num);
        }
      }
    }

  num=int_nums_c[wherec];
  for (j=0; j < _scf_info->num_bufs_c ; j++) {
    c_temp = &spkc->pk[0];

    for (i=num; i ; i--,c_temp++) {
      ij = (*c_temp).ij;
      kl = (*c_temp).kl;
      tmpval = (*c_temp).pval;

      _gmat->d[ij] += _dpmat->d[kl]*tmpval;
      _gmat->d[kl] += _dpmat->d[ij]*tmpval;
      }
   
    if (_scf_info->readflgc && j < _scf_info->num_bufs_c-1) {
      if(wherec==_scf_info->num_bufs_c) {
        bio_rewind(scr_c);
        wherec=0;
        }
      wherec++;
      num=int_nums_c[wherec];
      bio_read(scr_c,spkc->pk,sizeof(scf_pkint_c_t)*num);
      }
    }

  return 0;
  }

LOCAL_FUNCTION VOID_PTR
init_array(n)
int n;
{
  void *arr;

  arr = (void *) malloc(n);
  if(arr==NULL) {
    fprintf(stderr,"scf_make_pk: init_array:\n");
    fprintf(stderr,"malloc trouble\n");
    exit(99);
    }
  bzero(arr,n);
  return arr;
  }
