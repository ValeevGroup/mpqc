
/* This function will perform the diis extrapolation a la Pulay.
 * For more information read Hamilton and Pulay, J.Chem.Phys 84 (1986), 5728.
 */

/* $Log$
 * Revision 1.1  1993/12/29 12:53:04  etseidl
 * Initial revision
 *
 * Revision 1.1.1.1  1992/03/17  17:08:28  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  17:08:27  seidl
 * Initial revision
 *
 * Revision 1.5  1992/01/16  19:53:07  seidl
 * use new integral routines
 *
 * Revision 1.4  1992/01/13  19:11:10  seidl
 * add some comments
 *
 * Revision 1.3  1992/01/09  11:37:07  seidl
 * no longer include util/ipv2/ip_libv2.h
 *
 * Revision 1.2  1992/01/02  16:24:54  seidl
 * fix Log
 * */

static char rcsid[]="$Id$";

#include <stdio.h>
#include <tmpl.h>
#include <math.h>
#include <math/array/math_lib.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <chemistry/qc/sym/sym_lib.h>
#include "scf.h"

#include "scf_diis.gbl"
#include "scf_diis.lcl"

static double_vector_t btemp;
static double_matrix_t bold,bmat;

static struct diis_mats {
  sym_d_vector_t fock_c;
  sym_d_vector_t fock_o;
  sym_d_vector_t error;
  } *diism,dtemp;

GLOBAL_FUNCTION int
scf_diis(_irreps,_scf_info,_fock,_focko,_scf_vec,_occ_num,iter,
 _scr1,_scr2,_scrv1,_scrv2,_outfile)
scf_irreps_t *_irreps;
scf_struct_t *_scf_info;
sym_d_vector_t *_fock;
sym_d_vector_t *_focko;
sym_d_matrix_t *_scf_vec;
sym_d_vector_t *_occ_num;
int iter;
double_matrix_t *_scr1;
double_matrix_t *_scr2;
double_vector_t *_scrv1;
double_vector_t *_scrv2;
FILE *_outfile;
{
  int i,j,k,ij;
  int m,nn;
  int errcod;
  double occi, occj;
  int try = 0;
  int last = iter;
  int col = iter+2;
  int iopen = _scf_info->iopen;
  double etemp, dotp, norm, determ, etempo;
  double_matrix_t *scfv;
  double_vector_t *fc,*fo;
  double scale;
  scf_irrep_t *s;
  struct diis_mats *d;

  iter++;

/* Allocate memory for the fock and error matrices for the last "ndiis"
 * iterations.  This is going to lead to memory problems eventually,
 * so it will probably become necessary to transfer all this to disk
 */

  if(diism == NULL) {
    errcod = allocbn_double_matrix(&bmat,"n1 n2",_scf_info->ndiis+1,
                                                      _scf_info->ndiis+1);
    errcod = allocbn_double_matrix(&bold,"n1 n2",_scf_info->ndiis,
                                                      _scf_info->ndiis);
    errcod = allocbn_double_vector(&btemp,"n",_scf_info->ndiis+1);

    diism = 
      (struct diis_mats *) malloc(sizeof(struct diis_mats)*_scf_info->ndiis);

    for(i=0; i < _scf_info->ndiis ; i++) {
      errcod = alloc_s_d_vector(&diism[i].error,_irreps);
      if(errcod!=0) {
        fprintf(_outfile,"scf_diis: ");
        fprintf(_outfile,"could not allocate memory for error matrix\n");
        return(-1);
        }
      }
    }

  scale = 1.0 + _scf_info->diisdamp;

  if (iter > _scf_info->ndiis) {
    last = _scf_info->ndiis-1;
    col = _scf_info->ndiis+1;
    dtemp = diism[0];
    for (i=0; i < last ; i++) {
      diism[i] = diism[i+1];
      }
    diism[last] = dtemp;
    }
      
 /* save ao fock matrices in fock_save */
   
  d = &diism[last];
  if(iter > _scf_info->ndiis) {
    free_sym_d_vector(&d->fock_c);
    init_sym_d_vector(&d->fock_c);
    if(iopen) {
      free_sym_d_vector(&d->fock_o);
      init_sym_d_vector(&d->fock_o);
      }
    }
  assign_sym_d_vector(&d->fock_c,_fock);
  if(iopen) assign_sym_d_vector(&d->fock_o,_focko);

  _scf_info->diis_er=0.0;
  for (m=0; m < _irreps->nirrep ; m++) {
    s = &_irreps->ir[m];
    scfv = &_scf_vec->ir[m];
    fc = &_fock->ir[m];
    fo = &_focko->ir[m];

    if(nn=s->num_so) {

 /* Form error matrix in mo basis, place in _scr1.
  * This error matrix is just the off-diagonal elements of the fock matrix 
  * rather than FDS-SDF recommended in the 1982 DIIS paper by Pulay.
  */

      math_dmxdv_dm(scfv,1,fc,0,_scr1,0,nn,nn,nn,0);
      math_dmxdm_dv(_scr1,0,scfv,0,_scrv1,0,nn,nn,nn,0);
      if(_scf_info->iopen) {
        math_dmxdv_dm(scfv,1,fo,0,_scr1,0,nn,nn,nn,0);
        math_dmxdm_dv(_scr1,0,scfv,0,_scrv2,0,nn,nn,nn,0);
        }

      for (i=ij=0; i < nn; i++) {
        occi = _occ_num->ir[m].d[i];
        for (j=0; j <= i ; j++,ij++) {
          occj = _occ_num->ir[m].d[j];
          if (!iopen) {
            if (occi == 0.0 && occj != 0.0 ) {
              _scr1->d[i][j]= _scrv1->d[ij];
              _scr1->d[j][i]= _scrv1->d[ij];
              }
            else {
              _scr1->d[i][j]=_scr1->d[j][i]=0.0;
              }
            }
          else if(!_scf_info->twocon) {
            if ((occi == 1.0 || occi == 0.5) && occj == 2.0 )
              etemp = _scrv1->d[ij]-0.5*_scrv2->d[ij];
            else if(occi == 0.0) {
              if (occj == 2.0) etemp = _scrv1->d[ij];
              else if (occj != 0.0) etemp = 0.5*_scrv2->d[ij];
              else etemp = 0.0;
              }
            else etemp = 0.0;
            _scr1->d[i][j] = _scr1->d[j][i] = etemp;
            }
          else {
            if (occi != 2.0 && occi != 0.0 && occj == 2.0 )

  /* this should be changed if twocon is ever added
   * 2.0 should be c1[m] and 1.0 should be c2[m] */

              etemp = 2.0*_scrv1->d[ij]-1.0*_scrv2->d[ij];
            else if(occi == 0.0) {
              if (occj == 2.0) etemp = _scrv1->d[ij];
              else if (occj != 0.0) etemp = _scrv2->d[ij];
              else etemp = 0.0;
              }
            else etemp = 0.0;
            _scr1->d[i][j] = _scr1->d[j][i] = etemp;
            }
          }
        }

    /* transform error matrix into ao basis */
      
      math_dmxdm_dm(scfv,0,_scr1,0,_scr2,0,nn,nn,nn,0);
      math_dmxdm_dv(_scr2,0,scfv,1,&d->error.ir[m],0,nn,nn,nn,0);

      for(i=ij=0; i < nn ; i++) {
        for(j=0; j <= i ; j++,ij++) {
          etemp=fabs(d->error.ir[m].d[ij]);
          _scf_info->diis_er = (_scf_info->diis_er > etemp) ? 
                                _scf_info->diis_er : etemp;
          }
        }
      }
    }
               
  /* then set up B matrix, where B(i,j) = <ei|ej> */

  if (iter > _scf_info->ndiis) {
    for (i=0; i < last ; i++) {
      for (j=0; j <= i ; j++) {
        bold.d[i][j]=bold.d[j][i]=bold.d[i+1][j+1];
        }
      }
    }
  for (i=0; i <= last ; i++)
    bold.d[i][last]=bold.d[last][i] = sdot(&diism[i].error,&diism[last].error);

  bmat.d[0][0] = 0.0;
  btemp.d[0] = -1.0;
  norm = 1.0/bold.d[0][0];
  for (i=1; i <= last+1 ; i++) {
    bmat.d[i][0]=bmat.d[0][i] = -1.0;
    btemp.d[i] = 0.0;
    for (j=1; j <= i ; j++) {
      bmat.d[i][j]=bmat.d[j][i] = bold.d[i-1][j-1]*norm;
      if(i==j) bmat.d[i][j] *= scale;
      }
    }

 /* finally, solve the set of linear equations, obtain the coefficients,
  * and form the new fock matrix F= sum(i=1,n) ci*Fi
  */

  if (iter-1) {
    determ = math_lin(&bmat,&btemp,col,1);

  /* test for poorly conditioned equations */
    while (fabs(determ) < 1.0e-19 && try < last) {

      try++;
      col--;

      bmat.d[0][0] = 0.0;
      btemp.d[0] = -1.0;
      norm=1.0/bold.d[try][try];
      for (i=1; i <= _scf_info->ndiis-try ; i++) {
        bmat.d[i][0]=bmat.d[0][i] = -1.0;
        for (j=1; j <= i ; j++) {
          bmat.d[i][j]=bmat.d[j][i]=bold.d[i+try-1][j+try-1]*norm;
          if(i==j) bmat.d[i][j] *= scale;
          }
        btemp.d[i] = 0.0;
        }

      determ = math_lin(&bmat,&btemp,col,1);
      }

    if(fabs(determ) < 10.0e-20) {
      printf(" try %d no good\n",try);
      return(-1);
      }

    if(iter >= _scf_info->it_diis) {
      for(m=0; m < _irreps->nirrep ; m++) {
        s = &_irreps->ir[m];
        if(nn=s->num_so) {
          for (i=ij=0; i < nn ; i++) {
            for (j=0; j <= i ; j++,ij++) {
              int kk=1;
              etemp=0.0;
              etempo=0.0;
              for (k=try; k < last+1 ; k++) {
                if(iopen) 
                  etempo += btemp.d[kk]*diism[k].fock_o.ir[m].d[ij];
                etemp += btemp.d[kk]*diism[k].fock_c.ir[m].d[ij];
                kk++;
                }
              if(iopen) _focko->ir[m].d[ij] = etempo;
              _fock->ir[m].d[ij] = etemp;
              }
            }
          }
        }
      }
    }
  return(0);
  }

LOCAL_FUNCTION double
sdot(a,b)
sym_d_vector_t *a;
sym_d_vector_t *b;
{
  register int i,j,n,m;
  double *ta, *tb, tval;

  tval = 0.0;
  for(m=0; m < a->nirrep; m++) {
    if(n=a->ir[m].n) {
      ta = a->ir[m].d;
      tb = b->ir[m].d;
      for (i=0; i < n; i++,ta++,tb++) tval += (*ta) * (*tb);
      }
    }
  return(tval);
  }
