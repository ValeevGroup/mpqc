
/* This function will perform the diis extrapolation a la Pulay.
 * For more information read Hamilton and Pulay, J.Chem.Phys 84 (1986), 5728.
 */

/* $Log$
 * Revision 1.3  1994/06/08 01:15:03  cljanss
 * Many changes.  These include: newmat7 and nihmatrix -> scmat
 * and mpqcic -> MPSCF and updated optimize stuff.
 *
 * Revision 1.2  1993/12/30  13:31:15  etseidl
 * merge in clj changes, do global sum of exchange energy in scf_ex.c
 *
 * Revision 1.5  1992/06/23  20:04:21  seidl
 * change dmt matrices to uppercase,
 * get rid of unnecessary matrice
 *
 * Revision 1.4  1992/06/17  21:54:07  jannsen
 * cleaned up for saber-c
 *
 * Revision 1.3  1992/05/26  20:17:36  jannsen
 * use mtype_get to get message types for global operations
 * check results of memory allocations
 *
 * Revision 1.2  1992/04/06  12:35:28  seidl
 * set diism to null on exit from libscf
 *
 * Revision 1.1.1.1  1992/03/17  16:25:56  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  16:25:55  seidl
 * Initial revision
 *
 * Revision 1.3  1992/02/21  14:48:11  seidl
 * fix small bug in creation of FOCKC and FOCKO
 *
 * Revision 1.2  1992/02/13  00:41:00  seidl
 * add flag for freeing memory when done
 * return mo fock matrices when done
 *
 * Revision 1.1  1992/02/04  23:48:08  seidl
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
#include <stdlib.h>
#include <tmpl.h>
#include <math.h>
#include <math/dmt/libdmt.h>
#include <util/misc/libmisc.h>
#include <math/array/math_lib.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <util/keyval/ipv2c.h>
#include "scf.h"

#include "scf_diis.gbl"
#include "scf_diis.lcl"

static double_vector_t btemp;
static double_matrix_t bold,bmat;

static dmt_matrix FOCKC,FOCKO;

static struct diis_mats {
  dmt_matrix fock_c;
  dmt_matrix fock_o;
  dmt_matrix error;
  } *diism,dtemp;

GLOBAL_FUNCTION int
scf_diis(_scf_info,FOCK,FOCKO,SCF_VEC,_occ_num,iter,
 SCR1,SCR2,SCR3,done,_outfile)
scf_struct_t *_scf_info;
dmt_matrix FOCK;
dmt_matrix FOCKO;
dmt_matrix SCF_VEC;
double_vector_t *_occ_num;
int iter;
dmt_matrix SCR1;
dmt_matrix SCR2;
dmt_matrix SCR3;
int done;
FILE *_outfile;
{
  int i,j,k;
  int ib,jb,isz,jsz,ist,jst;
  int nl,nlocal=dmt_nlocal(FOCK);
  int try = 0;
  int last = iter;
  int col = iter+2;
  int iopen = _scf_info->iopen;
  double occi, occj;
  double norm, determ;
  double *Fblk,*FOblk,*Eblk;
  double scale;
  struct diis_mats *d;

  iter++;

/* if we're done, put the last mo fock matrices back in FOCK and FOCKO
 * and then release all static arrays and matrices 
 */
  if(done) {
    dmt_mult(FOCKC,SCF_VEC,SCR1);
    dmt_mult(SCF_VEC,SCR1,SCR2);
    dmt_copy(SCR2,FOCK);
    if(iopen) {
      dmt_mult(FOCKO,SCF_VEC,SCR1);
      dmt_mult(SCF_VEC,SCR1,SCR2);
      dmt_copy(SCR2,FOCKO);
      }

    dmt_free(FOCKC);
    if(iopen) dmt_free(FOCKO);

    free_double_matrix(&bmat);
    free_double_matrix(&bold);
    free_double_vector(&btemp);

    for(i=0; i < _scf_info->ndiis ; i++) {
      dmt_free(diism[i].fock_c);
      dmt_free(diism[i].error);
      if(iopen) dmt_free(diism[i].fock_o);
      }
    free(diism);
    diism=NULL;
    return(0);
    }


/* Allocate memory for the fock and error matrices for the last "ndiis"
 * iterations.  This is going to lead to memory problems eventually,
 * so it will probably become necessary to transfer all this to disk
 */

  if(diism == NULL) {
    char arg[80];

    if (  (allocbn_double_matrix(&bmat,"n1 n2",_scf_info->ndiis+1,
                                 _scf_info->ndiis+1) != IPE_OK)
        ||(allocbn_double_matrix(&bold,"n1 n2",_scf_info->ndiis,
                                 _scf_info->ndiis) != IPE_OK)
        ||(allocbn_double_vector(&btemp,"n",
                                 _scf_info->ndiis+1) != IPE_OK)) {
      printf("alloc of bmat, bold, and btemp failed\n");
      return -1;
      }

    FOCKC = dmt_create("libscfv3: diis fockc",_scf_info->nbfao,SCATTERED);
    if(iopen) 
      FOCKO = dmt_create("libscfv3: diis focko",_scf_info->nbfao,SCATTERED);
    else
      FOCKO = dmt_nil();

    diism = 
      (struct diis_mats *) malloc(sizeof(struct diis_mats)*_scf_info->ndiis);
    check_alloc(diism,"diism");

    for(i=0; i < _scf_info->ndiis ; i++) {
      sprintf(arg,"scf_diis: error matrix %d",i);
      diism[i].error = dmt_create(arg,_scf_info->nbfao,SCATTERED);

      sprintf(arg,"scf_diis: closed fock matrix %d",i);
      diism[i].fock_c = dmt_create(arg,_scf_info->nbfao,SCATTERED);

      if(iopen) {
        sprintf(arg,"scf_diis: open fock matrix %d",i);
        diism[i].fock_o = dmt_create(arg,_scf_info->nbfao,SCATTERED);
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
  dmt_copy(FOCK,d->fock_c);
  if(iopen) dmt_copy(FOCKO,d->fock_o);

 /* Form error matrix in mo basis
  * This error matrix is just the off-diagonal elements of the fock matrix 
  * rather than FDS-SDF recommended in the 1982 DIIS paper by Pulay.
  */

  dmt_mult(FOCK,SCF_VEC,SCR1);
  dmt_mult(SCF_VEC,SCR1,SCR2);
  dmt_copy(SCR2,FOCKC);

  if(iopen) {
    dmt_mult(FOCKO,SCF_VEC,SCR1);
    dmt_mult(SCF_VEC,SCR1,SCR2);
    dmt_copy(SCR2,FOCKO);
    }

  for(nl=0; nl < nlocal ; nl++) {
    dmt_get_block(FOCKC,nl,&ib,&jb,&Fblk);
    dmt_get_block(d->error,nl,&ib,&jb,&Eblk);
    if(iopen) dmt_get_block(FOCKO,nl,&ib,&jb,&FOblk);
    dmt_describe_block(FOCKC,ib,&ist,&isz);
    dmt_describe_block(FOCKC,jb,&jst,&jsz);

    for (i=0; i < isz; i++) {
      occi = _occ_num->d[ist+i];
      for (j=0; j < jsz ; j++) {
        occj = _occ_num->d[jst+j];
        if (!iopen) {
          if ((occi+occj)==2.0)
            Eblk[i*jsz+j] = Fblk[i*jsz+j];
          else
            Eblk[i*jsz+j] = 0.0;
          }
        else if(!_scf_info->twocon) {
          if(occi == occj)
            Eblk[i*jsz+j] = 0.0;
          else if(occi && occj)
            Eblk[i*jsz+j] = 2.0*Fblk[i*jsz+j]-FOblk[i*jsz+j];
          else if(occi==2.0 || occj==2.0)
            Eblk[i*jsz+j] = Fblk[i*jsz+j];
          else
            Eblk[i*jsz+j] = FOblk[i*jsz+j];
          }
        else {
          /* twocon code here eventually */
          printf("error: twocon code not ready\n");
          return -1;
          }
        }
      }
    }

  /* transform error matrix into ao basis */
      
  dmt_copy(SCF_VEC,SCR3);
  dmt_transpose(SCR3);
  dmt_mult(d->error,SCR3,SCR1);
  dmt_mult(SCR3,SCR1,SCR2);
  dmt_copy(SCR2,d->error);

  _scf_info->diis_er = dmt_max_abs(d->error);
               
  /* then set up B matrix, where B(i,j) = <ei|ej> */

  if (iter > _scf_info->ndiis) {
    for (i=0; i < last ; i++) {
      for (j=0; j <= i ; j++) {
        bold.d[i][j]=bold.d[j][i]=bold.d[i+1][j+1];
        }
      }
    }
  for (i=0; i <= last ; i++)
    bold.d[i][last]=bold.d[last][i] = 
      dmt_adotb(diism[i].error,diism[last].error);

  bmat.d[0][0] = 0.0;
  btemp.d[0] = -1.0;
  if (bold.d[0][0] > 1.e-10) norm = 1.0/bold.d[0][0];
  else norm = 1.0;
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
      if (bold.d[try][try] > 1.e-10) norm=1.0/bold.d[try][try];
      else norm = 1.0;
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
      int kk=1;

      dmt_fill(FOCK,0.0);
      if(iopen) dmt_fill(FOCKO,0.0);

      for (k=try; k < last+1 ; k++) {
        dmt_sum_scaled(diism[k].fock_c,btemp.d[kk],FOCK);
        if(iopen) dmt_sum_scaled(diism[k].fock_o,btemp.d[kk],FOCKO);
        kk++;
        }
      }
    }

  dmt_copy(FOCK,FOCKC);
  if(iopen) dmt_copy(FOCKO,FOCKO);

  return(0);
  }
