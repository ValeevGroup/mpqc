
/* This function will perform the diis extrapolation a la Pulay.
 * For more information read Hamilton and Pulay, J.Chem.Phys 84 (1986), 5728.
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <tmpl.h>
#include <math.h>
#include <math/dmt/libdmt.h>
#include <util/misc/libmisc.h>
#include <math/array/math_lib.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <chemistry/qc/dmtscf/scf.h>

#include <chemistry/qc/dmtscf/scf_diis.gbl>
#include <chemistry/qc/dmtscf/scf_diis.lcl>

/* keep around the bmatrix so we don't have to recompute the whole thing
 * every iteration
 */
static double_vector_t btemp;
static double_matrix_t bold,bmat;

/* these store the extrapolated fock matrices */
static dmt_matrix LFockC,LFockO;

/* this struct is for holding the last ndiis fock matrices and error vectors */
static struct diis_mats {
  dmt_matrix fock_c;
  dmt_matrix fock_o;
  dmt_matrix error;
} *diism,dtemp;

/***************************************************************************
 *
 * perform the diis extrapolation.
 *
 * input:
 *   scf_info = scf struct
 *   Fock     = the AO Fock matrix from this iteration (scattered)
 *   FockO    = the AO open-shell Fock matrix from this iteration (scattered)
 *   Scf_Vec  = the scf vector from this iteration (column dist)
 *   occ_num  = double pointer to a vector containing MO occupations
 *   iter     = the current iteration
 *   Scr1, Scr2, Scr3 = scratch dmt matrices (column distributed)
 *
 * on return:
 *   Fock and FockO contain the extrapolated Fock matrices in the AO basis
 *
 * return 0 on success, -1 on failure
 */

GLOBAL_FUNCTION int
scf_diis(scf_info, Fock, FockO, Scf_Vec, occ_num, iter, Scr1, Scr2, Scr3)
scf_struct_t *scf_info;
dmt_matrix Fock;
dmt_matrix FockO;
dmt_matrix Scf_Vec;
double *occ_num;
int iter;
dmt_matrix Scr1;
dmt_matrix Scr2;
dmt_matrix Scr3;
{
  int i,j,k;
  int ib,jb,isz,jsz,ist,jst;
  int nl,nlocal=dmt_nlocal(Fock);
  int try = 0;
  int last = iter;
  int col = iter+2;
  int iopen = scf_info->iopen;
  double occi, occj;
  double norm, determ;
  double *Fblk,*FOblk,*Eblk;
  double scale;
  struct diis_mats *d;

  assert(dmt_distribution(Scf_Vec) == COLUMNS);
  assert(dmt_distribution(Scr1) == COLUMNS);
  assert(dmt_distribution(Scr2) == COLUMNS);
  assert(dmt_distribution(Scr3) == COLUMNS);
  assert(dmt_distribution(Fock) == SCATTERED);
  if (scf_info->iopen) assert(dmt_distribution(FockO) == SCATTERED);

  iter++;

 /* Allocate memory for the fock and error matrices for the last "ndiis"
  * iterations.  This is going to lead to memory problems eventually,
  * so it will probably become necessary to transfer all this to disk
  */
  if (!diism) {
    if (scf_init_diis(scf_info) < 0) {
      fprintf(stderr,"scf_diis:  scf_init_diis failed\n");
      return -1;
    }
  }

  scale = 1.0 + scf_info->diisdamp;

  if (iter > scf_info->ndiis) {
    last = scf_info->ndiis-1;
    col = scf_info->ndiis+1;
    dtemp = diism[0];
    for (i=0; i < last ; i++) diism[i] = diism[i+1];
    diism[last] = dtemp;
  }
      
 /* save ao fock matrices in fock_save */
   
  d = &diism[last];
  dmt_copy(Fock,d->fock_c);
  if (iopen) dmt_copy(FockO,d->fock_o);

 /* Form error matrix in mo basis
  * This error matrix is just the off-diagonal blocks of the fock matrix 
  * rather than FDS-SDF recommended in the 1982 DIIS paper by Pulay.
  */

  dmt_mult(Fock,Scf_Vec,Scr1);
  dmt_mult(Scf_Vec,Scr1,Scr2);
  dmt_copy(Scr2,LFockC);

  if (iopen) {
    dmt_mult(FockO,Scf_Vec,Scr1);
    dmt_mult(Scf_Vec,Scr1,Scr2);
    dmt_copy(Scr2,LFockO);
  }

  for (nl=0; nl < nlocal ; nl++) {
    dmt_get_block(LFockC,nl,&ib,&jb,&Fblk);
    dmt_get_block(d->error,nl,&ib,&jb,&Eblk);
    if (iopen) dmt_get_block(LFockO,nl,&ib,&jb,&FOblk);
    dmt_describe_block(LFockC,ib,&ist,&isz);
    dmt_describe_block(LFockC,jb,&jst,&jsz);

    for (i=0; i < isz; i++) {
      occi = occ_num[ist+i];

      for (j=0; j < jsz ; j++) {
        occj = occ_num[jst+j];

        if (!iopen) {
          if ((occi+occj)==2.0)
            Eblk[i*jsz+j] = Fblk[i*jsz+j];
          else
            Eblk[i*jsz+j] = 0.0;
        } else if (!scf_info->twocon) {
          if(occi == occj)
            Eblk[i*jsz+j] = 0.0;
          else if(occi && occj)
            Eblk[i*jsz+j] = 2.0*Fblk[i*jsz+j]-FOblk[i*jsz+j];
          else if(occi==2.0 || occj==2.0)
            Eblk[i*jsz+j] = Fblk[i*jsz+j];
          else
            Eblk[i*jsz+j] = FOblk[i*jsz+j];
        } else {
          /* twocon code here eventually */
          fprintf(stderr,"scf_diis: error: twocon code not ready\n");
          return -1;
        }
      }
    }
  }

  /* transform error matrix into ao basis */
      
  dmt_copy(Scf_Vec,Scr3);
  dmt_transpose(Scr3);
  dmt_mult(d->error,Scr3,Scr1);
  dmt_mult(Scr3,Scr1,Scr2);
  dmt_copy(Scr2,d->error);

  scf_info->diis_er = dmt_max_abs(d->error);
               
 /* then set up B matrix, where B(i,j) = <ei|ej> */

 /* move bold(i+1,j+1) to bold(i,j) */
  if (iter > scf_info->ndiis) {
    for (i=0; i < last ; i++) {
      for (j=0; j <= i ; j++) {
        bold.d[i][j]=bold.d[j][i]=bold.d[i+1][j+1];
      }
    }
  }

 /* and set the current rows of bold */
  for (i=0; i <= last ; i++)
    bold.d[i][last]=bold.d[last][i] = 
      dmt_adotb(diism[i].error,diism[last].error);

  bmat.d[0][0] = 0.0;
  btemp.d[0] = -1.0;

  if (bold.d[0][0] > 1.e-10) {
    norm = 1.0/bold.d[0][0];
  } else {
    norm = 1.0;
  }

  for (i=1; i <= last+1 ; i++) {
    bmat.d[i][0]=bmat.d[0][i] = -1.0;
    btemp.d[i] = 0.0;
    for (j=1; j <= i ; j++) {
      bmat.d[i][j]=bmat.d[j][i] = bold.d[i-1][j-1]*norm;
      if (i==j) bmat.d[i][j] *= scale;
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

      if (bold.d[try][try] > 1.e-10) {
        norm=1.0/bold.d[try][try];
      } else {
        norm = 1.0;
      }

      for (i=1; i <= scf_info->ndiis-try ; i++) {
        bmat.d[i][0]=bmat.d[0][i] = -1.0;
        for (j=1; j <= i ; j++) {
          bmat.d[i][j]=bmat.d[j][i]=bold.d[i+try-1][j+try-1]*norm;
          if (i==j) bmat.d[i][j] *= scale;
        }
        btemp.d[i] = 0.0;
      }

      determ = math_lin(&bmat,&btemp,col,1);
    }

    if (fabs(determ) < 10.0e-20) {
      fprintf(stderr,"scf_diis:  try %d no good\n",try);
      return -1;
    }

    if (iter >= scf_info->it_diis) {
      int kk=1;

      dmt_fill(Fock,0.0);
      if (iopen) dmt_fill(FockO,0.0);

      for (k=try; k < last+1 ; k++) {
        dmt_sum_scaled(diism[k].fock_c,btemp.d[kk],Fock);
        if(iopen) dmt_sum_scaled(diism[k].fock_o,btemp.d[kk],FockO);
        kk++;
      }
    }
  }

  return(0);
}

/*************************************************************************
 *
 * free all static arrays and matrices 
 *
 */

GLOBAL_FUNCTION void
scf_done_diis(scf_info)
scf_struct_t *scf_info;
{
  int i;

  dmt_free(LFockC);
  if (scf_info->iopen) dmt_free(LFockO);

  free_double_matrix(&bmat);
  free_double_matrix(&bold);
  free_double_vector(&btemp);

  for (i=0; i < scf_info->ndiis ; i++) {
    dmt_free(diism[i].fock_c);
    dmt_free(diism[i].error);
    if (scf_info->iopen) dmt_free(diism[i].fock_o);
  }

  free(diism);
  diism=NULL;
}

/****************************************************************************
 *
 * given an scf struct, initialize the diis matrices
 *
 * return 0 on success, -1 on failure
 */

GLOBAL_FUNCTION int
scf_init_diis(scf_info)
scf_struct_t *scf_info;
{
  int i;
  int ndiis = scf_info->ndiis;
  int dim = ndiis+1;
  char arg[80];

  if ((allocbn_double_matrix(&bmat,"n1 n2",dim,dim) != 0) ||
      (allocbn_double_matrix(&bold,"n1 n2",ndiis,ndiis) != 0) ||
      (allocbn_double_vector(&btemp,"n",dim) != 0)) {
    fprintf(stderr,"scf_diis_init:  alloc of bmat, bold, and btemp failed\n");
    return -1;
  }

  LFockC = dmt_create("scf_diis: diis fockc",scf_info->nbfao,SCATTERED);

  if (scf_info->iopen) {
    LFockO = dmt_create("scf_diis: diis focko",scf_info->nbfao,SCATTERED);
  } else {
    LFockO = dmt_nil();
  }

  diism = (struct diis_mats *) malloc(sizeof(struct diis_mats)*ndiis);
  if (!diism) {
    fprintf(stderr,"scf_diis_init: could not malloc diism\n");
    return -1;
  }

  for(i=0; i < scf_info->ndiis ; i++) {
    sprintf(arg,"scf_diis: error matrix %d",i);
    diism[i].error = dmt_create(arg,scf_info->nbfao,SCATTERED);

    sprintf(arg,"scf_diis: closed fock matrix %d",i);
    diism[i].fock_c = dmt_create(arg,scf_info->nbfao,SCATTERED);

    if (scf_info->iopen) {
      sprintf(arg,"scf_diis: open fock matrix %d",i);
      diism[i].fock_o = dmt_create(arg,scf_info->nbfao,SCATTERED);
    }
  }

  return 0;
}
