
/* $Log$
 * Revision 1.1  1993/12/29 12:53:29  etseidl
 * Initial revision
 *
 * Revision 1.4  1992/08/12  11:18:39  seidl
 * kludge for matrices with a row=0.0
 *
 * Revision 1.3  1992/06/17  22:10:10  jannsen
 * clean up for saber-c
 *
 * Revision 1.2  1992/04/06  12:47:24  seidl
 * include math.h
 *
 * Revision 1.1.1.1  1992/03/17  16:35:05  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  16:35:04  seidl
 * Initial revision
 *
 * Revision 1.1  1992/01/06  11:34:27  seidl
 * Initial revision
 * */

static char *rcsid = "$Id$";

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <tmpl.h>
#include "matrix.h"

#include "flin.gbl"
#include "flin.lcl"

/* solves linear equations a * x = b
 * on input "a" contains coefficients, and "a" contains
 * known vectors.  "in" is the dimension of a(in*in)
 * "im" is the number of b vectors
 *
 * math_lin returns determinant of matrix a */

GLOBAL_FUNCTION double
math_lin(a,b,in,im)
double_matrix_t *a;
double_vector_t *b;
int in;
int im;
{
  int i,*indx;
  double det=0;

  indx = (int *) malloc(sizeof(int)*in);
  if (!indx) {fprintf(stderr,"math_lin: !indx\n"); exit(1);}

  ludcmp(a->d,in,indx,&det);

  for (i=0; i < in ; i++) det *= a->d[i][i];

  lubksb(a->d,in,indx,b->d);

  free(indx);

  return(det);
  }

LOCAL_FUNCTION VOID
ludcmp(a,n,indx,d)
double **a;
int n;
int *indx;
double *d;

{
  int i,imax,j,k;
  double big,dum,sum,temp;
  double *vv;

  vv = (double *) malloc(sizeof(double)*n);
  if (!vv) {fprintf(stderr,"ludcmp: !vv\n"); exit(1);}

  *d = 1.0;

  for (i=0; i < n ; i++) {
    big=0.0;
    for (j=0; j < n; j++) if ((temp=fabs(a[i][j])) > big) big=temp;
#if 0
    if (big == 0.0) {
      *d = 0.0;
      return;
      }
#else
    if(big==0.0) big=1.0e-16;
#endif
    vv[i] = 1.0/big;
    }

  for (j=0; j < n ; j++) {
    for (i=0; i < j ; i++) {
      sum = a[i][j];
      for (k=0; k < i ; k++) sum -= a[i][k]*a[k][j];
      a[i][j] = sum;
      }

    big = 0.0;
    for (i=j ; i < n ; i++) {
      sum=a[i][j];
      for (k=0; k < j ; k++) sum -= a[i][k]*a[k][j];
      a[i][j] = sum;
      if ((dum=vv[i]*fabs(sum)) >= big) {
        big = dum;
        imax = i;
        }
      }

    if (j != imax) {
      for (k=0; k < n; k++) {
        dum=a[imax][k];
        a[imax][k]=a[j][k];
        a[j][k]=dum;
        }
      *d = -(*d);
      vv[imax]=vv[j];
      }

    indx[j]=imax;
    if (a[j][j] == 0.0) a[j][j] = 1.0e-20;
    if (j != n-1) {
      dum = 1.0/a[j][j];
      for (i=j+1; i < n ; i++) a[i][j] *= dum;
      }
    }
  free(vv);
  }

LOCAL_FUNCTION VOID
lubksb(a,n,indx,b)
double **a;
int n;
int *indx;
double *b;
{
  int i,ii,ip,j;
  int t=0;
  double sum;

  for (i=0; i < n ; i++) {
    ip = indx[i];
    sum = b[ip];
    b[ip]=b[i];

    if(t) for (j=ii; j <= i-1 ; j++) sum -= a[i][j]*b[j];
    else if(sum) {
       ii=i;
       t++;
       }

    b[i]=sum;
    }

  for (i=n-1; i >= 0 ; i--) {
    sum = b[i];
    for (j=i+1; j < n ; j++) sum -= a[i][j]*b[j];
    b[i] = sum/a[i][i];
    }
  }
