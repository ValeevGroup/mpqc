
/* $Log$
 * Revision 1.2  1994/08/26 17:57:40  etseidl
 * get rid of rcs id's and fix a few warnings
 *
 * Revision 1.1.1.1  1993/12/29  12:53:29  etseidl
 * SC source tree 0.1
 *
 * Revision 1.2  1992/06/17  22:10:07  jannsen
 * clean up for saber-c
 *
 * Revision 1.1.1.1  1992/03/17  16:35:03  seidl
 * DOE-NIH Quantum Chemistry Library 0.0
 *
 * Revision 1.1  1992/03/17  16:35:01  seidl
 * Initial revision
 *
 * Revision 1.1  1991/12/20  16:01:16  seidl
 * Initial revision
 *
 * Revision 1.2  91/12/03  00:24:50  etseidl
 * change alloc_matrix to matrixallc, free_matrix to matrixfree
 * 
 * Revision 1.1  1991/12/02  23:40:01  etseidl
 * Initial revision
 * */

#include <stdio.h>
#include <math.h>
#include <tmpl.h>
#include "matrix.h"
#include "matrixallc.h"
#include "matrixfree.h"

#include "diag.gbl"
#include "diag.lcl"


/* translation into c of a translation into FORTRAN77 of the EISPACK */
/* matrix diagonalization routines */

/* array is a square double_matrix to be diagonalized
 * e_vals is a double_vector at least array->n1 in size
 *        on return e_vals holds the eigenvalues of array
 * e_vecs is a double_matrix at least array->n1 x array->n1 in size
 *        on return e_vecs holds the eigenvectors of array if matz==1
 * matz is a boolean, if matz=0, only eigenvalues are calculated
 * toler is a double, how close to zero you want the off-diagonal elements
 * sort is a boolean, if sort=1, the eigenvectors are sorted from smallest to
 *            largest
 */

GLOBAL_FUNCTION VOID 
math_diag_dm(array,e_vals,e_vecs,matz,toler,sort)
double_matrix_t *array;
double_vector_t *e_vals;
double_matrix_t *e_vecs;
int matz;
double toler;
int sort;
{
  int i, j;
  int n;
  int errcod;
  double_vector_t fv1;
  double_matrix_t temp;

  if(array->n1 != array->n2) {
    fprintf(stderr,"diag_dm:\n");
    fprintf(stderr,"array is not a square matrix\n");
    exit(-1);
    }
  if(array->n1 > e_vals->n) {
    fprintf(stderr,"diag_dm:\n");
    fprintf(stderr,"e_vals vector too small: want %d, is %d\n",
                                            array->n1,e_vals->n);
    exit(-1);
    }
  if(array->n1 > e_vecs->n1 || array->n2 > e_vecs->n2) {
    fprintf(stderr,"diag_dm:\n");
    fprintf(stderr,"e_vecs vector too small: want %dx%d, is %dx%d\n",
                         array->n1,array->n2,e_vecs->n1,e_vecs->n2);
    exit(-1);
    }

  n = array->n1;

  errcod = allocbn_double_vector(&fv1,"n",n);
  if(errcod != 0) {
    fprintf(stderr,"diag_dm:\n");
    fprintf(stderr,"could not allocate memory for temp array fv1\n");
    exit(-1);
    }

  errcod = allocbn_double_matrix(&temp,"n1 n2",n,n);
  if(errcod != 0) {
    fprintf(stderr,"diag_dm:\n");
    fprintf(stderr,"could not allocate memory for temp array temp\n");
    exit(-1);
    }

  for(i=0; i < n; i++)
    for(j=0; j < n; j++)
      e_vecs->d[i][j] = array->d[i][j];

  tred2(n,e_vecs,e_vals,&fv1,matz);

  for(i=0; i < n; i++)
    for(j=0; j < n; j++)
      temp.d[i][j]=e_vecs->d[j][i];

  tqli(n,e_vals,&temp,&fv1,matz,toler);

  for(i=0; i < n; i++)
    for(j=0; j < n; j++)
      e_vecs->d[i][j]=temp.d[j][i];

  if(sort) eigsort(n,e_vals,e_vecs);

  free_double_vector(&fv1);
  free_double_matrix(&temp);
  }
            
/* math_diag_dv is the same thing as math_diag_dm, except array is a 
 * double_vector which holds the lower triangle of a square, symmetric
 * matrix */

GLOBAL_FUNCTION VOID 
math_diag_dv(array,e_vals,e_vecs,matz,toler,sort)
double_vector_t *array;
double_vector_t *e_vals;
double_matrix_t *e_vecs;
int matz;
double toler;
int sort;
{
  int i,j,ij;
  int errcod;
  int n,nt;
  double_vector_t fv1;
  double_matrix_t temp;

  nt = array->n;
  n = (((int) sqrt(1.0+8.0*(double)nt))-1)/2;

  if(n > e_vals->n) {
    fprintf(stderr,"diag_dv:\n");
    fprintf(stderr,"e_vals vector too small: want %d, is %d\n",n,e_vals->n);
    exit(-1);
    }
  if(n > e_vecs->n1 || n > e_vecs->n2) {
    fprintf(stderr,"diag_dv:\n");
    fprintf(stderr,"e_vecs vector too small: want %dx%d, is %dx%d\n",
                                                 n,n,e_vecs->n1,e_vecs->n2);
    exit(-1);
    }

  errcod = allocbn_double_vector(&fv1,"n",n);
  if(errcod != 0) {
    fprintf(stderr,"diag_dv:\n");
    fprintf(stderr,"could not allocate memory for temp array fv1\n");
    exit(-1);
    }

  errcod = allocbn_double_matrix(&temp,"n1 n2",n,n);
  if(errcod != 0) {
    fprintf(stderr,"diag_dv:\n");
    fprintf(stderr,"could not allocate memory for temp array temp\n");
    exit(-1);
    }


  for(i=ij=0; i < n; i++)
    for(j=0; j <= i; j++,ij++)
      e_vecs->d[i][j] = e_vecs->d[j][i] = array->d[ij];

  tred2(n,e_vecs,e_vals,&fv1,matz);

  for (i=0; i < n; i++)
    for (j=0; j < n; j++)
      temp.d[i][j]=e_vecs->d[j][i];
            
  tqli(n,e_vals,&temp,&fv1,matz,toler);

  for(i=0; i < n; i++)
    for(j=0; j < n; j++)
      e_vecs->d[i][j]=temp.d[j][i];

  if(sort) eigsort(n,e_vals,e_vecs);

  free_double_vector(&fv1);
  free_double_matrix(&temp);
  }
            

/* converts symmetric matrix to a tridagonal form for use in tqli
 * if matz = 0, only find eigenvalues, else find both eigenvalues and
 * eigenvectors */

#ifdef DSIGN
#undef DSIGN
#endif
#define DSIGN(a,b) ((b) >= 0.0) ? (fabs(a)) : (-fabs(a))

LOCAL_FUNCTION VOID
tred2(n,aa,dd,ee,matz)
int n;
double_matrix_t *aa;
double_vector_t *dd;
double_vector_t *ee;
int matz;
{
  int i,j,k,l;
  double f,g,h,hh,scale,scale_inv,h_inv;
  double **a,*d,*e;

  a = aa->d;
  d = dd->d;
  e = ee->d;
  
  if (n == 1) return;

  for(i=n-1; i > 0; i--) {
    l = i-1;
    h = 0.0;
    scale = 0.0;
    if(l) {
      for(k=0; k <= l; k++) scale += fabs(a[i][k]);
      if (scale == 0.0) e[i] = a[i][l];
      else {
        scale_inv=1.0/scale;
        for (k=0; k <= l; k++) {
          a[i][k] *= scale_inv;
          h += a[i][k]*a[i][k];
          }
        f=a[i][l];
        g= -(DSIGN(sqrt(h),f));
        e[i] = scale*g;
        h -= f*g;
        a[i][l] = f-g;
        f = 0.0;
        h_inv=1.0/h;
        for (j=0; j <= l; j++) {
          if (matz) a[j][i] = a[i][j]*h_inv;
          g = 0.0;
          for (k=0; k <= j; k++) g += a[j][k]*a[i][k];
          if (l > j) for (k=j+1; k <= l; k++) g += a[k][j]*a[i][k];
          e[j] = g*h_inv;
          f += e[j]*a[i][j];
          }
        hh = f/(h+h);
        for (j=0; j <= l; j++) {
          f = a[i][j];
          g = e[j] - hh*f;
          e[j] = g;
          for (k=0; k <= j; k++) a[j][k] -= (f*e[k] + g*a[i][k]);
          }
        }
      }
    else {
      e[i] = a[i][l];
      }
    d[i] = h;
    }
  if(matz) d[0] = 0.0;
  e[0] = 0.0;

  for(i=0; i < n; i++) {
    l = i-1;
    if (matz) {
      if(d[i]) {
        for(j=0; j <= l; j++) {
          g = 0.0;
          for(k=0; k <= l; k++) g += a[i][k]*a[k][j];
          for(k=0; k <= l; k++) a[k][j] -= g*a[k][i];
          }
        }
      }
    d[i] = a[i][i];
    if(matz) {
      a[i][i] = 1.0;
      if(l >= 0) for (j=0; j<= l; j++) a[i][j] = a[j][i] = 0.0;
      }
    }
  }

/* diagonalizes tridiagonal matrix output by tred2 */
/* gives only eigenvalues if matz = 0, both eigenvalues and eigenvectors */
/* if matz = 1 */

LOCAL_FUNCTION VOID
tqli(n,dd,zz,ee,matz,toler)
int n;
double_vector_t *dd;
double_matrix_t *zz;
double_vector_t *ee;
int matz;
double toler;
{
  register int k;
  int i,l,m,iter;
  double g,r,s,c,p,f,b;
  double *d, **z, *e;
  double azi;

  d = dd->d;
  e = ee->d;
  z = zz->d;

  f=0.0;
  if (n == 1) {
    d[0]=z[0][0];
    z[0][0] = 1.0;
    return;
    }

  for (i=1; i < n ; i++) e[i-1] = e[i];
  e[n-1] = 0.0;
  for (l=0; l < n; l++) {
    iter = 0;
L1:
    for (m=l; m < n-1;m++) if (fabs(e[m]) < toler) goto L2;
    m=n-1;
L2:
    if (m != l) {
      if (iter++ == 35) {
        fprintf (stderr,"tqli not converging\n");
        continue;
        }
      else if (iter > 30) {
        fprintf(stderr,"tqli: iter=%d, l=%d, m=%d, e[m]=% f\n",
                iter,l,m,e[m]);
        }

      g = (d[l+1]-d[l])/(2.0*e[l]);
      r = sqrt(g*g + 1.0);
      g = d[m] - d[l] + e[l]/((g + DSIGN(r,g)));
      s=1.0;
      c=1.0;
      p=0.0;
      for (i=m-1; i >= l; i--) {
        f = s*e[i];
        b = c*e[i];
        if (fabs(f) >= fabs(g)) {
          c = g/f;
          r = sqrt(c*c + 1.0);
          e[i+1] = f*r;
          s=1.0/r;
          c *= s;
          }
        else {
          s = f/g;
          r = sqrt(s*s + 1.0);
          e[i+1] = g*r;
          c = 1.0/r;
          s *= c;
          }
        g = d[i+1] - p;
        r = (d[i]-g)*s + 2.0*c*b;
        p = s*r;
        d[i+1] = g+p;
        g = c*r-b;

        if (matz) {
          double *zi = z[i];
          double *zi1 = z[i+1];
          for (k=n; k ; k--,zi++,zi1++) {
            azi = *zi;
            f = *zi1;
            *zi1 = azi*s + c*f;
            *zi = azi*c - s*f;
            }
          }
        }

      d[l] -= p;
      e[l] = g;
      e[m] = 0.0;
      goto L1;
      }
    }
  }


LOCAL_FUNCTION VOID
eigsort(n,d,v)
int n;
double_vector_t *d;
double_matrix_t *v;
{
  int i,j,k;
  double p;

  for(i=0; i < n-1 ; i++) {
    k=i;
    p=d->d[i];
    for(j=i+1; j < n; j++) {
      if(d->d[j] < p) {
        k=j;
        p=d->d[j];
        }
      }
    if(k != i) {
      d->d[k]=d->d[i];
      d->d[i]=p;
      for(j=0; j < n; j++) {
        p=v->d[j][i];
        v->d[j][i]=v->d[j][k];
        v->d[j][k]=p;
        }
      }
    }
  }
