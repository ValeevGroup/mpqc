
/* 
 * These routines are based on the work of Edward T. Seidl at the
 * National Institutes of Health.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <math/scmat/cmatrix.h>

static void ludcmp(double**, int dim, int*, double*);
static void lubksb(double**, int dim, int*, double*);

static void tred2(int dim,double**,double*,double*,int);
static void tqli(int dim,double*,double**,double*,int,double);
static void eigsort(int dim,double*,double**);

double**
cmat_new_square_matrix(int n)
{
  double* mat = (double*) malloc(sizeof(double)*n*n);
  double** r = (double**) malloc(sizeof(double*)*n);
  cmat_matrix_pointers(r,mat,n,n);
  return r;
}

double**
cmat_new_rect_matrix(int n,int m)
{
  double* mat = (double*) malloc(sizeof(double)*n*m);
  double** r = (double**) malloc(sizeof(double*)*n);
  cmat_matrix_pointers(r,mat,n,m);
  return r;
}

/* this deletes both square and triangular matrices */
void
cmat_delete_matrix(double**m)
{
  free(m[0]);
  free(m);
}

void
cmat_transpose_square_matrix(double**matrix, int n)
{
  int i,j;
  for (i=0; i<n; i++) {
      for (j=0; j<i; j++) {
          double tmp = matrix[i][j];
          matrix[i][j] = matrix[j][i];
          matrix[j][i] = tmp;
        }
    }
}

void
cmat_matrix_pointers(double**ptrs,double*matrix,int nrow, int ncol)
{
  int i;
  for (i=0; i<nrow; i++) ptrs[i] = &matrix[i*ncol];
}

/*
 * a contains pointers to the an area of contiguous storage.
 * Its dimensions are nr by nc.  On exit it will be transposed,
 * however the a vector of double* is itself unchanged.  Another
 * vector is needed to access the storage or a must be updated
 * after this routine is called.
 */
void
cmat_transpose_matrix(double**a, int nr, int nc)
{
  int i,j;
  double* tmpp;
  double* tmp;

  if (nr == nc) {
      cmat_transpose_square_matrix(a,nr);
      return;
    };

  tmp = (double*) malloc(sizeof(double)*nr*nc);

  tmpp = tmp;
  for (i=0; i<nc; i++) {
      for (j=0; j<nr; j++) {
          *tmpp = a[j][i];
          tmpp++;
        }
    }

  memcpy(a[0],tmp,sizeof(double)*nr*nc);

  free(tmp);
}

/* a is symmetric if sym is true */
double
cmat_determ(double** a, int sym, int dim)
{
  int i;
  double det=0;
  int *indx= (int*) malloc(sizeof(int)*dim);

  if (sym) {
    fprintf(stderr,"cmat_determ: can't handle sym=1 yet\n");
    abort();
  }

  ludcmp(a,dim,indx,&det);
  free(indx);

  for (i=0; i < dim; i++) det *= a[i][i];

  return det;
}

/* a is symmetric if sym is true */
double
cmat_solve_lin(double** a, int sym, double* b, int dim)
{
  int i;
  double det=0;
  int *indx= (int*) malloc(sizeof(int)*dim);

  if (sym) {
      fprintf(stderr,"cmat_solve_lin: can't handle sym=1 yet\n");
      abort();
    }

  ludcmp(a,dim,indx,&det);

  for(i=0; i < dim; i++) det *= a[i][i];

  lubksb(a,dim,indx,b);
  free(indx);

  return det;
  }

double 
cmat_invert(double**a, int sym, int dim)
{
  int i,j;
  double det=0;
  double **a_orig;
  int *indx= (int*) malloc(sizeof(int)*dim);
  double **y;
  double* b;

  /* if a is a symmetric matrix, then copy it to a nonsymmetric matrix */
  if (sym) {
      a_orig = a;
      a = cmat_new_square_matrix(dim);
      for (i=0; i<dim; i++) {
          for (j=0; j<=i; j++) {
              a[i][j] = a[j][i] = a_orig[i][j];
            }
        }
    }

  ludcmp(a,dim,indx,&det);

  for(i=0; i < dim; i++) det *= a[i][i];

  b = (double*) malloc(sizeof(double)*dim);

  /* y=a */
  y = cmat_new_square_matrix(dim);
  for (i=0; i<dim; i++) {
      for (j=0; j<dim; j++) y[i][j] = a[i][j];
    }

  for(i=0; i < dim; i++) {
    for(j=0; j < dim; j++) b[j]=0;
    b[i]=1;
    lubksb(y,dim,indx,b);
    for(j=0; j < dim; j++) a[j][i]=b[j];
    }

  free(b);
  cmat_delete_matrix(y);
  free(indx);

  /* copy a back to the original a, if necessary */
  if (sym) {
      int i,j;
      for (i=0; i<dim; i++) {
          for (j=0; j<=i; j++) {
              a_orig[i][j] = 0.5*(a[i][j]+a[j][i]);
            }
        }
      cmat_delete_matrix(a);
    }
  return det;
  }

static void
ludcmp(double** a, int n, int *indx, double *d)
{
  int i,imax,j,k;
  double big,dum,sum,temp;

  double* vv = (double*) malloc(sizeof(double)*n);

  *d = 1.0;

  for (i=0; i < n ; i++) {
    big=0.0;
    for (j=0; j < n; j++) if ((temp=fabs(a[i][j])) > big) big=temp;
#if 0
    if (big == 0.0) {
      *d = 0.0;
      free(vv);
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

static void
lubksb(double** a, int n, int *indx, double* b)
{
  int i,ii,ip,j;
  int t=0;
  double sum;

  for (i=0; i < n ; i++) {
    ip = indx[i];
    sum = b[ip];
    b[ip]=b[i];

    if(t) {
      for (j=ii; j <= i-1 ; j++)
        sum -= a[i][j]*b[j];
      }
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

/*
 * This does c(t) (+)= a(t) * b(t), where the (t) means the transpose
 * of the matrix can be optionally used and the (+) means that accumulation
 * is optional.  The dimensions of the matrices is as follows:
 * a(nr,nl) (if ta then a(nl,nr))
 * b(nl,nc) (if tb then b(nc,nl))
 * c(nr,nc) (if tc then c(nc,nr))
 */
void
cmat_mxm(double** a, int ta, double** b, int tb, double** c, int tc,
         int nr, int nl, int nc, int add)
{
  int odd_nr,odd_nc;
  int i,j,k;
  double t00,t01,t10,t11;
  double *att,*bt;
  double *at1,*bt1;
  double** old_a = 0;
  double** old_b = 0;

  odd_nr = (nr)%2;
  odd_nc = (nc)%2;

  if(ta) {
      cmat_transpose_matrix(a,nl,nr);
      if (nr > nl) {
          old_a = a;
          a = (double**) malloc(nr*sizeof(double*));
          a[0] = old_a[0];
        }
      cmat_matrix_pointers(a,a[0],nr,nl);
    }
  if(!tb) {
      cmat_transpose_matrix(b,nl,nc);
      if (nc > nl) {
          old_b = b;
          b = (double**) malloc(nc*sizeof(double*));
          b[0] = old_b[0];
        }
      cmat_matrix_pointers(b,b[0],nc,nl);
    }

  for(j=0; j < nc-1 ; j+=2) {
    for(i=0; i < nr-1 ; i+=2) {
      att=a[i]; bt=b[j];
      at1=a[i+1]; bt1=b[j+1];
      if(add) {
        if(tc) {
          t00 = c[j][i];
          t01 = c[j+1][i];
          t10 = c[j][i+1];
          t11 = c[j+1][i+1];
          }
        else {
          t00 = c[i][j];
          t01 = c[i][j+1];
          t10 = c[i+1][j];
          t11 = c[i+1][j+1];
          }
        }
      else
        t00=t01=t10=t11=0.0;
      for(k=nl; k ; k--,att++,bt++,at1++,bt1++) {
        t00 += *att * *bt;
        t01 += *att * *bt1;
        t10 += *at1 * *bt;
        t11 += *at1 * *bt1;
        }
      if(tc) {
        c[j][i]=t00;
        c[j+1][i]=t01;
        c[j][i+1]=t10;
        c[j+1][i+1]=t11;
        }
      else {
        c[i][j]=t00;
        c[i][j+1]=t01;
        c[i+1][j]=t10;
        c[i+1][j+1]=t11;
        }
      }
    if(odd_nr) {
      att=a[i]; bt=b[j];
      bt1=b[j+1];
      if(add) {
        if(tc) {
          t00 = c[j][i];
          t01 = c[j+1][i];
          }
        else {
          t00 = c[i][j];
          t01 = c[i][j+1];
          }
        }
      else t00=t01=0.0;
      for(k= nl; k ; k--,att++,bt++,bt1++) {
        t00 += *att * *bt;
        t01 += *att * *bt1;
        }
      if(tc) {
        c[j][i]=t00;
        c[j+1][i]=t01;
        }
      else {
        c[i][j]=t00;
        c[i][j+1]=t01;
        }
      }
    }
  if(odd_nc) {
    for(i=0; i < nr-1 ; i+=2) {
      att=a[i]; bt=b[j];
      at1=a[i+1];
      if(add) {
        if(tc) {
          t00 = c[j][i];
          t10 = c[j][i+1];
          }
        else {
          t00 = c[i][j];
          t10 = c[i+1][j];
          }
        }
      else t00=t10=0.0;
      for(k= nl; k ; k--,att++,bt++,at1++) {
        t00 += *att * *bt;
        t10 += *at1 * *bt;
        }
      if(tc) {
        c[j][i]=t00;
        c[j][i+1]=t10;
        }
      else {
        c[i][j]=t00;
        c[i+1][j]=t10;
        }
      }
    if(odd_nr) {
      att=a[i]; bt=b[j];
      if(add) t00 = (tc) ? c[j][i] : c[i][j];
      else t00=0.0;
      for(k=nl; k ; k--,att++,bt++) t00 += *att * *bt;
      if(tc) c[j][i]=t00;
      else c[i][j]=t00;
      }
    }

  if(ta) {
      cmat_transpose_matrix(a,nr,nl);
      if (old_a) {
          free(a);
          a = old_a;
        }
      cmat_matrix_pointers(a,a[0],nr,nl);
    }
  if(!tb) {
      cmat_transpose_matrix(b,nc,nl);
      if (old_b) {
          free(b);
          b = old_b;
        }
      cmat_matrix_pointers(b,b[0],nl,nc);
    }
  }

/*
 * a is symmetric (na,na) in a triangular storage format
 * b is rectangular (na,nb)
 * a (+)= b * transpose(b) (+= if add)
 */
void
cmat_symmetric_mxm(double**a,int na, /* a is (na,na) */
                   double**b,int nb, /* b is (na,nb) */
                   int add)
{
  int i,j,k;
  for (i=0; i<na; i++) {
      double*ai=a[i];
      for (j=0; j<=i; j++) {
          double*bi=b[i];
          double*bj=b[j];
          double tmp;
          if (add) tmp = 0.0;
          else tmp = ai[j];
          for (k=nb; k; k--,bi++,bj++) {
              tmp += *bi * *bj;
            }
          ai[j] = tmp;
        }
    }
}

/*
 * a is symmetric (na,na) in a triangular storage format
 * b is symmetric (nb,nb) in a triangular storage format
 * a (+)= c * b * transpose(c) (+= if add)
 */
void
cmat_transform_symmetric_matrix(double**a,int na, /* a is (na,na) */
                                double**b,int nb, /* b is (nb,nb) */
                                double**c,        /* c is (na,nb) */
                                int add)
{
  int i,j,k;
  double**t;
  double* brow;

  /* create a temporary matrix, t */
  t = cmat_new_rect_matrix(na,nb);

  /* t = transpose(b * transpose(c)) */
  brow = (double*) malloc(sizeof(double)*nb);
  for (i=0; i<nb; i++) {
      for (k=0; k<=i; k++) brow[k] = b[i][k];
      for (   ; k<nb; k++) brow[k] = b[k][i];
      for (j=0; j<na; j++) {
          double*bi = brow;
          double*cj = c[j];
          double tmp = 0.0;
          for (k=nb; k; k--,bi++,cj++) tmp += *bi * *cj;
          t[j][i] = tmp;
        }
    }
  free(brow);

  /* a = c * transpose(t) */
  for (i=0; i<na; i++) {
      for (j=0; j<=i; j++) {
          double*ci = c[i];
          double*tj = t[j];
          double tmp;
          if (add) tmp = a[i][j];
          else tmp = 0.0;
          for (k=nb; k; k--,ci++,tj++) tmp += *ci * *tj;
          a[i][j] = tmp;
        }
    }

  /* delete the temporary */
  cmat_delete_matrix(t);
}

/*
 * a is symmetric (na,na) in a triangular storage format
 * b is diagonal (nb,nb) in a vector storage format
 * a (+)= c * b * transpose(c) (+= if add)
 */
void
cmat_transform_diagonal_matrix(double**a,int na, /* a is (na,na) */
                               double*b,int nb,  /* b is (nb,nb) */
                               double**c,        /* c is (na,nb) */
                               int add)
{
  int i,j,k;
  double t;

  for (i=0; i < na; i++) {
    for (j=0; j <= i; j++) {
      t=0;
      for (k=0; k < nb; k++)
        t += c[i][k] * c[j][k] * b[k];
      if (add)
        a[i][j] += t;
      else
        a[i][j] = t;
    }
  }
}

/*
 * Argument a contains pointers to the rows of a symmetrix matrix.  The
 * in each row is the row number + 1.  These rows are stored in
 * contiguous memory starting with 0.  Evecs also contains pointers to
 * contiguous memory.  N is the dimension.
 */
void
cmat_diag(double**a, double*evals, double**evecs, int n,
              int matz, double tol)
{
  int i,j;
  int diagonal=1;
  double*fv1;

 /* I'm having problems with diagonalizing matrices which are already
  * diagonal.  So let's first check to see if _a_ is diagonal, and if it
  * is, then just return the diagonal elements in evals and a unit matrix
  * in evecs
  */

  for (i=1; i < n; i++) {
    for (j=0; j < i; j++) {
      if (fabs(a[i][j]) > tol) diagonal=0;
      }
    }

  if (diagonal) {
    for(i=0; i < n; i++) {
      evals[i] = a[i][i];
      evecs[i][i] = 1.0;

      for(j=0; j < i; j++) {
        evecs[i][j] = evecs[j][i] = 0.0;
        }
      }
    eigsort(n,evals,evecs);
    return;
    }

  fv1 = (double*) malloc(sizeof(double)*n);

  for(i=0; i < n; i++) {
      for(j=0; j <= i; j++) {
          evecs[i][j] = evecs[j][i] = a[i][j];
        }
    }

  tred2(n,evecs,evals,fv1,1);

  cmat_transpose_square_matrix(evecs,n);
  tqli(n,evals,evecs,fv1,1,tol);
  cmat_transpose_square_matrix(evecs,n);

  eigsort(n,evals,evecs);

  free(fv1);
  }

#ifndef __GNUC__
static double dsign(double a, double b)
#else
inline static double dsign(double a, double b)
#endif
{
  return (b >= 0.0) ? fabs(a) : -fabs(a);
}

static void
tred2(int n,double** a,double* d,double* e,int matz)
{
  int i,j,k,l;
  double f,g,h,hh,scale,scale_inv,h_inv;
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
        g= -(dsign(sqrt(h),f));
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

static void
tqli(int n, double* d, double** z, double* e, int matz, double toler)
{
  register int k;
  int i,l,m,iter;
  double g,r,s,c,p,f,b;
  double azi;

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
      if (iter++ == 30) {
        fprintf (stderr,"tqli not converging %d %g\n",l,e[l]);
        continue;
        }

      g = (d[l+1]-d[l])/(2.0*e[l]);
      r = sqrt(g*g + 1.0);
      g = d[m] - d[l] + e[l]/((g + dsign(r,g)));
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

static void
eigsort(int n, double* d, double** v)
{
  int i,j,k;
  double p;

  for(i=0; i < n-1 ; i++) {
    k=i;
    p=d[i];
    for(j=i+1; j < n; j++) {
      if(d[j] < p) {
        k=j;
        p=d[j];
        }
      }
    if(k != i) {
      d[k]=d[i];
      d[i]=p;
      for(j=0; j < n; j++) {
        p=v[j][i];
        v[j][i]=v[j][k];
        v[j][k]=p;
        }
      }
    }
  }
