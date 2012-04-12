
/* 
 * These routines are based on the work of Edward T. Seidl at the
 * National Institutes of Health.
 */

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string.h>
#include <math/scmat/cmatrix.h>
#include <math/scmat/lapack.h>
#include <util/misc/consumableresources.h>

using namespace sc;

extern "C" {

static void ludcmp(double**, int, int*, double*);
static void lubksb(double**, int, int*, double*);
static void symm_lu_decomp(double**, int, double*);
static void symm_lu_back_sub(double**, int, double*);

static void tred2(int dim,double**,double*,double*,int);
static void tqli(int dim,double*,double**,double*,int,double);
static void eigsort(int dim,double*,double**);

double**
cmat_new_square_matrix(int n)
{
  double *mat;
  double **r;
  if (n == 0) return 0;
  mat = allocate<double>(n*n);
  r = new double*[n];
  cmat_matrix_pointers(r,mat,n,n);
  return r;
}

double**
cmat_new_rect_matrix(int n,int m)
{
  double *mat;
  double **r;
  if (n == 0 || m == 0) return 0;
  mat = allocate<double>(n*m);
  r = new double*[n];
  cmat_matrix_pointers(r,mat,n,m);
  return r;
}

/* this deletes both square and triangular matrices */
void
cmat_delete_matrix(double**m)
{
  if (m) {
      deallocate(m[0]);
      delete[] m;
    }
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

  if (nr == 0 || nc == 0) return;

  if (nr == nc) {
      cmat_transpose_square_matrix(a,nr);
      return;
    };

  tmp = allocate<double>(nr*nc);
  tmpp = tmp;
  for (i=0; i<nc; i++) {
      for (j=0; j<nr; j++) {
          *tmpp = a[j][i];
          tmpp++;
        }
    }

  memcpy(a[0],tmp,sizeof(double)*nr*nc);

  deallocate(tmp);
}

/* a is symmetric if sym is true */
double
cmat_determ(double** a, int sym, int dim)
{
  int i;
  double det=0;

  if (sym) {
    symm_lu_decomp(a,dim,&det);
  } else {
    int *indx= new int[dim];
    ludcmp(a,dim,indx,&det);
    delete[] indx;
  }

  if (fabs(det) < 1.0e-16) return 0;

  for (i=0; i < dim; i++) det *= a[i][i];

  return det;
}

/* a is symmetric if sym is true */
double
cmat_solve_lin(double** a, int sym, double* b, int dim)
{
  int i;
  double det=0;

  if (sym) {
    symm_lu_decomp(a,dim,&det);
    if (fabs(det) < 1.0e-16) return 0;
    symm_lu_back_sub(a,dim,b);
  } else {
    int *indx= new int[dim];
    ludcmp(a,dim,indx,&det);
    if (fabs(det) < 1.0e-16) return 0;
    lubksb(a,dim,indx,b);
    delete[] indx;
  }

  for(i=0; i < dim; i++) det *= a[i][i];
  if (fabs(det) < 1.0e-16) return 0;

  return det;
}

double 
cmat_invert(double**a, int sym, int dim)
{
  int i,j;
  double det=0;
  double **y;
  double *b;

  b = new double[dim];
  y = cmat_new_square_matrix(dim);

  if (sym) {
    symm_lu_decomp(a,dim,&det);
    if (fabs(det) < 1.0e-16) return 0;

    for (i=0; i < dim; i++) det *= a[i][i];
    if (fabs(det) < 1.0e-16) return 0;

    for (i=0; i < dim; i++) {
      for (j=0; j < dim; j++) b[j]=0;
      b[i]=1;
      symm_lu_back_sub(a,dim,b);
      for (j=0; j < dim; j++) y[j][i]=b[j];
    }

    for (i=0; i < dim; i++)
      for (j=0; j <= i; j++)
        a[i][j] = y[i][j];

  } else {
    int *indx= new int[dim];

    ludcmp(a,dim,indx,&det);
    if (fabs(det) < 1.0e-16) return 0;

    for (i=0; i < dim; i++) det *= a[i][i];
    if (fabs(det) < 1.0e-16) return 0;

    for (i=0; i < dim; i++) {
      memset(b,0,sizeof(double)*dim);
      b[i]=1;
      lubksb(a,dim,indx,b);
      for (j=0; j < dim; j++) y[j][i]=b[j];
    }

    for (i=0; i < dim; i++)
      for (j=0; j < dim; j++)
        a[i][j] = y[i][j];
    delete[] indx;
  }

  delete[] b;
  cmat_delete_matrix(y);

  return det;
}

static void
ludcmp(double** a, int n, int *indx, double *d)
{
  int i,imax=0,j,k;
  double big,dum,sum,temp;

  double* vv = new double[n];

  *d = 1.0;

  for (i=0; i < n ; i++) {
    big=0.0;
    for (j=0; j < n; j++) if ((temp=fabs(a[i][j])) > big) big=temp;
#if 1
    if (big == 0.0) {
      *d = 0.0;
      delete[] vv;
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
  delete[] vv;
  }

static void
lubksb(double** a, int n, int *indx, double* b)
{
  int i,ii=0,ip,j;
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
 * this is LU decomposition where A is a symmetric matrix
 * when A is symmetric, then 
 *   beta(i,j) = A(i,j) - sum_k(i-1) beta(k,i)*beta(k,j)/beta(k,k)
 *   alpha(i,j) = beta(j,i)/beta(j,j)
 *
 * since we're storing beta in a, the indices of beta will be switched
 * since alpha is expressed in terms of beta, we don't store it
 *
 * so we have
 *   beta(i,j) = A(i,j) - sum_k(i-1) beta(i,k)*beta(j,k)/beta(k,k)
 *   alpha(i,j) = beta(i,j)/beta(j,j)
 */

static void
symm_lu_decomp(double** a, int n, double *d)
{
  int i,j,k;
  double tmp;

  std::vector<double> v(n, 0.0);

  /* check for singular matrix */
  for (i=0; i < n; i++) {
    for (j=0; j < i; j++) {
      v[i] = ((tmp=fabs(a[i][j])) > v[i]) ? tmp : v[i];
      v[j] = (tmp > v[j]) ? tmp : v[j];
    }
    v[i] = ((tmp=fabs(a[i][i])) > v[i]) ? tmp : v[i];
  }

  for (i=0; i < n; i++) {
    if (fabs(v[i]) < 1.0e-16) {
      fprintf(stderr,"\n  warning: singular matrix in symm_lu_decomp\n");
      *d = 0.0;
      return;
    }
  }

  *d = 1.0;

  for (i=0; i < n ; i++) {
    /* check to make sure we're not going to blow up */
    if (i < n-1) {
      tmp = 0; 
      for (k=0; k < i-1; k++)
        tmp += a[i][k]*a[i][k]/a[k][k];
      if (fabs(a[i][i]-tmp) < 1.0e-16) {
        fprintf(stderr,"\n  warning: singular matrix in symm_lu_decomp 2\n");
        *d = 0;
        return;
      }
    }
    for (j=i; j < n; j++) {
      tmp = 0;
      for (k=0; k <= i-1; k++)
        tmp -= a[i][k]*a[j][k]/a[k][k];
      a[j][i] += tmp;
    }
  }
}

static void
symm_lu_back_sub(double** a, int n, double* b)
{
  int i,j;
  double sum;

 /* form y(i) = bi - sum_j(i-1) alpha(i,j)*y(j)
  * alpha(i,j) = beta(j,i)/beta(j,j), but beta is stored lower instead of 
  * upper triangle, so alpha(i,j) = beta(i,j)/beta(j,j)
  */
  for (i=0; i < n ; i++) {
    sum = 0;
    for (j=0; j < i; j++)
      sum += (a[i][j]/a[j][j]) * b[j];
    b[i] -= sum;
  }

 /* now form x(i) = 1/beta(i,i)*[y(i) - sum_j=i+1(N) beta(i,j)*x(j)]
  * is really ...[...beta(j,i)*x(j)]
  */
  for (i=n-1; i >= 0 ; i--) {
    sum = b[i];
    for (j=i+1; j < n ; j++) sum -= a[j][i]*b[j];
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
          a = new double*[nr];
          a[0] = old_a[0];
        }
      cmat_matrix_pointers(a,a[0],nr,nl);
    }
  if(!tb) {
      cmat_transpose_matrix(b,nl,nc);
      if (nc > nl) {
          old_b = b;
          b = new double*[nc];
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
          delete[] a;
          a = old_a;
        }
      cmat_matrix_pointers(a,a[0],nr,nl);
    }
  if(!tb) {
      cmat_transpose_matrix(b,nc,nl);
      if (old_b) {
          delete[] b;
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
          if (add) tmp = ai[j];
          else tmp = 0.0;
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
  brow = new double[nb];
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
  delete[] brow;

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

  fv1 = new double[n];

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

  delete[] fv1;
  }

#define dsign(a,b) (((b) >= 0.0) ? fabs(a) : -fabs(a))

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

void
cmat_schmidt(double **C, double *S, int nrow, int nc)
{
  int i,j,ij;
  int m;
  double vtmp;
  std::vector<double> v(nrow);

  for (m=0; m < nc; m++) {
    v[0] = C[0][m] * S[0];

    for (i=ij=1; i < nrow; i++) {
      for (j=0,vtmp=0.0; j < i; j++,ij++) {
        vtmp += C[j][m]*S[ij];
        v[j] += C[i][m]*S[ij];
      }
      v[i] = vtmp + C[i][m]*S[ij];
      ij++;
    }

    for (i=0,vtmp=0.0; i < nrow; i++)
      vtmp += v[i]*C[i][m];

    if (!vtmp) {
      fprintf(stderr,"cmat_schmidt: bogus\n");
      abort();
    }

    if (vtmp < 1.0e-15)
      vtmp = 1.0e-15;

    vtmp = 1.0/sqrt(vtmp);
    
    for (i=0; i < nrow; i++) {
      v[i] *= vtmp;
      C[i][m] *= vtmp;
    }

    if (m < nc-1) {
      for (i=m+1,vtmp=0.0; i < nc; i++) {
        for (j=0,vtmp=0.0; j < nrow; j++)
          vtmp += v[j] * C[j][i];
        for (j=0; j < nrow; j++)
          C[j][i] -= vtmp * C[j][m];
      }
    }
  }
}

/* Returns the number of linearly independent vectors
   orthogonal wrt S. */
int
cmat_schmidt_tol(double **C, double *S, int nrow, int ncol,
                 double tolerance, double *res)
{
  int i,j,ij;
  int m;
  double vtmp;
  int northog = 0;
  std::vector<double> v(nrow);

  if (res) *res = 1.0;

  /* Orthonormalize the columns of C wrt S. */
  for (m=0; m < ncol; m++) {
    v[0] = C[0][m] * S[0];

    for (i=ij=1; i < nrow; i++) {
      for (j=0,vtmp=0.0; j < i; j++,ij++) {
        vtmp += C[j][m]*S[ij];
        v[j] += C[i][m]*S[ij];
      }
      v[i] = vtmp + C[i][m]*S[ij];
      ij++;
    }

    for (i=0,vtmp=0.0; i < nrow; i++)
      vtmp += v[i]*C[i][m];

    if (vtmp < tolerance) continue;

    if (res && (m == 0 || vtmp < *res)) *res = vtmp;

    vtmp = 1.0/sqrt(vtmp);
    
    for (i=0; i < nrow; i++) {
      v[i] *= vtmp;
      C[i][northog] = C[i][m] * vtmp;
    }

    for (i=m+1,vtmp=0.0; i < ncol; i++) {
        for (j=0,vtmp=0.0; j < nrow; j++)
            vtmp += v[j] * C[j][i];
        for (j=0; j < nrow; j++)
            C[j][i] -= vtmp * C[j][northog];
      }
    northog++;
  }
  return northog;
}

void
cmat_eigensystem(/*const*/ double**atri, /*const*/ double**stri, double*evals, double**evecs, int n,
                 int matz)
{
  if (n == 0) return;

  // Convert to the lapack storage format.
  const int n2 = n*n;
  double *a = new double[n2];
  double *s = new double[n2];
  for (int i=0; i<n; i++) {
    for (int j=0; j<=i; j++) {
      a[i*n+j] = atri[i][j];
      a[j*n+i] = atri[i][j];
      s[i*n+j] = stri[i][j];
      s[j*n+i] = stri[i][j];
    }
  }

  // solve generalized eigenvalue problem with DSYGV
  const blasint itype = 1;
  const char jobz_V = 'V';
  const char uplo_U = 'U';  // lower triangle in C -> upper triange in Fortran
  blasint lwork = -1;
  blasint info;
  double optlwork;
  const blasint nn = n;
  F77_DSYGV(&itype,&jobz_V,&uplo_U, &nn,
            a, &nn,
            s, &nn,
            evals,
            &optlwork,&lwork,
            &info);
  if (info) {
    ExEnv::outn() << "dsygv could not determine work size: info = "
                  << info << std::endl;
    abort();
  }
  lwork = (blasint)optlwork;
  double *work = new double[lwork];
  F77_DSYGV(&itype,&jobz_V,&uplo_U, &nn,
              a, &nn,
              s, &nn,
              evals,
              work,&lwork,
              &info);
  if (info) {
    std::ostringstream oss;
    oss << "dsygv could not solve eigensystem: info = " << info;
    throw AlgorithmException(oss.str().c_str(),
                             __FILE__, __LINE__);
  }

  // the vector is placed in a -> transpose to C order
  int ij=0;
  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++, ++ij) {
      evecs[j][i] = a[ij];
    }
  }

  // cleanup
  delete[] a;
  delete[] s;
  delete[] work;
}

} // end of extern "C"
