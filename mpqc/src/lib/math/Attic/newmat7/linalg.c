/* 
 * Blas 1 routines written in C
 * Implemented: ddot, dnrm2, dscal
 * Blas 1: dasum, daxpy, dcopy, ddot, dnrm2, drot, drotg, dscal, dswap
 *
 * J.C. Meza
 * Sandia National Laboratories
 * meza@california.sandia.gov
 */

#include "cblas.h"
double ddot(int n, double *dx, int incx, double *dy, int incy)

{
  /* Local variables */
  int i, m;
  double dtemp = 0.;
  int ix, iy;
  
  /* Function Body */
  /*     forms the dot product of two vectors. */
  /*     uses unrolled loops for increments equal to one. */
  /*     jack dongarra, linpack, 3/11/78. */
  

  dtemp = 0.;
  if (n <= 0) return 0.;

  if (incx == 1 && incy == 1) {/* code for both increments equal to 1 */
    
    /* clean-up loop */
    
    m = n % 5;
    if (m != 0) {
      for (i=0; i<m; ++i) dtemp += dx[i]*dy[i];
      if (n < 5) return dtemp;
    }
    
    for (i=m; i<n; i+=5) {
      dtemp = dtemp + dx[i]*dy[i] + dx[i+1]*dy[i+1] + dx[i+2]*dy[i+2] + 
	dx[i+3]*dy[i+3] + dx[i+4]*dy[i+4];
    }
  }
  
  else { /* code for unequal increments or equal increments not equal to 1 */
  
    ix = 1;    iy = 1;
    if (incx < 0) ix = (-n + 1) * incx + 1;
    if (incy < 0) iy = (-n + 1) * incy + 1;

    for (i=0; i<n; ++i) {
      dtemp += dx[ix] * dy[iy];
      ix += incx;
      iy += incy;
    }
  }
  return dtemp;
  
}

void dscal(int n, double da, double *dx, int incx)

{

  /* Local variables */
  static int i, m, ix;

  /* Function Body */
  /*     scales a vector by a constant. */
  /*     uses unrolled loops for increment equal to one. */
  /*     jack dongarra, linpack, 3/11/78. */
  /*     modified to correct problem with negative increment, 8/21/90. */
  

  if (n <= 0) return;
  if (incx == 1) {/* code for increment equal to 1 */

    /* clean-up loop */
    m = n % 5;
    if (m != 0) {
      for (i = 0; i < m; ++i) dx[i] = da * dx[i];
      if (n < 5) return;
    }

    for (i = m; i < n; i += 5) {
      dx[i] = da * dx[i];
      dx[i + 1] = da * dx[i + 1];
      dx[i + 2] = da * dx[i + 2];
      dx[i + 3] = da * dx[i + 3];
      dx[i + 4] = da * dx[i + 4];
    }
  }
  
  else {/* code for increment not equal to 1 */
      
    ix = 1;
    if (incx < 0) ix = (-n + + 1) * incx + 1;
    for (i = 0; i < n; ++i) {
      dx[ix] = da * dx[ix];
      ix += incx;
    }
  }
}

double dnrm2(int n, double *dx, int incx)
{
  int i, ix;
  double sum;

  sum = 0.0;
  if (incx == 1) {
    for (i=0; i<n; ++i) {
      sum += dx[i]*dx[i];
    }
  }
  else {
    ix = 1;
    for (i=0; i<n; ++i) {
      sum += dx[ix];
      ix = ix + incx;
    }
  }
  
  return sqrt(sum);

}

