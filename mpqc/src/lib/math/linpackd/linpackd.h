#ifndef _math_linpackd_linpackd_h
#define _math_linpackd_linpackd_h

#ifdef __cplusplus
extern "C" {
#endif

int daxpy_(int *n, double *da, double *dx, int *incx, double *dy, int *incy);
int dcopy_(int *n, double *dx, int *incx, double *dy, int *incy);
int dscal_(int *n, double *da, double *dx, int *incx);
int dswap_(int *n, double *dx, int *incx, double *dy, int *incy);
int dsvdc_(double *x, int *ldx, int *n, int *p,
           double *s, double *e, double *u, int *ldu,
           double *v, int *ldv, double *work, int *job, int *info);
int drotg_(double *da, double *db, double *c, double *s);
int drot_(int *n, double *dx, int *incx, double *dy, int *incy,
          double *c, double *s);
double dnrm2_(int *n, double *dx, int *incx);
double ddot_(int *n, double *dx, int *incx, double *dy, int *incy);

#ifdef __cplusplus
}
#endif

#endif
