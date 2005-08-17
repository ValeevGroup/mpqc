extern "C" {
#include <chemistry/qc/mbptr12/f77sym.h>

extern void F77_DGEMM(const char*, const char*, const int*,
const int*, const int*, const double*, const double*, const int*,
const double*, const int*, const double*, double*, const int*);

extern void F77_DAXPY(const int* n, const double* da, const double* dx,
const int* incx, double* dy, const int* incy);

extern double F77_DDOT(const int* n, const double* dx, const int* incx,
const double* dy, const int* incy);

}

