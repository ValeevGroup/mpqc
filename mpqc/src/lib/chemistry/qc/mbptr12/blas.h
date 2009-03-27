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

namespace sc {

  //
  // taken from psi3/src/lib/libqt/blas_intfc.c
  //

  /*!
  ** C_DGEMM()
  **
  ** This function calculates C(m,n)=alpha*(opT)A(m,k)*(opT)B(k,n)+ beta*C(m,n)
  **
  ** These arguments mimic their Fortran conterparts; parameters have been
  ** reversed nca, ncb, ncc, A, B, C,  to make it correct for C.
  **
  ** \param char transa: On entry, specifies the form of (op)A used in the
  **                    matrix multiplication:
  **                    If transa = 'N' or 'n', (op)A = A.
  **                    If transa = 'T' or 't', (op)A = transp(A).
  **                    If transa = 'R' or 'r', (op)A = conjugate(A).
  **                    If transa = 'C' or 'c', (op)A = conjug_transp(A).
  **                    On exit, transa is unchanged.
  **
  ** \param char transb: On entry, specifies the form of (op)B used in the
  **                    matrix multiplication:
  **                    If transb = 'N' or 'n', (op)B = B.
  **                    If transb = 'T' or 't', (op)B = transp(B).
  **                    If transb = 'R' or 'r', (op)B = conjugate(B)
  **
  ** \param int m:      On entry, the number of rows of the matrix (op)A and of
  **                    the matrix C; m >= 0. On exit, m is unchanged.
  **
  ** \param int n:      On entry, the number of columns of the matrix (op)B and
  **                    of the matrix C; n >= 0. On exit, n is unchanged.
  **
  ** \param int k:      On entry, the number of columns of the matrix (op)A and
  **                    the number of rows of the matrix (op)B; k >= 0. On exit,
  **                    k is unchanged.
  **
  ** \param double alpha:  On entry, specifies the scalar alpha. On exit, alpha is
  **                       unchanged.
  **
  ** \param double* A:  On entry, a two-dimensional array A with dimensions ka
  **                    by nca. For (op)A = A  or  conjugate(A), nca >= k and the
  **                    leading m by k portion of the array A contains the matrix
  **                    A. For (op)A = transp(A) or conjug_transp(A), nca >= m
  **                    and the leading k by m part of the array A contains the
  **                    matrix A. On exit, a is unchanged.
  **
  ** \param int nca:    On entry, the second dimension of array A.
  **                    For (op)A = A  or conjugate(A), nca >= MAX(1,k).
  **                    For (op)A=transp(A) or conjug_transp(A), nca >= MAX(1,m).
  **                    On exit, nca is unchanged.
  **
  ** \param double* B:  On entry, a two-dimensional array B with dimensions kb
  **                  by ncb. For (op)B = B or conjugate(B), kb >= k and the
  **                  leading k by n portion of the array contains the matrix
  **                  B. For (op)B = transp(B) or conjug_transp(B), ncb >= k and
  **                  the leading n by k part of the array contains the matrix
  **                  B. On exit, B is unchanged.
  **
  ** \param int ncb:    On entry, the second dimension of array B.
  **                    For (op)B = B or <conjugate(B), ncb >= MAX(1,n).
  **                    For (op)B = transp(B) or conjug_transp(B), ncb >=
  **                    MAX(1,k). On exit, ncb is unchanged.
  **
  ** \param double beta: On entry, specifies the scalar beta. On exit, beta is
  **                    unchanged.
  **
  ** \param double* C:  On entry, a two-dimensional array with the dimension
  **                    at least m by ncc. On exit,  the leading  m by n part of
  **                    array C is overwritten by the matrix alpha*(op)A*(op)B +
  **                    beta*C.
  **
  ** \param int ncc:    On entry, the second dimension  of array C;
  **                    ncc >=MAX(1,n).  On exit, ncc is unchanged.
  */
  void C_DGEMM(char transa, char transb, int m, int n, int k, double alpha,
               const double *A, int nca, const double *B, int ncb, double beta,
               double *C, int ncc);

}
