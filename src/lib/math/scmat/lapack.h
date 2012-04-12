
#include <math/scmat/blas.h>

extern "C" {
#include <math/scmat/f77sym.h>

extern void F77_DGESVD(const char* jobu, const char* jobvt, const blas_f77_integer_t* m,
                       const blas_f77_integer_t* n, double* A, const blas_f77_integer_t* lda, double* S, double* U, const blas_f77_integer_t* ldu,
                       double* Vt, const blas_f77_integer_t* ldvt, double* work, blas_f77_integer_t* lwork, blas_f77_integer_t* info);

extern void F77_DSPSVX(const char* fact, const char* uplo, const blas_f77_integer_t* n, const blas_f77_integer_t* nrhs,
                       const double* AP, double* AFP, blas_f77_integer_t* ipiv, const double* BB, const blas_f77_integer_t* nb,
                       double* XX, const blas_f77_integer_t* nx, double* rcond, double* ferr, double* berr,
                       double* work, blas_f77_integer_t* iwork, blas_f77_integer_t* info);

extern void F77_DSYEVD(const char* need_evals, const char* uplo, const blas_f77_integer_t* n,
                       double* Asq, const blas_f77_integer_t* lda, double* evals, double* work, const blas_f77_integer_t* lwork,
                       blas_f77_integer_t* iwork, const blas_f77_integer_t* liwork, blas_f77_integer_t* info);

extern void F77_DSPTRF(const char* uplo, const blas_f77_integer_t* n, double* AP, blas_f77_integer_t* ipiv, blas_f77_integer_t* info);

extern void F77_DPPTRF(const char* uplo, const blas_f77_integer_t* n, double* AP, blas_f77_integer_t* info);

extern void F77_DSPTRI(const char* uplo, const blas_f77_integer_t* n, double* AP, const blas_f77_integer_t* ipiv, double* work, blas_f77_integer_t* info);

extern void F77_DPPTRI(const char* uplo, const blas_f77_integer_t* n, double* AP, blas_f77_integer_t* info);

extern double F77_DLANSP(const char* norm, const char* uplo, const blas_f77_integer_t* n, const double* A_packed, double* work);

extern void F77_DSPCON(const char* uplo, const blas_f77_integer_t* n, const double* A_packed, const blas_f77_integer_t* ipiv,
                       const double* anorm, double* rcond, double* work, blas_f77_integer_t* iwork, blas_f77_integer_t* info);

extern void F77_DPPCON(const char* uplo, const blas_f77_integer_t* n, const double* A_packed,
                       const double* anorm, double* rcond, double* work, blas_f77_integer_t* iwork, blas_f77_integer_t* info);

extern double F77_DLAMCH(const char* e);

extern void F77_DLACPY(const char* uplo, const blas_f77_integer_t* m, const blas_f77_integer_t* n, const double* A, const blas_f77_integer_t* lda,
                       double* B, const blas_f77_integer_t* ldb, blas_f77_integer_t* info);

extern void F77_DSPTRS(const char* uplo, const blas_f77_integer_t* n, const blas_f77_integer_t* nrhs, const double* AFP, const blas_f77_integer_t* ipiv,
                       const double* X, const blas_f77_integer_t* ldx, blas_f77_integer_t* info);

extern void F77_DPPTRS(const char* uplo, const blas_f77_integer_t* n, const blas_f77_integer_t* nrhs, const double* AFP,
                       const double* X, const blas_f77_integer_t* ldx, blas_f77_integer_t* info);

extern void F77_DSPRFS(const char* uplo, const blas_f77_integer_t* n, const blas_f77_integer_t* nrhs, const double* A, const double* AF,
                       const blas_f77_integer_t* ipiv, const double* B, const blas_f77_integer_t* ldb, const double* X,
                       const blas_f77_integer_t* ldx, double* ferr, double* berr, double* work, blas_f77_integer_t* iwork, blas_f77_integer_t* info);

extern void F77_DPPRFS(const char* uplo, const blas_f77_integer_t* n, const blas_f77_integer_t* nrhs, const double* A, const double* AF,
                       const double* B, const blas_f77_integer_t* ldb, const double* X,
                       const blas_f77_integer_t* ldx, double* ferr, double* berr, double* work, blas_f77_integer_t* iwork, blas_f77_integer_t* info);

extern void F77_DSYGV(const blas_f77_integer_t* itype, const char* jobz, const char* uplo, const blas_f77_integer_t* n,
                      double* Asq, const blas_f77_integer_t* lda, double* Bsq, const blas_f77_integer_t* ldb, double* evals,
                      double* work, const blas_f77_integer_t* lwork,
                      blas_f77_integer_t* info);

}

