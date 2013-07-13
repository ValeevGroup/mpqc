
#include <math/scmat/blas.h>

extern "C" {
#include <math/scmat/f77sym.h>

extern void F77_DGESVD(const char* jobu, const char* jobvt, const blasint* m,
                       const blasint* n, double* A, const blasint* lda, double* S, double* U, const blasint* ldu,
                       double* Vt, const blasint* ldvt, double* work, blasint* lwork, blasint* info);

extern void F77_DSPSVX(const char* fact, const char* uplo, const blasint* n, const blasint* nrhs,
                       const double* AP, double* AFP, blasint* ipiv, const double* BB, const blasint* nb,
                       double* XX, const blasint* nx, double* rcond, double* ferr, double* berr,
                       double* work, blasint* iwork, blasint* info);

extern void F77_DSYEVD(const char* need_evals, const char* uplo, const blasint* n,
                       double* Asq, const blasint* lda, double* evals, double* work, const blasint* lwork,
                       blasint* iwork, const blasint* liwork, blasint* info);

extern void F77_DSPTRF(const char* uplo, const blasint* n, double* AP, blasint* ipiv, blasint* info);

extern void F77_DPPTRF(const char* uplo, const blasint* n, double* AP, blasint* info);

extern void F77_DSPTRI(const char* uplo, const blasint* n, double* AP, const blasint* ipiv, double* work, blasint* info);

extern void F77_DPPTRI(const char* uplo, const blasint* n, double* AP, blasint* info);

extern double F77_DLANSP(const char* norm, const char* uplo, const blasint* n, const double* A_packed, double* work);

extern void F77_DSPCON(const char* uplo, const blasint* n, const double* A_packed, const blasint* ipiv,
                       const double* anorm, double* rcond, double* work, blasint* iwork, blasint* info);

extern void F77_DPPCON(const char* uplo, const blasint* n, const double* A_packed,
                       const double* anorm, double* rcond, double* work, blasint* iwork, blasint* info);

extern double F77_DLAMCH(const char* e);

extern void F77_DLACPY(const char* uplo, const blasint* m, const blasint* n, const double* A, const blasint* lda,
                       double* B, const blasint* ldb, blasint* info);

extern void F77_DSPTRS(const char* uplo, const blasint* n, const blasint* nrhs, const double* AFP, const blasint* ipiv,
                       const double* X, const blasint* ldx, blasint* info);

extern void F77_DPPTRS(const char* uplo, const blasint* n, const blasint* nrhs, const double* AFP,
                       const double* X, const blasint* ldx, blasint* info);

extern void F77_DSPRFS(const char* uplo, const blasint* n, const blasint* nrhs, const double* A, const double* AF,
                       const blasint* ipiv, const double* B, const blasint* ldb, const double* X,
                       const blasint* ldx, double* ferr, double* berr, double* work, blasint* iwork, blasint* info);

extern void F77_DPPRFS(const char* uplo, const blasint* n, const blasint* nrhs, const double* A, const double* AF,
                       const double* B, const blasint* ldb, const double* X,
                       const blasint* ldx, double* ferr, double* berr, double* work, blasint* iwork, blasint* info);

extern void F77_DSYGV(const blasint* itype, const char* jobz, const char* uplo, const blasint* n,
                      double* Asq, const blasint* lda, double* Bsq, const blasint* ldb, double* evals,
                      double* work, const blasint* lwork,
                      blasint* info);

}

