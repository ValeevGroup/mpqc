extern "C" {
#include <chemistry/qc/mbptr12/f77sym.h>

extern void F77_DGESVD(const char* jobu, const char* jobvt, const int* m,
const int* n, double* A, const int* lda, double* S, double* U, const int* ldu,
double* Vt, const int* ldvt, double* work, int* lwork, int* info);

extern void F77_DSPSVX(const char* fact, const char* uplo, const int* n, const int* nrhs,
                       const double* AP, double* AFP, int* ipiv, const double* BB, const int* nb,
                       double* XX, const int* nx, double* rcond, double* ferr, double* berr,
                       double* work, int* iwork, int* info);

}

