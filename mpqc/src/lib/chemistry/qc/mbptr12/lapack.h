extern "C" {
#include <chemistry/qc/mbptr12/f77sym.h>

extern void F77_DGESVD(const char* jobu, const char* jobvt, const int* m,
const int* n, double* A, const int* lda, double* S, double* U, const int* ldu,
double* Vt, const int* ldvt, double* work, int* lwork, int* info);

}

