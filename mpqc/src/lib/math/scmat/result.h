
#ifndef _math_scmat_result_h
#define _math_scmat_result_h
#ifdef __GNUC__
#pragma interface
#endif

#include <util/misc/compute.h>
#include <math/scmat/matrix.h>

typedef AccResult<RefSCMatrix> AccResultRefSCMatrix;
typedef AccResult<RefSymmSCMatrix> AccResultRefSymmSCMatrix;
typedef AccResult<RefDiagSCMatrix> AccResultRefDiagSCMatrix;
typedef AccResult<RefSCVector> AccResultRefSCVector;

typedef AccResult<RefSCMatrix> ResultRefSCMatrix;
typedef AccResult<RefSymmSCMatrix> ResultRefSymmSCMatrix;
typedef AccResult<RefDiagSCMatrix> ResultRefDiagSCMatrix;
typedef AccResult<RefSCVector> ResultRefSCVector;

#endif
