
#ifndef _math_scmat_result_h
#define _math_scmat_result_h
#ifdef __GNUC__
#pragma interface
#endif

#include <util/misc/compute.h>
#include <math/scmat/matrix.h>

SSAccResult_dec(RefSCMatrix);
SSAccResult_dec(RefSymmSCMatrix);
SSAccResult_dec(RefDiagSCMatrix);
SSAccResult_dec(RefSCVector);

typedef AccResultRefSCMatrix ResultRefSCMatrix;
typedef AccResultRefSymmSCMatrix ResultRefSymmSCMatrix;
typedef AccResultRefDiagSCMatrix ResultRefDiagSCMatrix;
typedef AccResultRefSCVector ResultRefSCVector;

#endif
