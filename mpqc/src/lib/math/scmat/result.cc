
#ifdef __GNUC__
#pragma implementation
#endif

#include <math/scmat/result.h>
#include <util/state/state.h>

#ifdef __GNUC__
template class AccResult<RefSCMatrix>;
template class AccResult<RefSymmSCMatrix>;
template class AccResult<RefDiagSCMatrix>;
template class AccResult<RefSCVector>;
#endif
