// CCSDT by Varun Rishi, September 2017
//

#include "mpqc/chemistry/qc/lcao/cc/ccsdt.h"
#include "mpqc/util/keyval/forcelink.h"

#if TA_DEFAULT_POLICY == 0
template class mpqc::lcao::CCSDT<TA::TensorD, TA::DensePolicy>;
MPQC_CLASS_EXPORT2("CCSDT", mpqc::lcao::CCSDT<TA::TensorD, TA::DensePolicy>);
#elif TA_DEFAULT_POLICY == 1
template class mpqc::lcao::CCSDT<TA::TensorD, TA::SparsePolicy>;
MPQC_CLASS_EXPORT2("CCSDT", mpqc::lcao::CCSDT<TA::TensorD, TA::SparsePolicy>);
#endif
