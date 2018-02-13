// CCSDT-1b by Varun Rishi, August 2017
//

#include "mpqc/chemistry/qc/lcao/cc/ccsdt1b.h"
#include "mpqc/util/keyval/forcelink.h"

#if TA_DEFAULT_POLICY == 0
template class mpqc::lcao::CCSDT1b<TA::TensorD, TA::DensePolicy>;
MPQC_CLASS_EXPORT2("CCSDT-1b", mpqc::lcao::CCSDT1b<TA::TensorD, TA::DensePolicy>);
#elif TA_DEFAULT_POLICY == 1
template class mpqc::lcao::CCSDT1b<TA::TensorD, TA::SparsePolicy>;
MPQC_CLASS_EXPORT2("CCSDT-1b", mpqc::lcao::CCSDT1b<TA::TensorD, TA::SparsePolicy>);
#endif
