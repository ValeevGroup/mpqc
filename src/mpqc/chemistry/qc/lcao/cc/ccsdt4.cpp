// CCSDT-4 by Varun Rishi, August 2017
//

#include "mpqc/chemistry/qc/lcao/cc/ccsdt4.h"
#include "mpqc/util/keyval/forcelink.h"

#if TA_DEFAULT_POLICY == 0
template class mpqc::lcao::CCSDT4<TA::TensorD, TA::DensePolicy>;
MPQC_CLASS_EXPORT2("CCSDT-4", mpqc::lcao::CCSDT4<TA::TensorD, TA::DensePolicy>);
#elif TA_DEFAULT_POLICY == 1
template class mpqc::lcao::CCSDT4<TA::TensorD, TA::SparsePolicy>;
MPQC_CLASS_EXPORT2("CCSDT-4", mpqc::lcao::CCSDT4<TA::TensorD, TA::SparsePolicy>);
#endif
