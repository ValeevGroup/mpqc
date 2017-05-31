#include "eom_ccsd.h"
#include "mpqc/util/keyval/forcelink.h"


#if TA_DEFAULT_POLICY == 0
template class mpqc::lcao::EOM_CCSD<TA::TensorD, TA::DensePolicy>;
MPQC_CLASS_EXPORT2("EOM-CCSD", mpqc::lcao::EOM_CCSD<TA::TensorD, TA::DensePolicy>);
#elif TA_DEFAULT_POLICY == 1
template class mpqc::lcao::EOM_CCSD<TA::TensorD, TA::SparsePolicy>;
MPQC_CLASS_EXPORT2("EOM-CCSD", mpqc::lcao::EOM_CCSD<TA::TensorD, TA::SparsePolicy>);
#endif
