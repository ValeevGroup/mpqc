#include "mpqc/chemistry/qc/lcao/cc/gamma_point_ccsd.h"
#include "mpqc/util/keyval/forcelink.h"

#if TA_DEFAULT_POLICY == 0
template class mpqc::lcao::GammaPointCCSD<TA::TensorZ, TA::DensePolicy>;
MPQC_CLASS_EXPORT2("GammaPointCCSD", mpqc::lcao::GammaPointCCSD<TA::TensorZ, TA::DensePolicy>);
#elif TA_DEFAULT_POLICY == 1
template class mpqc::lcao::GammaPointCCSD<TA::TensorZ, TA::SparsePolicy>;
MPQC_CLASS_EXPORT2("GammaPointCCSD", mpqc::lcao::GammaPointCCSD<TA::TensorZ, TA::SparsePolicy>);
#endif
