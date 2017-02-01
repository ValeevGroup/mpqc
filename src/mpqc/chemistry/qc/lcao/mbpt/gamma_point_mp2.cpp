#include "mpqc/chemistry/qc/lcao/mbpt/gamma_point_mp2.h"
#include "mpqc/util/keyval/forcelink.h"

#if TA_DEFAULT_POLICY == 0
template class mpqc::lcao::GammaPointMP2<TA::TensorZ, TA::DensePolicy>;
MPQC_CLASS_EXPORT2("GammaPointMP2", mpqc::lcao::GammaPointMP2<TA::TensorZ, TA::DensePolicy>);

#elif TA_DEFAULT_POLICY == 1
template class mpqc::lcao::GammaPointMP2<TA::TensorZ, TA::SparsePolicy>;
MPQC_CLASS_EXPORT2("GammaPointMP2", mpqc::lcao::GammaPointMP2<TA::TensorZ, TA::SparsePolicy>);

#endif

