#include "mpqc/chemistry/qc/lcao/cc/gamma_point_ccsd.h"
#include "mpqc/util/keyval/forcelink.h"

#if TA_DEFAULT_POLICY == 0

#elif TA_DEFAULT_POLICY == 1
template class mpqc::lcao::GammaPointCCSDVersion2<TA::TensorD, TA::SparsePolicy>;
MPQC_CLASS_EXPORT2("GammaPointCCSDVersion2", mpqc::lcao::GammaPointCCSDVersion2<TA::TensorD, TA::SparsePolicy>);
#endif
