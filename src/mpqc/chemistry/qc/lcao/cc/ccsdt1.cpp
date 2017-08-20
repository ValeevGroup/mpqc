// Adapted to CCSDT-1 _VR_July18_17
// Created by Chong Peng on 11/2/16.
//

#include "mpqc/chemistry/qc/lcao/cc/ccsdt1.h"
#include "mpqc/util/keyval/forcelink.h"

#if TA_DEFAULT_POLICY == 0
template class mpqc::lcao::CCSDT1<TA::TensorD, TA::DensePolicy>;
MPQC_CLASS_EXPORT2("CCSDT-1", mpqc::lcao::CCSDT1<TA::TensorD, TA::DensePolicy>);
#elif TA_DEFAULT_POLICY == 1
template class mpqc::lcao::CCSDT1<TA::TensorD, TA::SparsePolicy>;
MPQC_CLASS_EXPORT2("CCSDT-1", mpqc::lcao::CCSDT1<TA::TensorD, TA::SparsePolicy>);
#endif
