//
// Created by Chong Peng on 3/1/17.
//

#include "mpqc/chemistry/qc/lcao/ci/cis.h"
#include "mpqc/util/keyval/forcelink.h"

#if TA_DEFAULT_POLICY == 0
template class mpqc::lcao::CIS<TA::TensorD, TA::DensePolicy>;
MPQC_CLASS_EXPORT2("CIS", mpqc::lcao::CIS<TA::TensorD, TA::DensePolicy>);
#elif TA_DEFAULT_POLICY == 1
template class mpqc::lcao::CIS<TA::TensorD, TA::SparsePolicy>;
MPQC_CLASS_EXPORT2("CIS", mpqc::lcao::CIS<TA::TensorD, TA::SparsePolicy>);
#endif