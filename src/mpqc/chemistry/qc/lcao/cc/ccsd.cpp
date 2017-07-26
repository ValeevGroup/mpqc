//
// Created by Chong Peng on 11/2/16.
//

#include "mpqc/chemistry/qc/lcao/cc/ccsd.h"
#include "mpqc/util/keyval/forcelink.h"

#if TA_DEFAULT_POLICY == 0
template class mpqc::lcao::CCSD<TA::TensorD, TA::DensePolicy>;
MPQC_CLASS_EXPORT2("CCSD", mpqc::lcao::CCSD<TA::TensorD, TA::DensePolicy>);
#elif TA_DEFAULT_POLICY == 1
template class mpqc::lcao::CCSD<TA::TensorD, TA::SparsePolicy>;
MPQC_CLASS_EXPORT2("CCSD", mpqc::lcao::CCSD<TA::TensorD, TA::SparsePolicy>);
#endif
