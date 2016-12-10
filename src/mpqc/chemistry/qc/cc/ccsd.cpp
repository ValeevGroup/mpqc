//
// Created by Chong Peng on 11/2/16.
//

#include "mpqc/chemistry/qc/cc/ccsd.h"
#include "mpqc/util/keyval/forcelink.h"

template class mpqc::lcao::CCSD<TA::TensorD, TA::SparsePolicy>;

MPQC_CLASS_EXPORT2("CCSD", mpqc::lcao::CCSD<TA::TensorD, TA::SparsePolicy>);
