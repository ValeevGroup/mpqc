//
// Created by Chong Peng on 11/2/16.
//

#include "mpqc/chemistry/qc/cc/ccsd_t.h"
#include "mpqc/util/keyval/forcelink.h"

template class  mpqc::cc::CCSD_T<TA::TensorD, TA::SparsePolicy>;

MPQC_CLASS_EXPORT2("CCSD(T)", mpqc::cc::CCSD_T<TA::TensorD, TA::SparsePolicy>);
