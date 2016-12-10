//
// Created by Chong Peng on 11/2/16.
//

#include "mpqc/chemistry/qc/cc/dbccsd.h"
#include "mpqc/util/keyval/forcelink.h"

template class mpqc::lcao::DBCCSD<TA::TensorD, TA::SparsePolicy>;

MPQC_CLASS_EXPORT2("DBCCSD", mpqc::lcao::DBCCSD<TA::TensorD, TA::SparsePolicy>);
