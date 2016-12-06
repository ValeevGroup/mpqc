//
// Created by Chong Peng on 10/6/16.
//

#include "mp2.h"
#include "mpqc/util/keyval/forcelink.h"

template class mpqc::mbpt::RMP2<TA::TensorD, TA::SparsePolicy>;
template class mpqc::mbpt::RIRMP2<TA::TensorD, TA::SparsePolicy>;

MPQC_CLASS_EXPORT2("RMP2", mpqc::mbpt::RMP2<TA::TensorD, TA::SparsePolicy>);
MPQC_CLASS_EXPORT2("RI-RMP2", mpqc::mbpt::RIRMP2<TA::TensorD, TA::SparsePolicy>);


