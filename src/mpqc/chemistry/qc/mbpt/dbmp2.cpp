//
// Created by Chong Peng on 11/2/16.
//

#include "mpqc/chemistry/qc/mbpt/dbmp2.h"
#include "mpqc/util/keyval/forcelink.h"

template class mpqc::mbpt::DBRMP2<TA::TensorD, TA::SparsePolicy>;
template class mpqc::mbpt::RIDBRMP2<TA::TensorD, TA::SparsePolicy>;

MPQC_CLASS_EXPORT2("DBRMP2", mpqc::mbpt::DBRMP2<TA::TensorD, TA::SparsePolicy>);
MPQC_CLASS_EXPORT2("RI-DBRMP2", mpqc::mbpt::RIDBRMP2<TA::TensorD, TA::SparsePolicy>);
