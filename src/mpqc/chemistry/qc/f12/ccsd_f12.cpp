//
// Created by Chong Peng on 11/1/16.
//

#include "ccsd_f12.h"
#include "mpqc/util/keyval/forcelink.h"

template class mpqc::f12::CCSD_F12<TA::TensorD>;

MPQC_CLASS_EXPORT2("CCSD(F12)", mpqc::f12::CCSD_F12<TA::TensorD>);
