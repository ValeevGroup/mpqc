//
// Created by Chong Peng on 10/13/16.
//

#include "mp2f12.h"
#include "mpqc/util/keyval/forcelink.h"

template class mpqc::f12::RMP2F12<TA::TensorD>;
template class mpqc::f12::RIRMP2F12<TA::TensorD>;

MPQC_CLASS_EXPORT2("RMP2F12", mpqc::f12::RMP2F12<TA::TensorD>);
MPQC_CLASS_EXPORT2("RI-RMP2F12", mpqc::f12::RIRMP2F12<TA::TensorD>);


