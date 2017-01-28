//
// Created by Chong Peng on 10/13/16.
//

#include "mp2f12.h"
#include "mpqc/util/keyval/forcelink.h"

#if TA_DEFAULT_POLICY == 1
template class mpqc::lcao::RMP2F12<TA::TensorD>;
template class mpqc::lcao::RIRMP2F12<TA::TensorD>;

MPQC_CLASS_EXPORT2("RMP2F12", mpqc::lcao::RMP2F12<TA::TensorD>);
MPQC_CLASS_EXPORT2("RI-RMP2F12", mpqc::lcao::RIRMP2F12<TA::TensorD>);
#endif