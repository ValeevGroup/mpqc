//
// Created by Chong Peng on 11/1/16.
//

#include "mpqc/chemistry/qc/lcao/f12/gf2f12.h"
#include "mpqc/util/keyval/forcelink.h"

#if TA_DEFAULT_POLICY == 1
template class mpqc::lcao::GF2F12<TA::TensorD>;
MPQC_CLASS_EXPORT2("GF2F12", mpqc::lcao::GF2F12<TA::TensorD>);
#endif