//
// Created by Chong Peng on 11/1/16.
//

#include "ccsd_f12.h"
#include "mpqc/util/keyval/forcelink.h"

#if TA_DEFAULT_POLICY == 1
template class mpqc::lcao::CCSD_F12<TA::TensorD>;
MPQC_CLASS_EXPORT2("CCSD(F12)", mpqc::lcao::CCSD_F12<TA::TensorD>);
#endif
