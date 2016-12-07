//
// Created by Chong Peng on 11/15/16.
//

#include "ccsd_t_f12.h"
#include "mpqc/util/keyval/forcelink.h"

#if TA_DEFAULT_POLICY == 1
template class mpqc::f12::CCSD_T_F12<TA::TensorD>;
MPQC_CLASS_EXPORT2("CCSD(T)F12", mpqc::f12::CCSD_T_F12<TA::TensorD>);
#endif
