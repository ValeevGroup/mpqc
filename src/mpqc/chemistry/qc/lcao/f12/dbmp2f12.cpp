//
// Created by Chong Peng on 10/10/16.
//

#include "dbmp2f12.h"
#include "mpqc/util/keyval/forcelink.h"

#if TA_DEFAULT_POLICY == 1
template class mpqc::lcao::RIDBRMP2F12<TA::TensorD>;
MPQC_CLASS_EXPORT2("RI-DBRMP2F12", mpqc::lcao::RIDBRMP2F12<TA::TensorD>);
#endif


