//
// Created by Chong Peng on 11/1/16.
//

#include "mpqc/chemistry/qc/lcao/f12/dbccsd_f12.h"
#include "mpqc/util/keyval/forcelink.h"

#if TA_DEFAULT_POLICY == 1
template class mpqc::lcao::DBCCSD_F12<TA::TensorD>;
MPQC_CLASS_EXPORT2("DBCCSD(F12)", mpqc::lcao::DBCCSD_F12<TA::TensorD>);
#endif
