//
// Created by Chong Peng on 11/2/16.
//

#include "mpqc/chemistry/qc/cc/dbccsd.h"
#include "mpqc/util/keyval/forcelink.h"

#if TA_DEFAULT_POLICY == 0
template class mpqc::cc::DBCCSD<TA::TensorD, TA::DensePolicy>;
MPQC_CLASS_EXPORT2("DBCCSD", mpqc::cc::DBCCSD<TA::TensorD, TA::DensePolicy>);
#elif TA_DEFAULT_POLICY == 1
template class mpqc::cc::DBCCSD<TA::TensorD, TA::SparsePolicy>;
MPQC_CLASS_EXPORT2("DBCCSD", mpqc::cc::DBCCSD<TA::TensorD, TA::SparsePolicy>);
#endif