//
// Created by Chong Peng on 10/6/16.
//

#include "mp2.h"
#include "mpqc/util/keyval/forcelink.h"

#if TA_DEFAULT_POLICY == 0

template class mpqc::lcao::RMP2<TA::TensorD, TA::DensePolicy>;
MPQC_CLASS_EXPORT2("RMP2", mpqc::lcao::RMP2<TA::TensorD, TA::DensePolicy>);

template class mpqc::lcao::RIRMP2<TA::TensorD, TA::DensePolicy>;
MPQC_CLASS_EXPORT2("RI-RMP2", mpqc::lcao::RIRMP2<TA::TensorD, TA::DensePolicy>);

#elif TA_DEFAULT_POLICY == 1
template class mpqc::lcao::RMP2<TA::TensorD, TA::SparsePolicy>;
MPQC_CLASS_EXPORT2("RMP2", mpqc::lcao::RMP2<TA::TensorD, TA::SparsePolicy>);

template class mpqc::lcao::RIRMP2<TA::TensorD, TA::SparsePolicy>;
MPQC_CLASS_EXPORT2("RI-RMP2", mpqc::lcao::RIRMP2<TA::TensorD, TA::SparsePolicy>);

#endif
