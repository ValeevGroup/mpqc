//
// Created by Chong Peng on 11/1/16.
//

#include "ccsd_t.h"

MPQC_CLASS_EXPORT_KEY2("CCSD(T)", mpqc::cc::CCSD_T<TA::TensorD, TA::SparsePolicy>);

namespace mpqc {
namespace cc {


template class CCSD_T<TA::TensorD, TA::SparsePolicy>;
template <>
CCSD_T<TA::TensorD, TA::SparsePolicy>::~CCSD_T() = default;

// template class DBCCSD<TA::TensorD, TA::SparsePolicy>;
}
}