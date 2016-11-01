//
// Created by Chong Peng on 6/6/16.
//

#include "ccsd.h"

MPQC_CLASS_EXPORT_KEY2("CCSD", mpqc::cc::CCSD<TA::TensorD, TA::SparsePolicy>);

namespace mpqc {
namespace cc {


template class CCSD<TA::TensorD, TA::SparsePolicy>;

template <>
CCSD<TA::TensorD, TA::SparsePolicy>::~CCSD() = default;

}
}
