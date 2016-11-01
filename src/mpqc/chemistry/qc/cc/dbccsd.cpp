//
// Created by Chong Peng on 11/1/16.
//

#include "dbccsd.h"

MPQC_CLASS_EXPORT_KEY2("DBCCSD", mpqc::cc::DBCCSD<TA::TensorD, TA::SparsePolicy>);

namespace mpqc {
namespace cc {


template class DBCCSD<TA::TensorD, TA::SparsePolicy>;

template <>
DBCCSD<TA::TensorD, TA::SparsePolicy>::~DBCCSD() = default;

}
}