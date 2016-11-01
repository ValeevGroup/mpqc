//
// Created by Chong Peng on 11/1/16.
//

#include <mpqc/chemistry/qc/f12/dbccsdf12.h>

MPQC_CLASS_EXPORT_KEY2("DBCCSD(F12)", mpqc::f12::DBCCSDF12<TA::TensorD>);

namespace mpqc {
namespace f12 {

template class DBCCSDF12<TA::TensorD>;

template <>
DBCCSDF12<TA::TensorD>::~DBCCSDF12() = default;

//template class DBCCSDF12<TA::TensorD>;

}
}