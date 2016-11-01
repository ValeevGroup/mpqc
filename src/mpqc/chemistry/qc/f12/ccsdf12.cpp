//
// Created by Chong Peng on 11/1/16.
//

#include <mpqc/chemistry/qc/f12/ccsdf12.h>

MPQC_CLASS_EXPORT_KEY2("CCSD(F12)", mpqc::f12::CCSDF12<TA::TensorD>);

namespace mpqc {
namespace f12 {

template class CCSDF12<TA::TensorD>;
template <>
CCSDF12<TA::TensorD>::~CCSDF12() = default;


}
}