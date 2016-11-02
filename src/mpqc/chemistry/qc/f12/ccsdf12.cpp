//
// Created by Chong Peng on 11/1/16.
//

#include <mpqc/chemistry/qc/f12/ccsdf12.h>

namespace mpqc {
namespace f12 {

template class CCSDF12<TA::TensorD>;
template <>
CCSDF12<TA::TensorD>::~CCSDF12() = default;


}
}
