//
// Created by Chong Peng on 11/1/16.
//

#include <mpqc/chemistry/qc/f12/gf2f12.h>

namespace mpqc {
namespace f12 {

template class GF2F12<TA::TensorD>;
template <>
GF2F12<TA::TensorD>::~GF2F12() = default;


}
}
