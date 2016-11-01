//
// Created by Chong Peng on 11/1/16.
//

#include <mpqc/chemistry/qc/f12/gf2f12.h>

MPQC_CLASS_EXPORT_KEY2("GF2F12", mpqc::f12::GF2F12<TA::TensorD>);

namespace mpqc {
namespace f12 {

template class GF2F12<TA::TensorD>;
template <>
GF2F12<TA::TensorD>::~GF2F12() = default;


}
}