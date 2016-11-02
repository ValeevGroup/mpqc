//
// Created by Chong Peng on 11/1/16.
//

#include <mpqc/chemistry/qc/f12/dbccsdf12.h>

namespace mpqc {
namespace f12 {

template class DBCCSDF12<TA::TensorD>;

template <>
DBCCSDF12<TA::TensorD>::~DBCCSDF12() = default;

//template class DBCCSDF12<TA::TensorD>;

}
}
