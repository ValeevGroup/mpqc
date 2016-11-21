//
// Created by Chong Peng on 4/26/16.
//

#include "mpqc/chemistry/qc/integrals/direct_ao_factory.h"
#include "mpqc/chemistry/qc/integrals/lcao_factory.h"

namespace mpqc {
namespace integrals {

template class AOFactory<TA::TensorD, TA::SparsePolicy>;

template class DirectAOFactory<TA::TensorD, TA::SparsePolicy>;

template class LCAOFactory<TA::TensorD, TA::SparsePolicy>;
}
}