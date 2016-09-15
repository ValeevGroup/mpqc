//
// Created by Chong Peng on 4/26/16.
//

#include <mpqc/chemistry/qc/integrals/lcao_factory.h>
#include <mpqc/chemistry/qc/integrals/direct_atomic_integral.h>

namespace mpqc {
namespace integrals {

template class AtomicIntegral<TA::TensorD, TA::SparsePolicy>;

template class DirectAtomicIntegral<TA::TensorD, TA::SparsePolicy>;

template class LCAOFactory<TA::TensorD, TA::SparsePolicy>;
}
}