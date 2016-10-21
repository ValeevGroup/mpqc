//
// Created by Chong Peng on 4/26/16.
//

#include <mpqc/chemistry/qc/integrals/lcao_factory.h>
#include <mpqc/chemistry/qc/integrals/direct_atomic_integral.h>


MPQC_CLASS_EXPORT_KEY2("AtomicIntegralTensorD", mpqc::integrals::AtomicIntegral<TA::TensorD,TA::SparsePolicy>);
MPQC_CLASS_EXPORT_KEY2("LCAOFactoryTensorD", mpqc::integrals::LCAOFactory<TA::TensorD,TA::SparsePolicy>);
MPQC_CLASS_EXPORT_KEY2("DirectAtomicIntegralTensorD", mpqc::integrals::DirectAtomicIntegral<TA::TensorD,TA::SparsePolicy>);

namespace mpqc {
namespace integrals {


template class AtomicIntegral<TA::TensorD, TA::SparsePolicy>;

template class DirectAtomicIntegral<TA::TensorD, TA::SparsePolicy>;

template class LCAOFactory<TA::TensorD, TA::SparsePolicy>;
}
}