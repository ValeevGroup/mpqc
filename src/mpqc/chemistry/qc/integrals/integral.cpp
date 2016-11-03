//
// Created by Chong Peng on 4/26/16.
//

#include <mpqc/chemistry/qc/integrals/direct_ao_factory.h>
#include <mpqc/chemistry/qc/integrals/lcao_factory.h>

MPQC_CLASS_EXPORT_KEY2(
    "AOFactoryTensorD",
    mpqc::integrals::AOFactory<TA::TensorD, TA::SparsePolicy>);
MPQC_CLASS_EXPORT_KEY2(
    "LCAOFactoryTensorD",
    mpqc::integrals::LCAOFactory<TA::TensorD, TA::SparsePolicy>);
MPQC_CLASS_EXPORT_KEY2(
    "DirectAOFactoryTensorD",
    mpqc::integrals::DirectAOFactory<TA::TensorD, TA::SparsePolicy>);

namespace mpqc {
namespace integrals {

template class AOFactory<TA::TensorD, TA::SparsePolicy>;

template class DirectAOFactory<TA::TensorD, TA::SparsePolicy>;

template class LCAOFactory<TA::TensorD, TA::SparsePolicy>;
}
}