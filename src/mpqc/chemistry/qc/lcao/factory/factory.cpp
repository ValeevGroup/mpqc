//
// Created by Chong Peng on 4/26/16.
//

#include "mpqc/chemistry/qc/lcao/factory/lcao_factory.h"
#include "mpqc/chemistry/qc/lcao/factory/periodic_lcao_factory.h"

namespace mpqc {
namespace lcao {
namespace gaussian {

#if TA_DEFAULT_POLICY == 0
template class AOFactory<TA::TensorD, TA::DensePolicy>;
#elif TA_DEFAULT_POLICY == 1
template class AOFactory<TA::TensorD, TA::SparsePolicy>;
template class PeriodicAOFactory<TA::TensorZ, TA::SparsePolicy>;
#endif

}  // namespace gaussian

#if TA_DEFAULT_POLICY == 0
template class LCAOFactory<TA::TensorD, TA::DensePolicy>;
#elif TA_DEFAULT_POLICY == 1
template class LCAOFactory<TA::TensorD, TA::SparsePolicy>;
template class PeriodicLCAOFactory<TA::TensorZ, TA::SparsePolicy>;
#endif

}  // namespace lcao
}  // namespace mpqc
