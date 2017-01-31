//
// Created by Chong Peng on 10/6/16.
//

#include "lcao_wfn.h"

namespace mpqc {
namespace lcao {


#if TA_DEFAULT_POLICY == 0
template class LCAOWavefunction<TA::TensorD, TA::DensePolicy>;
#elif TA_DEFAULT_POLICY == 1
template class LCAOWavefunction<TA::TensorD, TA::SparsePolicy>;
template class PeriodicLCAOWavefunction<TA::TensorZ, TA::SparsePolicy>;
#endif

}
}
