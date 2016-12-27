//
// Created by Chong Peng on 10/6/16.
//

#include "lcao_wfn.h"

namespace mpqc {
namespace lcao {

template class LCAOWavefunction<TA::TensorD, TA::DensePolicy>;
template class LCAOWavefunction<TA::TensorD, TA::SparsePolicy>;

}
}
