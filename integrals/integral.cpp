//
// Created by Chong Peng on 4/26/16.
//

#include "./molecular_integral.h"


namespace mpqc{
namespace integrals{

    template class AtomicIntegral<TA::TensorD, TA::SparsePolicy>;

    template class MolecularIntegral<TA::TensorD, TA::SparsePolicy>;

}
}