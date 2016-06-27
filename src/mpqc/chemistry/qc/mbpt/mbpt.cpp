//
// Created by Chong Peng on 6/6/16.
//

#include <mpqc/chemistry/qc/mbpt/mp2.h>
#include <mpqc/chemistry/qc/mbpt/dbmp2.h>

namespace mpqc{
namespace mbpt{

template class MP2<TA::TensorD, TA::SparsePolicy>;
template class DBMP2<TA::TensorD, TA::SparsePolicy>;

}
}