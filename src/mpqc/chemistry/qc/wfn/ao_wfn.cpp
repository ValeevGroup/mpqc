#include "mpqc/chemistry/qc/properties/propertybase.h"
#include "mpqc/chemistry/qc/wfn/ao_wfn.h"

namespace mpqc{
namespace qc{


template class AOWavefunction<TA::TensorD, TA::SparsePolicy>;
template class AOWavefunction<TA::TensorD, TA::DensePolicy>;
template class PeriodicAOWavefunction<TA::TensorZ, TA::SparsePolicy>;

}//namespace qc
}//namespace mpqc

