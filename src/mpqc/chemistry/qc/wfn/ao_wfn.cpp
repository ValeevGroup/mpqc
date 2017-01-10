#include "mpqc/chemistry/qc/wfn/ao_wfn.h"
#include "mpqc/chemistry/qc/properties/property.h"

namespace mpqc {
namespace lcao {

template class AOWavefunction<TA::TensorD, TA::SparsePolicy>;
template class AOWavefunction<TA::TensorD, TA::DensePolicy>;
template class PeriodicAOWavefunction<TA::TensorZ, TA::SparsePolicy>;

}  // namespace lcao
}  // namespace mpqc
