#include "mpqc/chemistry/qc/lcao/wfn/ao_wfn.h"
#include "mpqc/chemistry/qc/properties/property.h"

namespace mpqc {
namespace lcao {

#if TA_DEFAULT_POLICY == 0
template class AOWavefunction<TA::TensorD, TA::DensePolicy>;
#elif TA_DEFAULT_POLICY == 1
template class AOWavefunction<TA::TensorD, TA::SparsePolicy>;
template class PeriodicAOWavefunction<TA::TensorZ, TA::SparsePolicy>;
#endif

}  // namespace lcao
}  // namespace mpqc
