#include "mpqc/chemistry/qc/lcao/wfn/periodic_ao_wfn.h"
#include "mpqc/chemistry/qc/properties/property.h"

namespace mpqc {
namespace lcao {

#if TA_DEFAULT_POLICY == 1
template class PeriodicAOWavefunction<TA::TensorD, TA::SparsePolicy>;
#endif

}  // namespace lcao
}  // namespace mpqc
