#include <mpqc/chemistry/qc/properties/propertybase.h>
#include <mpqc/chemistry/qc/wfn/ao_wfn.h>

namespace mpqc{
namespace qc{

AOWavefunction::AOWavefunction(const KeyVal &kv) : Wavefunction(kv)
{
  ao_ints_ = integrals::detail::construct_atomic_integral<TA::TensorD, TA::SparsePolicy>(kv);
  ao_ints_->set_orbital_basis_registry(this->wfn_world()->basis_registry());

}

AOWavefunction::~AOWavefunction() = default;

void AOWavefunction::compute(qc::PropertyBase *pb) { pb->apply(this); }

void AOWavefunction::obsolete() {
  ao_integrals().registry().purge(wfn_world()->world());
}

}//namespace qc
}//namespace mpqc

