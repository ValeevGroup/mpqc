#include <mpqc/chemistry/qc/properties/propertybase.h>
#include <mpqc/chemistry/qc/wfn/ao_wfn.h>

namespace mpqc{
namespace qc{

AOWavefunction::AOWavefunction(KeyVal const &kv) : Wavefunction(kv) {}

AOWavefunction::~AOWavefunction() = default;

void AOWavefunction::compute(qc::PropertyBase *pb) { pb->apply(this); }

void AOWavefunction::obsolete() {
  wfn_world()->ao_integrals().registry().purge(wfn_world()->world());
}

}//namespace qc
}//namespace mpqc

