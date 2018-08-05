
#include "mpqc/util/keyval/forcelink.h"
#include "mpqc/chemistry/qc/lcao/wfn/wfn.h"

namespace mpqc{
namespace lcao{

WavefunctionWorld::WavefunctionWorld(KeyVal const &kv)
    : ::mpqc::WavefunctionWorld(kv) {
  basis_registry_ = std::make_shared<gaussian::OrbitalBasisRegistry>(kv);
}

Wavefunction::Wavefunction(const KeyVal &kv) : ::mpqc::Wavefunction(kv) {
  // demand (and construct if needed) an LCAO WavefunctionWorld
  if (std::dynamic_pointer_cast<WavefunctionWorld>(this->::mpqc::Wavefunction::wfn_world()) == nullptr) {
    reset_wfn_world(std::make_shared<WavefunctionWorld>(kv));
  }
}

Wavefunction::~Wavefunction() { }

}  // namespace lcao
}  // namespace mpqc

MPQC_CLASS_EXPORT2("LCAOWfnWorld", mpqc::lcao::WavefunctionWorld);