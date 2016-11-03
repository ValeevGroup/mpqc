#include <mpqc/chemistry/qc/wfn/wfn_world.h>

using namespace mpqc;
qc::WavefunctionWorld::WavefunctionWorld(KeyVal const &kv)
    : world_(*kv.value<madness::World *>("$:world"))
{
  mol_ = kv.keyval("molecule").class_ptr<Molecule>();
  basis_registry_ = std::make_shared<basis::OrbitalBasisRegistry>(kv);
}
qc::WavefunctionWorld::~WavefunctionWorld() = default;
