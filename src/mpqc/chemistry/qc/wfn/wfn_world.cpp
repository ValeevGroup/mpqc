#include "mpqc/chemistry/qc/wfn/wfn_world.h"

using namespace mpqc;
lcao::WavefunctionWorld::WavefunctionWorld(KeyVal const &kv)
    : world_(*kv.value<madness::World *>("$:world"))
{
  atoms_ = kv.keyval("molecule").class_ptr<Molecule>();
  basis_registry_ = std::make_shared<gaussian::OrbitalBasisRegistry>(kv);
}
lcao::WavefunctionWorld::~WavefunctionWorld() = default;
