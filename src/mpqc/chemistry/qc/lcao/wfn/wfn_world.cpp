#include "mpqc/chemistry/qc/lcao/wfn/wfn_world.h"

using namespace mpqc;
lcao::WavefunctionWorld::WavefunctionWorld(KeyVal const &kv)
    : world_(*kv.value<madness::World *>("$:world")) {
  atoms_ = kv.class_ptr<Molecule>("atoms");

  if (!atoms_)
    atoms_ = kv.class_ptr<Molecule>("molecule");
  if (!atoms_)
    throw InputError(
        "WavefunctionWorld did not find keyword \"atoms\"",
        __FILE__, __LINE__);

  basis_registry_ = std::make_shared<gaussian::OrbitalBasisRegistry>(kv);
}
