#include <mpqc/chemistry/qc/wfn/wfn_world.h>

MPQC_CLASS_EXPORT_KEY2(mpqc::qc::WfnWorld, "WfnWorld");

using namespace mpqc;
qc::WfnWorld::WfnWorld(KeyVal const &kv)
    : world_(*kv.value<madness::World *>("$:world")), ao_ints_(kv) {
  mol_ = kv.keyval("molecule").class_ptr<molecule::Molecule>();
  basis_registry_ = ao_ints_.orbital_basis_registry();
}
qc::WfnWorld::~WfnWorld() = default;
