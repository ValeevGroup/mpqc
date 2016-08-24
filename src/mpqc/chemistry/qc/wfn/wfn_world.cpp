#include <mpqc/chemistry/qc/wfn/wfn_world.h>

MPQC_CLASS_EXPORT_KEY2(mpqc::qc::WfnWorld, "WfnWorld");

using namespace mpqc;
qc::WfnWorld::WfnWorld(KeyVal const &kv)
    : world_(*kv.value<madness::World *>("world")),
      ao_ints_(
          *kv.value<integrals::AtomicIntegral<TA::TensorD, TA::SparsePolicy> *>(
              "ao_integrals")) {
  mol_ = kv.class_ptr<molecule::Molecule>("molecule");
}
