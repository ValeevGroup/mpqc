#include <mpqc/chemistry/qc/wfn/wfn_world.h>

MPQC_CLASS_EXPORT_KEY2(mpqc::qc::WfnWorld, "WfnWorld");

using namespace mpqc;
qc::WfnWorld::WfnWorld(KeyVal const &kv) {
  mol_ = kv.class_ptr<molecule::Molecule>("molecule");
//   ao_ints_ =
//       kv.class_ptr<integrals::AtomicIntegral<TA::TensorD, TA::SparseShape>>(
//           "ao_integrals");
}
