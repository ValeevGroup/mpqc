//
// Created by Eduard Valeyev on 8/3/18.
//

#include "mpqc_task.h"

#include "mpqc/chemistry/qc/lcao/wfn/wfn.h"
#include "mpqc/chemistry/qc/properties/property.h"

namespace mpqc {

MPQCTask::MPQCTask(madness::World& world, std::shared_ptr<KeyVal> kv)
    : world_(world), keyval_(kv) {}
MPQCTask::~MPQCTask() {}

madness::World& MPQCTask::world() const { return world_; }
const std::shared_ptr<KeyVal>& MPQCTask::keyval() const { return keyval_; }

void MPQCTask::run() {
  auto world_popper = TA::push_default_world(world_);

  // set the sparse_threshold
  const double threshold = keyval_->value<double>("sparse_threshold", 1e-20);
  TiledArray::SparseShape<float>::threshold(threshold);

  auto property = keyval_->class_ptr<Property>("property");
  if (property != nullptr) {
    property->evaluate();
    auto kv_prop = keyval_->keyval("property");
    property->write(kv_prop);
  } else {
    throw InputError("invalid property", __FILE__, __LINE__, "property");
  }
}

}  // namespace mpqc
