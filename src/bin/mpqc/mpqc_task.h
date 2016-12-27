/*
 * mpqc_task.h
 *
 *  Created on: Nov 14, 2016
 *      Author: evaleev
 */

#ifndef MPQC4_SRC_BIN_MPQC_MPQC_TASK_H_
#define MPQC4_SRC_BIN_MPQC_MPQC_TASK_H_

#include <memory>

#include <madness/world/world.h>

#include "mpqc/util/keyval/keyval.h"
#include "mpqc/chemistry/qc/wfn/wfn.h"

namespace mpqc {

/// \brief An MPQC computation
///
/// A computation is specified by a KeyVal object and a World object
class MPQCTask {
 public:
  MPQCTask(madness::World &world, std::shared_ptr<KeyVal> kv)
      : world_(world), keyval_(kv) {}
  ~MPQCTask() = default;

  madness::World& world() const {
    return world_;
  }
  const std::shared_ptr<KeyVal>& keyval() const {
    return keyval_;
  }

  void run() {
    auto world_popper = TA::push_default_world(world_);

    auto wfn = keyval_->keyval("wfn").class_ptr<lcao::Wavefunction>();

    //  auto energy_prop = qc::Energy(kv);
    //  auto energy_prop_ptr = &energy_prop;

    // set the sparse_threshold
    const double threshold = keyval_->value<double>("sparse_threshold", 1e-20);
    TiledArray::SparseShape<float>::threshold(threshold);

    const auto value = wfn->value();
    keyval_->assign("wfn:energy", value);
  }

 private:
  madness::World& world_;
  std::shared_ptr<KeyVal> keyval_;
};

}  // namespace mpqc


#endif  // MPQC4_SRC_BIN_MPQC_MPQC_TASK_H_
