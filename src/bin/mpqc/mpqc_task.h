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

namespace mpqc {

/// \brief An MPQC computation
///
/// A computation is specified by a KeyVal object and a World object
class MPQCTask {
 public:
  MPQCTask(madness::World &world, std::shared_ptr<KeyVal> kv);
  ~MPQCTask();

  madness::World& world() const;
  const std::shared_ptr<KeyVal>& keyval() const;

  void run();

 private:
  madness::World& world_;
  std::shared_ptr<KeyVal> keyval_;
};

}  // namespace mpqc


#endif  // MPQC4_SRC_BIN_MPQC_MPQC_TASK_H_
