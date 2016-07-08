//
// task_integrals_helper.h
//
// Copyright (C) 2015 Drew Lewis
// Maintainer Drew Lewis
//

#ifndef MPQC_INTEGRALS_TASKINTEGRALSHELPER_H
#define MPQC_INTEGRALS_TASKINTEGRALSHELPER_H

#include "../../../../../common/typedefs.h"
#include "../../../../../include/tiledarray.h"
#include <mpqc/chemistry/qc/integrals/task_integrals_common.h>

#include <mpqc/chemistry/qc/integrals/screening/screen_base.h>

#include <libint2/engine.h>

namespace mpqc {
namespace integrals {
namespace detail {

extern double integral_engine_precision;

using Engine = libint2::Engine;

inline void set_eng_precision(Engine &eng) {
  eng.set_precision(integral_engine_precision);
}

template <typename... Shells>
inline const auto& shell_set(Engine &e, Shells &&... shells)
    -> decltype(e.compute(std::forward<Shells>(shells)...)) {
  return e.compute(std::forward<Shells>(shells)...);
}

TA::TensorD integral_kernel(Engine &eng, TA::Range &&rng,
                            std::array<ShellVec const *, 2> shell_ptrs,
                            Screener &);

TA::TensorD integral_kernel(Engine &eng, TA::Range &&rng,
                            std::array<ShellVec const *, 3> shell_ptrs,
                            Screener &screen);

TA::TensorD integral_kernel(Engine &eng, TA::Range &&rng,
                            std::array<ShellVec const *, 4> shell_ptrs,
                            Screener &screen);

}  // namespace detail
}  // namespace integrals
}  // namespace mpqc

#endif  // MPQC_INTEGRALS_TASKINTEGRALSHELPER_H
