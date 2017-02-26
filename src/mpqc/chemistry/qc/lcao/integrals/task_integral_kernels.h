//
// task_integrals_helper.h
//
// Copyright (C) 2015 Drew Lewis
// Maintainer Drew Lewis
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_TASK_INTEGRAL_KERNELS_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_TASK_INTEGRAL_KERNELS_H_

#include <libint2/engine.h>
#include <tiledarray.h>


#include "mpqc/chemistry/qc/lcao/integrals/task_integrals_common.h"
#include "mpqc/chemistry/qc/lcao/integrals/screening/screen_base.h"

namespace mpqc {
namespace lcao {
namespace gaussian {
namespace detail {

extern double integral_engine_precision;

using Engine = libint2::Engine;
using data_pointer = TA::TensorD::pointer;

inline void set_eng_precision(Engine &eng) {
  eng.set_precision(integral_engine_precision);
}

template <typename... Shells>
inline void shell_set(Engine &e, Shells &&... shells) {
  e.compute(std::forward<Shells>(shells)...);
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
}  // namespace gaussian
}  // namespace lcao
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_TASK_INTEGRAL_KERNELS_H_
