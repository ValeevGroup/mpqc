//
// task_integrals_helper.h
//
// Copyright (C) 2015 Drew Lewis
// Maintainer Drew Lewis
//

#pragma once

#ifndef MPQC_INTEGRALS_TASKINTEGRALSHELPER_H
#define MPQC_INTEGRALS_TASKINTEGRALSHELPER_H

#include "../common/typedefs.h"
#include "../include/tiledarray.h"
#include "../include/libint.h"

//  #include "../basis/basis.h"
// #include "../basis/basis_fwd.h"

namespace mpqc {
namespace integrals {
namespace detail {

using OneE_Engine = libint2::OneBodyEngine;
using TwoE_Engine = libint2::TwoBodyEngine<libint2::Coulomb>;

inline const double *shell_set(TwoE_Engine &e, Shell const &s0, Shell const &s1,
                               Shell const &s2, Shell const &s3);

inline const double *
shell_set(TwoE_Engine &e, Shell const &s0, Shell const &s1); 

inline const double *
shell_set(TwoE_Engine &e, Shell const &s0, Shell const &s1, Shell const &s2);

inline const double *
shell_set(OneE_Engine &e, Shell const &s0, Shell const &s1);

template <typename Engine>
TA::TensorD integral_kernel(Engine &eng, TA::Range &&rng,
                            std::array<ShellVec const *, 2> shell_ptrs); 

template <typename Engine>
TA::TensorD integral_kernel(Engine &eng, TA::Range &&rng,
                            std::array<ShellVec const *, 3> shell_ptrs,
                           MatrixD const &X, MatrixD const &ab); 

template <typename Engine>
TA::TensorD integral_kernel(Engine &eng, TA::Range &&rng,
                            std::array<ShellVec const *, 3> shell_ptrs); 

template <typename Engine>
TA::TensorD integral_kernel(Engine &eng, TA::Range &&rng,
                            std::array<ShellVec const *, 4> shell_ptrs);


} // namespace detail
} // namespace integrals
} // namespace mpqc

#endif // MPQC_INTEGRALS_TASKINTEGRALSHELPER_H
