//
// Copyright (C) 2015 Drew Lewis
// Maintainer Drew Lewis
//
//

#pragma once

#ifndef MPQC_INTEGRALS_INTEGRALSCREENINGMATRICES_H
#define MPQC_INTEGRALS_INTEGRALSCREENINGMATRICES_H

#include "task_integrals_common.h"

#include "../integrals/integral_engine_pool.h"
#include "../include/tiledarray.h"
#include "../include/eigen.h"


namespace mpqc {
namespace integrals {

namespace detail {

struct ScreeningMatrices {
    std::vector<std::vector<MatrixD>> shell_screenings;
    MatrixD cluster_screening;
};

// Depends on the integrals being 1. DF, 2 obs, 3 obs
ScreeningMatrices
screening_matrix_X(mad::World &world, ShrPool<TwoE_Engine> &engines,
                   basis::Basis const &basis);

ScreeningMatrices
screening_matrix_ab(mad::World &world, ShrPool<TwoE_Engine> &engines,
                    basis::Basis const &basis);

} // namespace detail
} // namespace integrals
} // namespace mpqc 

#endif //  MPQC_INTEGRALS_INTEGRALSCREENINGMATRICES_H
