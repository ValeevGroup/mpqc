// Copyright (C) 2015 Drew Lewis
// Maintainer Drew Lewis
//
// Based on file task_integrals.hpp from mpqc
//

#pragma once

#ifndef MPQC_INTEGRALS_SCREENEDTASKINTEGRALS_H
#define MPQC_INTEGRALS_SCREENEDTASKINTEGRALS_H

#include "task_integrals_common.h"
#include "integral_screening_matrices.h"

#include "../integrals/integral_engine_pool.h"
#include "../include/eigen.h"

#include <memory>
#include <algorithm>

namespace mpqc {
namespace integrals {

namespace detail {

template <typename E, typename Op>
DArray<3, Ttype<Op>, SpPolicy>
compute_screened_integrals(mad::World &world, ShrPool<E> &engines,
                           Barray<3> const &bases, Op op) {
    // Depends on auxilary basis being in position 0.
    const auto Q_X = screening_matrix_X(engines, bases[0].cluster_shells());
    const auto Q_ab = screening_matrix_ab(engines, bases[1].cluster_shells());

    auto trange = create_trange(bases);
    const auto tvolume = trange.tiles().volume();

    auto pmap = SpPolicy::default_pmap(world, tvolume);

    // TODO actually put tasks here
    for (auto const &ord : *pmap) {
        auto const &idx = trange.tiles().idx(ord);
        std::cout << "(" << idx[0] << ", " << idx[1] << ", " << idx[2]
                  << ") = " << Q_X.cluster_screening(idx[0]) * Q_ab.cluster_screening(idx[1], idx[2]) << std::endl;
    }
}

template <typename E, typename Op>
DArray<4, Ttype<Op>, SpPolicy>
compute_screened_integrals(mad::World &world, ShrPool<E> &engines,
                           Barray<4> const &bases, Op op) {
    // Assumes basis 1 and basis 2 are the same
    auto Q_ab = screening_matrix_ab(engines, bases[0].cluster_shells());
}

} // namespace detail

/*! \brief Construct integral tensors in parallel with screening.
 *
 */
template <typename E, unsigned long N, typename Op>
DArray<N, detail::Ttype<Op>, SpPolicy>
ScreenedTaskInts(mad::World &world, ShrPool<E> &engines, Barray<N> const &bases,
                 Op op) {
    static_assert(N == 3 || N == 4,
                  "Screening only avalible for 3 and 4 center integrals.");
    return detail::compute_screened_integrals(world, engines, bases, op);
}


} // namespace integrals
} // namespace mpqc


#endif // MPQC_INTEGRALS_SCREENEDTASKINTEGRALS_H
