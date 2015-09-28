// Copyright (C) 2015 Drew Lewis
// Maintainer Drew Lewis
//
// Based on file task_integrals.hpp from mpqc
//

#pragma once

#ifndef MPQC_INTEGRALS_SCREENEDTASKINTEGRALS_H
#define MPQC_INTEGRALS_SCREENEDTASKINTEGRALS_H

#include "task_integrals_common.h"

#include "../integrals/integral_engine_pool.h"

#include <memory>
#include <algorithm>

namespace mpqc {
namespace integrals {

namespace detail {

int64_t nshells_in_basis(std::vector<tcc::basis::ClusterShells> const &cs) {
    return std::accumulate(cs.begin(), cs.end(), 0,
                           [](int nsh, tcc::basis::ClusterShells const &c) {
        return nsh + c.nshells();
    });
}

struct ScreeningMats {
    MatrixD cluster_screening;
    MatrixD shell_screening;
};

// Depends on the integrals being 1. DF, 2 obs, 3 obs
ScreeningMats
screening_matrix_X(ShrPool<TwoE_Engine> &engines,
                   std::vector<tcc::basis::ClusterShells> const &cs) {

    ScreeningMats sc_mats;
    // value for clusters
    sc_mats.cluster_screening = VectorD(cs.size());
    // value for shells
    sc_mats.shell_screening = VectorD(nshells_in_basis(cs));

    // TODO try and fix this
    engines->set_precision(0.);

    const auto unit = Shell::unit();

    // Loop over clusters
    auto shell_count = 0;
    auto cl_count = 0;
    for (auto const &cl : cs) {
        auto cl_norm2 = 0;

        // Loop over shells
        for (auto const &sh : cl.flattened_shells()) {
            auto nsh = sh.size();
            const auto *buf = engines->local().compute(sh, unit, sh, unit);

            auto sh_norm2 = 0.0;
            for (auto i = 0ul; i < nsh * nsh; ++i) {
                const auto bufi2 = buf[i] * buf[i];
                sh_norm2 += bufi2;
                cl_norm2 += bufi2;
            }
            sc_mats.shell_screening[shell_count++] = std::sqrt(sh_norm2);
        }
        sc_mats.cluster_screening[cl_count++] = std::sqrt(cl_norm2);
    }

    return sc_mats;
}

ScreeningMats
screening_matrix_ab(ShrPool<TwoE_Engine> &engines,
                    std::vector<tcc::basis::ClusterShells> const &cs) {
    // value for clusters
    ScreeningMats sc_mats;
    sc_mats.cluster_screening = MatrixD(cs.size(), cs.size());

    // value for shells
    const auto nshells = nshells_in_basis(cs);
    sc_mats.shell_screening = MatrixD(nshells, nshells);

    // TODO pick back up here tomorrow, this part will need focus
    for (auto const &cl : cs) {
        auto const &cl_shells = cl.flattened_shells();
        const auto cl_nshells = cl_shells.size();

        for(auto i = 0ul; i < cl_nshells; ++i){
            auto const &sh0 = cl_shells[i];

            for(auto j = 0ul; j <= i; ++j){
                auto const &sh1 = cl_shells[j];

                const auto *buf = engines->local().compute(sh0, sh1, sh0, sh1);

                // extract ints into an Eigen Matrix
                Matrix shblk = Matrix::Zero(n1, n2);
                for (size_t f1 = 0, f12 = 0; f1 != n1; ++f1)
                    for (size_t f2 = 0; f2 != n2; ++f2, ++f12) {
                        const auto int1212 = buf[f12 * n12 + f12];
                        shblk(f1, f2) = int1212;
                    }

                K(s1, s2) = K(s2, s1)
                      = std::sqrt(shblk.lpNorm<Eigen::Infinity>());
            }
        }
    }

    return sc_mats;
}

template <typename E, typename Op>
DArray<3, Ttype<Op>, SpPolicy>
compute_screened_integrals(mad::World &world, ShrPool<E> &engines,
                           Barray<3> const &bases, Op op) {
    // Depends on auxilary basis being in position 0.
    auto Q_X = screening_matrix_X(engines, bases[0]);

    // Assumes basis 1 and basis 2 are the same
    auto Q_ab = screening_matrix_ab(engines, bases[1]);
}

template <typename E, typename Op>
DArray<4, Ttype<Op>, SpPolicy>
compute_screened_integrals(mad::World &world, ShrPool<E> &engines,
                           Barray<4> const &bases, Op op) {
    // Assumes basis 1 and basis 2 are the same
    auto Q_ab = screening_matrix_ab(engines, bases[0]);
}

} // namespace detail

/*! \brief Construct integral tensors in parallel with screening.
 *
 */
template <typename E, unsigned long N, typename Op>
DArray<N, detail::Ttype<Op>, SpPolicy>
ScreenedTaskInts(mad::World &world, ShrPool<E> const &engines,
                 Barray<N> const &bases, Op op) {
    static_assert(N == 3 || N == 4,
                  "Screening only avalible for 3 and 4 center integrals.");
    return detail::compute_screened_integrals(world, engines, bases, op);
}


} // namespace integrals
} // namespace mpqc


#endif // MPQC_INTEGRALS_SCREENEDTASKINTEGRALS_H
