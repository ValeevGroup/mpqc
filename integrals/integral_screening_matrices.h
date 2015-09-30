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
#include "../include/eigen.h"

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

struct ScreeningMatrices {
    std::vector<std::vector<MatrixD>> shell_screenings;
    MatrixD cluster_screening;
};

// Depends on the integrals being 1. DF, 2 obs, 3 obs
ScreeningMatrices
screening_matrix_X(ShrPool<TwoE_Engine> &engines,
                   std::vector<tcc::basis::ClusterShells> const &cs) {

    const auto nclusters = cs.size();
    ScreeningMatrices sc_mats;
    sc_mats.cluster_screening = VectorD(nclusters);

    auto &cl_mat = sc_mats.cluster_screening;
    auto &sh_vec = sc_mats.shell_screenings;
    sh_vec.reserve(nclusters);

    auto &eng = engines->local();
    eng.set_precision(0.);

    const auto unit = Shell::unit();

    // Loop over clusters
    for (auto c = 0ul; c < nclusters; ++c) {
        auto const &cl_shells = cs[c].flattened_shells();
        const auto cl_size = cl_shells.size();

        sh_vec.emplace_back(std::vector<MatrixD>{VectorD(cl_size)});
        auto &sh_mat = sh_vec.back().back();

        for (auto s = 0ul; s < cl_size; ++s) {

            auto const &sh = cl_shells[s];
            auto nsh = sh.size();
            const auto *buf = engines->local().compute(sh, unit, sh, unit);

            const auto bmap = Eig::Map<const MatrixD>(buf, nsh, nsh);
            sh_mat(s) = std::sqrt(bmap.lpNorm<2>());
        }

        cl_mat(c) = sh_mat.norm();
    }

    return sc_mats;
}

ScreeningMatrices
screening_matrix_ab(ShrPool<TwoE_Engine> &engines,
                    std::vector<tcc::basis::ClusterShells> const &cs) {

    const auto nclusters = cs.size();

    ScreeningMatrices sc_mats;
    sc_mats.cluster_screening = MatrixD(nclusters, nclusters);
    auto &cl_mat = sc_mats.cluster_screening;

    // value for shells
    auto &sh_vecs = sc_mats.shell_screenings;
    sh_vecs.reserve(nclusters);

    auto &eng = engines->local();
    eng.set_precision(0.);

    for (auto c0 = 0ul; c0 < nclusters; ++c0) {
        auto const &shells0 = cs[c0].flattened_shells();
        const auto nshells0 = shells0.size();

        sh_vecs.emplace_back(std::vector<MatrixD>{});
        auto &current_vec = sh_vecs.back();
        current_vec.reserve(nclusters);

        for (auto c1 = 0ul; c1 < nclusters; ++c1) {
            auto const &shells1 = cs[c1].flattened_shells();
            const auto nshells1 = shells1.size();
            current_vec.emplace_back(MatrixD(nshells0, nshells1));
            auto &sh_mat = current_vec.back();

            for (auto s0 = 0ul; s0 < nshells0; ++s0) {
                auto const &sh0 = shells0[s0];
                const auto nsh0 = sh0.size();

                for (auto s1 = 0ul; s1 < nshells1; ++s1) {
                    auto const &sh1 = shells1[s1];
                    const auto nsh1 = sh1.size();

                    const auto *buf = eng.compute(sh0, sh1, sh0, sh1);
                    const auto bmap = Eig::Map<const MatrixD>(buf, nsh0 * nsh1,
                                                              nsh0 * nsh1);

                    sh_mat(s0, s1) = std::sqrt(bmap.lpNorm<2>());
                }
            }
            cl_mat(c0, c1) = std::sqrt(sh_mat.lpNorm<2>());
        }
    }

    return sc_mats;
}

} // namespace detail
} // namespace integrals
} // namespace mpqc

#endif //  MPQC_INTEGRALS_INTEGRALSCREENINGMATRICES_H
