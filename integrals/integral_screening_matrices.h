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

    auto &eng = engines->local();
    eng.set_precision(0.);

    const auto unit = Shell::unit();

    // Loop over clusters
    auto shell_count = 0;
    auto cl_count = 0;
    for (auto const &cl : cs) {
        auto cl_norm = 0.0;

        // Loop over shells
        for (auto const &sh : cl.flattened_shells()) {
            auto nsh = sh.size();
            const auto *buf = engines->local().compute(sh, unit, sh, unit);

            const auto sh_norm = *std::max_element(buf, buf + nsh * nsh);
            cl_norm = std::max(sh_norm, cl_norm);
            sc_mats.shell_screening(shell_count++) = std::sqrt(sh_norm);
        }
        sc_mats.cluster_screening(cl_count++) = std::sqrt(cl_norm);
    }

    return sc_mats;
}

ScreeningMats
screening_matrix_ab(ShrPool<TwoE_Engine> &engines,
                    std::vector<tcc::basis::ClusterShells> const &cs) {
    // value for clusters
    ScreeningMats sc_mats;
    sc_mats.cluster_screening = MatrixD(cs.size(), cs.size());
    auto &cl_mat = sc_mats.cluster_screening;

    // value for shells
    const auto nshells = nshells_in_basis(cs);
    sc_mats.shell_screening = MatrixD(nshells, nshells);
    auto &shell_mat = sc_mats.shell_screening;

    auto &eng = engines->local();
    eng.set_precision(0.);

    // Shell ord goes outside cluster loop because it is keeping track of
    // global shell index

    auto s0_start = 0;
    for (auto cl0 = 0ul; cl0 < cs.size(); ++cl0) {
        auto const &shells0 = cs[cl0].flattened_shells();

        auto cluster_norm = 0.0;
        auto s1_start = 0;
        for (auto cl1 = 0ul; cl1 <= cl0; ++cl1) {
            auto const &shells1 = cs[cl1].flattened_shells();

            for (auto s0 = 0ul; s0 < shells0.size(); ++s0) {
                auto const &sh0 = shells0[s0];
                const auto nsh0 = sh0.size();
                const auto s0_pos = s0 + s0_start;

                for (auto s1 = 0ul; s1 < shells1.size(); ++s1) {
                    auto const &sh1 = shells1[s1];
                    const auto nsh1 = sh1.size();
                    const auto s1_pos = s1 + s1_start;


                    const auto *buf = eng.compute(sh0, sh1, sh0, sh1);
                    const auto bmap = Eig::Map<const MatrixD>(buf, nsh0, nsh1);

                    const auto norm = bmap.lpNorm<Eigen::Infinity>();
                    cluster_norm = std::max(norm, cluster_norm);
                    const auto sqrt_norm = std::sqrt(norm);
                    shell_mat(s0_pos, s1_pos) = sqrt_norm;
                    shell_mat(s1_pos, s0_pos) = sqrt_norm;
                }
            }
            s1_start += shells1.size();

            const auto sqrt_cl_norm = std::sqrt(cluster_norm);
            cl_mat(cl0, cl1) = sqrt_cl_norm;
            cl_mat(cl1, cl0) = sqrt_cl_norm;
        }
        s0_start += shells0.size();
    }

    return sc_mats;
}

} // namespace detail
} // namespace integrals
} // namespace mpqc

#endif //  MPQC_INTEGRALS_INTEGRALSCREENINGMATRICES_H

