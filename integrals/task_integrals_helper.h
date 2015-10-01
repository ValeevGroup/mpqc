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

#include "../basis/basis.h"
#include "../basis/cluster_shells.h"

namespace mpqc {
namespace integrals {
namespace detail {

using Shell = libint2::Shell;
static const auto unit_shell = Shell::unit();

using OneE_Engine = libint2::OneBodyEngine;
using TwoE_Engine = libint2::TwoBodyEngine<libint2::Coulomb>;

inline const double *shell_set(TwoE_Engine &e, Shell const &s0, Shell const &s1,
                               Shell const &s2, Shell const &s3) {
    return e.compute(s0, s1, s2, s3);
}

inline const double *
shell_set(TwoE_Engine &e, Shell const &s0, Shell const &s1) {
    return e.compute(s0, unit_shell, s1, unit_shell);
}

inline const double *
shell_set(TwoE_Engine &e, Shell const &s0, Shell const &s1, Shell const &s2) {
    return e.compute(s0, unit_shell, s1, s2);
}

inline const double *
shell_set(OneE_Engine &e, Shell const &s0, Shell const &s1) {
    return e.compute(s0, s1);
}

template <typename Engine>
TA::TensorD
integral_kernel(Engine &eng, TA::Range &&rng,
                std::array<tcc::basis::ClusterShells const *, 2> shell_ptrs) {

    auto const &sh0 = shell_ptrs[0]->flattened_shells();
    auto const &sh1 = shell_ptrs[1]->flattened_shells();

    auto const &lobound = rng.lobound();
    std::array<unsigned long, 2> lb = {{lobound[0], lobound[1]}};
    std::array<unsigned long, 2> ub = lb;

    auto tile = TA::TensorD(std::move(rng));

    for (auto const &s0 : sh0) {
        const auto ns0 = s0.size();
        ub[0] += ns0;

        lb[1] = ub[1] = lobound[1];
        for (auto const &s1 : sh1) {
            const auto ns1 = s1.size();
            ub[1] += ns1;

            tile.block(lb, ub) = TA::make_map(shell_set(eng, s0, s1), lb, ub);

            lb[1] = ub[1];
        }
        lb[0] = ub[0];
    }

    return tile;
}

template <typename Engine>
TA::TensorD
integral_kernel(Engine &eng, TA::Range &&rng,
                std::array<tcc::basis::ClusterShells const *, 3> shell_ptrs) {

    auto const &sh0 = shell_ptrs[0]->flattened_shells();
    auto const &sh1 = shell_ptrs[1]->flattened_shells();
    auto const &sh2 = shell_ptrs[2]->flattened_shells();

    auto const &lobound = rng.lobound();
    std::array<unsigned long, 3> lb = {{lobound[0], lobound[1], lobound[2]}};
    std::array<unsigned long, 3> ub = lb;

    auto tile = TA::TensorD(std::move(rng));

    // init map
    const double dummy = 0.0;
    auto map = TA::make_map(&dummy, {0,0,0}, {1,1,1});

    for (auto const &s0 : sh0) {
        const auto ns0 = s0.size();
        ub[0] += ns0;

        lb[1] = ub[1] = lobound[1];
        for (auto const &s1 : sh1) {
            const auto ns1 = s1.size();
            ub[1] += ns1;

            lb[2] = ub[2] = lobound[2];
            for (auto const &s2 : sh2) {
                const auto ns2 = s2.size();
                ub[2] += ns2;

                TA::remap(map, shell_set(eng, s0, s1, s2), lb, ub);
                // map.range().resize(lb, ub);
                // map.reset_data(shell_set(eng, s0, s1, s2));
                tile.block(lb, ub) = map;

                lb[2] = ub[2];
            }
            lb[1] = ub[1];
        }
        lb[0] = ub[0];
    }

    return tile;
}

// For screening
template <typename Engine>
TA::TensorD
integral_kernel(Engine &eng, TA::Range &&rng,
                std::array<tcc::basis::ClusterShells const *, 3> shell_ptrs,
                MatrixD const &X, MatrixD const &ab) {

    auto const &sh0 = shell_ptrs[0]->flattened_shells();
    auto const &sh1 = shell_ptrs[1]->flattened_shells();
    auto const &sh2 = shell_ptrs[2]->flattened_shells();

    const auto nsh0 = sh0.size();
    const auto nsh1 = sh1.size();
    const auto nsh2 = sh2.size();

    auto const &lobound = rng.lobound();
    std::array<unsigned long, 3> lb = {{lobound[0], lobound[1], lobound[2]}};
    std::array<unsigned long, 3> ub = lb;

    auto tile = TA::TensorD(std::move(rng), 0.0);

    // init map
    const double dummy = 0.0;
    auto map = TA::make_map(&dummy, {0,0,0}, {1,1,1});

    for (auto idx0 = 0ul; idx0 < nsh0; ++idx0) {
        auto const &s0 = sh0[idx0];
        ub[0] += s0.size();

        const auto X_norm_est = X(idx0);
        lb[1] = ub[1] = lobound[1];
        for (auto idx1 = 0ul; idx1 < nsh1; ++idx1) {
            auto const &s1 = sh1[idx1];
            ub[1] += s1.size();

            lb[2] = ub[2] = lobound[2];
            for (auto idx2 = 0ul; idx2 < nsh2; ++idx2) {
                auto const &s2 = sh2[idx2];
                const auto ns2 = s2.size();
                ub[2] += ns2;

                // TODO make this configurable
                if (X_norm_est * ab(idx1, idx2) > 1e-10) {
                    map.range().resize(lb, ub);
                    map.reset_data(shell_set(eng, s0, s1, s2));
                    tile.block(lb, ub) = map;
                }

                lb[2] = ub[2];
            }
            lb[1] = ub[1];
        }
        lb[0] = ub[0];
    }

    return tile;
}

template <typename Engine>
TA::TensorD
integral_kernel(Engine &eng, TA::Range &&rng,
                std::array<tcc::basis::ClusterShells const *, 4> shell_ptrs) {

    auto const &sh0 = shell_ptrs[0]->flattened_shells();
    auto const &sh1 = shell_ptrs[1]->flattened_shells();
    auto const &sh2 = shell_ptrs[2]->flattened_shells();
    auto const &sh3 = shell_ptrs[3]->flattened_shells();

    auto const &lobound = rng.lobound();
    std::array<unsigned long, 4> lb
          = {{lobound[0], lobound[1], lobound[2], lobound[3]}};
    std::array<unsigned long, 4> ub = lb;

    auto tile = TA::TensorD(std::move(rng));

    for (auto const &s0 : sh0) {
        ub[0] += s0.size();

        lb[1] = ub[1] = lobound[1];
        for (auto const &s1 : sh1) {
            ub[1] += s1.size();

            lb[2] = ub[2] = lobound[2];
            for (auto const &s2 : sh2) {
                ub[2] += s2.size();

                lb[3] = ub[3] = lobound[3];
                for (auto const &s3 : sh3) {
                    ub[3] += s3.size();

                    tile.block(lb, ub)
                          = TA::make_map(shell_set(eng, s0, s1, s2, s3), lb,
                                         ub);

                    lb[3] = ub[3];
                }
                lb[2] = ub[2];
            }
            lb[1] = ub[1];
        }
        lb[0] = ub[0];
    }

    return tile;
}

} // namespace detail
} // namespace integrals
} // namespace mpqc

#endif // MPQC_INTEGRALS_TASKINTEGRALSHELPER_H
