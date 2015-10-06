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
#include "../basis/basis_fwd.h"

namespace mpqc {
namespace integrals {
namespace detail {

using OneE_Engine = libint2::OneBodyEngine;
using TwoE_Engine = libint2::TwoBodyEngine<libint2::Coulomb>;
using ShellVec = basis::ShellVec;

inline const double *shell_set(TwoE_Engine &e, Shell const &s0, Shell const &s1,
                               Shell const &s2, Shell const &s3) {
    return e.compute(s0, s1, s2, s3);
}

inline const double *
shell_set(TwoE_Engine &e, Shell const &s0, Shell const &s1) {
    const auto unit = Shell::unit();
    return e.compute(s0, unit, s1, unit);
}

inline const double *
shell_set(TwoE_Engine &e, Shell const &s0, Shell const &s1, Shell const &s2) {
    const auto unit = Shell::unit();
    return e.compute(s0, unit, s1, s2);
}

inline const double *
shell_set(OneE_Engine &e, Shell const &s0, Shell const &s1) {
    return e.compute(s0, s1);
}

template <typename Engine>
TA::TensorD
integral_kernel(Engine &eng, TA::Range &&rng,
                std::array<ShellVec const *, 2> shell_ptrs) {

    auto const &sh0 = *shell_ptrs[0];
    auto const &sh1 = *shell_ptrs[1];

    auto const &lowbound = rng.lobound();
    const auto start0 = lowbound[0];
    const auto start1 = lowbound[1];

    auto tile = TA::TensorD(std::move(rng));

    auto bf0 = start0;
    for (auto const &s0 : sh0) {
        const auto ns0 = s0.size();

        auto bf1 = start1;
        for (auto const &s1 : sh1) {
            const auto ns1 = s1.size();

            const auto lb = {bf0, bf1};
            const auto ub = {bf0 + ns0, bf1 + ns1};
            tile.block(lb, ub) = TA::make_map(shell_set(eng, s0, s1), lb, ub);

            bf1 += ns1;
        }
        bf0 += ns0;
    }

    return tile;
}

template <typename Engine>
TA::TensorD
integral_kernel(Engine &eng, TA::Range &&rng,
                std::array<ShellVec const *, 3> shell_ptrs) {

    auto const &sh0 = *shell_ptrs[0];
    auto const &sh1 = *shell_ptrs[1];
    auto const &sh2 = *shell_ptrs[2];

    auto const &lowbound = rng.lobound();
    const auto start0 = lowbound[0];
    const auto start1 = lowbound[1];
    const auto start2 = lowbound[2];

    auto tile = TA::TensorD(std::move(rng));

    auto bf0 = start0;
    for (auto const &s0 : sh0) {
        const auto ns0 = s0.size();

        auto bf1 = start1;
        for (auto const &s1 : sh1) {
            const auto ns1 = s1.size();

            auto bf2 = start2;
            for (auto const &s2 : sh2) {
                const auto ns2 = s2.size();

                const auto lb = {bf0, bf1, bf2};
                const auto ub = {bf0 + ns0, bf1 + ns1, bf2 + ns2};
                tile.block(lb, ub)
                      = TA::make_map(shell_set(eng, s0, s1, s2), lb, ub);

                bf2 += ns2;
            }
            bf1 += ns1;
        }
        bf0 += ns0;
    }

    return tile;
}

// For screening
template <typename Engine>
TA::TensorD
integral_kernel(Engine &eng, TA::Range &&rng,
                std::array<ShellVec const *, 3> shell_ptrs,
                MatrixD const &X, MatrixD const &ab) {

    auto const &sh0 = *shell_ptrs[0];
    auto const &sh1 = *shell_ptrs[1];
    auto const &sh2 = *shell_ptrs[2];

    const auto nsh0 = sh0.size();
    const auto nsh1 = sh1.size();
    const auto nsh2 = sh2.size();

    auto const &lowbound = rng.lobound();
    const auto start0 = lowbound[0];
    const auto start1 = lowbound[1];
    const auto start2 = lowbound[2];

    // Since screening it's important to zero this.
    auto tile = TA::TensorD(std::move(rng), 0.0);

    auto bf0 = start0;
    for (auto idx0 = 0ul; idx0 < nsh0; ++idx0) {
        auto const &s0 = sh0[idx0];
        const auto ns0 = s0.size();
        const auto X_norm_est = X(idx0);

        auto bf1 = start1;
        for (auto idx1 = 0ul; idx1 < nsh1; ++idx1) {
            auto const &s1 = sh1[idx1];
            const auto ns1 = s1.size();

            auto bf2 = start2;
            for (auto idx2 = 0ul; idx2 < nsh2; ++idx2) {
                auto const &s2 = sh2[idx2];
                const auto ns2 = s2.size();

                // Screen that bad boy
                if (X_norm_est * ab(idx1, idx2) > 1e-10) {
                    const auto lb = {bf0, bf1, bf2};
                    const auto ub = {bf0 + ns0, bf1 + ns1, bf2 + ns2};
                    tile.block(lb, ub)
                          = TA::make_map(shell_set(eng, s0, s1, s2), lb, ub);
                }

                bf2 += ns2;
            }
            bf1 += ns1;
        }
        bf0 += ns0;
    }

    return tile;
}

template <typename Engine>
TA::TensorD
integral_kernel(Engine &eng, TA::Range &&rng,
                std::array<ShellVec const *, 4> shell_ptrs) {

    auto const &sh0 = *shell_ptrs[0];
    auto const &sh1 = *shell_ptrs[1];
    auto const &sh2 = *shell_ptrs[2];
    auto const &sh3 = *shell_ptrs[3];

    auto const &lowbound = rng.lobound();
    const auto start0 = lowbound[0];
    const auto start1 = lowbound[1];
    const auto start2 = lowbound[2];
    const auto start3 = lowbound[3];

    auto tile = TA::TensorD(std::move(rng));

    auto bf0 = start0;
    for (auto const &s0 : sh0) {
        const auto ns0 = s0.size();

        auto bf1 = start1;
        for (auto const &s1 : sh1) {
            const auto ns1 = s1.size();

            auto bf2 = start2;
            for (auto const &s2 : sh2) {
                const auto ns2 = s2.size();

                auto bf3 = start3;
                for (auto const &s3 : sh3) {
                    const auto ns3 = s3.size();

                    const auto lb = {bf0, bf1, bf2, bf3};
                    const auto ub
                          = {bf0 + ns0, bf1 + ns1, bf2 + ns2, bf3 + ns3};

                    tile.block(lb, ub)
                          = TA::make_map(shell_set(eng, s0, s1, s2, s3), lb,
                                         ub);

                    bf3 += ns3;
                }
                bf2 += ns2;
            }
            bf1 += ns1;
        }
        bf0 += ns0;
    }

    return tile;
}


} // namespace detail
} // namespace integrals
} // namespace mpqc

#endif // MPQC_INTEGRALS_TASKINTEGRALSHELPER_H
