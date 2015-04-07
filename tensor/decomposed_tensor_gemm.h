#pragma once
#ifndef TCC_TENSOR_DECOMPOSEDTENSORGEMM_H
#define TCC_TENSOR_DECOMPOSEDTENSORGEMM_H

#include "../include/tiledarray.h"
#include "decomposed_tensor.h"
#include "decomposed_tensor_algebra.h"

namespace tcc {
namespace tensor {

template <typename T>
DecomposedTensor<T>
gemm(DecomposedTensor<T> const &a, DecomposedTensor<T> const &b, const T factor,
     TA::math::GemmHelper const &gemm_helper) {
    // Assume that b is never decomposed.
    if (gemm_helper.left_rank() == 3) {
        if (gemm_helper.right_rank() == 2
            && gemm_helper.result_rank() == 3) { // Eri3 * D
            if (a.ndecomp() == 1) {
                return DecomposedTensor<T>{a.cut(),
                                           a.tensor(0).gemm(b.tensor(0), factor,
                                                            gemm_helper)};
            } else if (a.ndecomp() == 2) {
                auto Rp = a.tensor(1).gemm(b.tensor(0), factor, gemm_helper);
                return DecomposedTensor<T>{a.cut(), a.tensor(0).clone(),
                                           std::move(Rp)};
            }
        }
    }
    assert(false);
    return DecomposedTensor<T>{};
}

template <typename T>
DecomposedTensor<T> &gemm(DecomposedTensor<T> &c, DecomposedTensor<T> const &a,
                          DecomposedTensor<T> const &b, const T factor,
                          TA::math::GemmHelper const &gemm_helper) {
    // Only density contraction for now.
    assert(gemm_helper.result_rank() == 3 && gemm_helper.left_rank() == 3
           && gemm_helper.right_rank() == 2);

    // assume b is never decomposed.
    if (c.ndecomp() == 1) {
        if (a.ndecomp() == 1) {
            c.tensor(0).gemm(a, b, factor, gemm_helper);
            return c;
        }

        auto Rp = a.tensor(1).gemm(b.tensor(0), factor, gemm_helper);
        auto NoT = gemm_helper.left_op();
        auto gh = TA::math::GemmHelper(NoT, NoT, c.tensor(0).range().dim(),
                                        a.tensor(0).range().dim(),
                                        Rp.range().dim());
        c.tensor(0).gemm(a.tensor(0), Rp, factor, gh);
        return c;
    } else {
        if (a.ndecomp() == 1) {
            auto ab = gemm(a, b, factor, gemm_helper);
            c = algebra::combine(c).add(ab);
            return c;
        }

        auto ab = gemm(a, b, factor, gemm_helper);

        auto const &c_left_extent = c.tensor(0).range().size();
        const auto out_dim = c.rank() + ab.rank();
        TA::Range l_range(c_left_extent[0], out_dim);
        TA::Tensor<T> l_tensor(std::move(l_range));

        auto const &c_right_extent = c.tensor(1).range().size();
        TA::Range r_range(out_dim, c_right_extent[1], c_right_extent[2]);
        TA::Tensor<T> r_tensor(std::move(r_range));

        // Fill L
        auto Lmap = TA::eigen_map(l_tensor, c_left_extent[0], out_dim);
        auto c_lmap = TA::eigen_map(c.tensor(0), c_left_extent[0], c.rank());
        auto ab_lmap = TA::eigen_map(ab.tensor(0), c_left_extent[0], ab.rank());
        Lmap.leftCols(c.rank()) = c_lmap;
        Lmap.rightCols(ab.rank()) = ab_lmap;

        // Fill R
        const auto long_dim = c_right_extent[1] * c_right_extent[2];
        auto Rmap = TA::eigen_map(r_tensor, out_dim, long_dim);
        auto c_rmap = TA::eigen_map(c.tensor(1), c.rank(), long_dim);
        auto ab_rmap = TA::eigen_map(ab.tensor(1), ab.rank(), long_dim);

        Rmap.topRows(c.rank()) = c_rmap;
        Rmap.bottomRows(ab.rank()) = ab_rmap;
        c = DecomposedTensor<T>{c.cut(), std::move(l_tensor),
                                std::move(r_tensor)};
        return c;
    }

    assert(false);
}

} // namespace tensor
} // namespace tcc

#endif // TCC_TENSOR_DECOMPOSEDTENSORGEMM_H
