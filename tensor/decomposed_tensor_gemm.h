#pragma once
#ifndef TCC_TENSOR_DECOMPOSEDTENSORGEMM_H
#define TCC_TENSOR_DECOMPOSEDTENSORGEMM_H

#include "../include/tiledarray.h"
#include "decomposed_tensor.h"
#include "decomposed_tensor_algebra.h"
#include "decomposed_tensor_addition.h"

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
            if (a.ndecomp() == 1) {              // Reg gemm
                return DecomposedTensor<T>{a.cut(),
                                           a.tensor(0).gemm(b.tensor(0), factor,
                                                            gemm_helper)};
            } else if (a.ndecomp() == 2) { // LR gemm
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
            c.tensor(0).gemm(a.tensor(0), b.tensor(0), factor, gemm_helper);
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
            auto NoT = gemm_helper.left_op();
            auto gh = TA::math::GemmHelper(NoT, NoT, 3, 2, 3);
            auto temp = a.tensor(0)
                              .gemm(b.tensor(0), 1.0, gemm_helper)
                              .gemm(c.tensor(0), c.tensor(1), 1.0, gh);
            c = DecomposedTensor<T>(c.cut(), std::move(temp));

            return c;
        }

        auto ab = gemm(a, b, factor, gemm_helper);
        auto const &c_left_extent = c.tensor(0).range().size();
        auto const &c_right_extent = c.tensor(1).range().size();
        const auto long_dim = c_right_extent[1] * c_right_extent[2];
        const auto out_dim = c.rank() + ab.rank();

        if (out_dim >= 0.50 * std::min(c_left_extent[0], long_dim)) {
            auto NoT = gemm_helper.left_op();
            auto gh = TA::math::GemmHelper(NoT, NoT, 3, 2, 3);
            auto temp = algebra::combine(c);
            c = DecomposedTensor<T>{c.cut(), temp.gemm(ab.tensor(0),
                                                       ab.tensor(1), 1.0, gh)};
            return c;
        }

        c = add(c, ab);
        return c;
    }

    assert(false);
}

} // namespace tensor
} // namespace tcc

#endif // TCC_TENSOR_DECOMPOSEDTENSORGEMM_H
