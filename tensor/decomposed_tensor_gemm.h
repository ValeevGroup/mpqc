#pragma once
#ifndef TCC_TENSOR_DECOMPOSEDTENSORGEMM_H
#define TCC_TENSOR_DECOMPOSEDTENSORGEMM_H

#include "../include/tiledarray.h"
#include "decomposed_tensor.h"
#include "decomposed_tensor_algebra.h"
#include "decomposed_tensor_addition.h"
#include "decomposed_tensor_gemm_helper.h"

namespace tcc {
namespace tensor {

template <typename T>
DecomposedTensor<T>
gemm(DecomposedTensor<T> const &a, DecomposedTensor<T> const &b, const T factor,
     TA::math::GemmHelper const &gh) {

    auto out = DecomposedTensor<T>{};
    if (gh.result_rank() == 3) {
        if (gh.left_rank() == 3) {
            if (gh.right_rank() == 2) {
                out = detail::low_rank_gemm<3, 3, 2>{}(a, b, factor, gh);
            }
        } else if (gh.left_rank() == 2) {
            if (gh.right_rank() == 3) {
                out = detail::low_rank_gemm<3, 2, 3>{}(a, b, factor, gh);
            }
        }
    } else if (gh.result_rank() == 2) {
        if (gh.left_rank() == 3) {
            if (gh.right_rank() == 3) {
                out = detail::low_rank_gemm<2,3,3>{}(a,b,factor,gh);
            }
        }
    } else {
        assert(false);
        return DecomposedTensor<T>{};
    }
    return out;
}

template <typename T>
DecomposedTensor<T> &gemm(DecomposedTensor<T> &c, DecomposedTensor<T> const &a,
                          DecomposedTensor<T> const &b, const T factor,
                          TA::math::GemmHelper const &gh) {
    if (gh.result_rank() == 3) {
        if (gh.left_rank() == 3) {
            if (gh.right_rank() == 2) { // Eri3 * D
                detail::low_rank_gemm<3, 3, 2>{}(c, a, b, factor, gh);
            }
        } else if (gh.left_rank() == 2) {
            if (gh.right_rank() == 3) {
                detail::low_rank_gemm<3, 2, 3>{}(c, a, b, factor, gh);
            }
        }
    } else if(gh.result_rank() == 2){
        if(gh.left_rank() == 3){
            if(gh.right_rank() == 3){
                detail::low_rank_gemm<2,3,3>{}(c, a, b, factor, gh);
            }
        }
    }else {
        assert(false);
        return c;
    }
    return c;
}

} // namespace tensor
} // namespace tcc

#endif // TCC_TENSOR_DECOMPOSEDTENSORGEMM_H
