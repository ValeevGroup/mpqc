#pragma once
#ifndef TCC_TENSOR_DECOMPOSEDTENSORGEMM_H
#define TCC_TENSOR_DECOMPOSEDTENSORGEMM_H

#include "decomposed_tensor.h"

namespace tcc {
namespace tensor {

template <typename T>
DecomposedTensor<T>
gemm(DecomposedTensor<T> const &a, DecomposedTensor<T> const &b, const T factor,
     TA::math::GemmHelper const &gemm_helper) {
    assert(false);
}

template <typename T>
DecomposedTensor<T> &gemm(DecomposedTensor<T> &c, DecomposedTensor<T> const &a,
                          DecomposedTensor<T> const &b, const T factor,
                          TA::math::GemmHelper const &gemm_helper) {
    assert(false);
}

} // namespace tensor
} // namespace tcc

#endif // TCC_TENSOR_DECOMPOSEDTENSORGEMM_H
