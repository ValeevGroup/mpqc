#pragma once
#ifndef TCC_TENSOR_DECOMPOSEDTENSORUNARY_H
#define TCC_TENSOR_DECOMPOSEDTENSORUNARY_H

#include "decomposed_tensor.h"

namespace tcc {
namespace tensor {

template <typename T>
T norm(DecomposedTensor<T> const &t) {
    assert(false);
}

template <typename T>
DecomposedTensor<T>
permute(DecomposedTensor<T> const &t, TA::Permutation const &p) {
    assert(false);
}

template <typename T>
DecomposedTensor<T> clone(DecomposedTensor<T> const &t) {
    assert(false);
}

template <typename T>
DecomposedTensor<T> compress(DecomposedTensor<T> const &t, double cut) {
    assert(false);
}

} // namespace tensor
} // namespace tcc

#endif // TCC_TENSOR_DECOMPOSEDTENSORUNARY_H
