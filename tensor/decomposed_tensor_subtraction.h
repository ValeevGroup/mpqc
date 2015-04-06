#pragma once
#ifndef TCC_TENSOR_DECOMPOSEDTENSORSUBTRACTION_H
#define TCC_TENSOR_DECOMPOSEDTENSORSUBTRACTION_H

#include "decomposed_tensor.h"

namespace tcc {
namespace tensor {

template <typename T>
DecomposedTensor<T>
subt(DecomposedTensor<T> const &l, DecomposedTensor<T> const &r) {
    assert(false);
}

template <typename T>
DecomposedTensor<T>
subt(DecomposedTensor<T> const &l, DecomposedTensor<T> const &r,
     TA::Permutation const &p) {
    assert(false);
}

template <typename T>
DecomposedTensor<T> subt(DecomposedTensor<T> const &l,
                         DecomposedTensor<T> const &r, const T factor) {
    assert(false);
}

template <typename T>
DecomposedTensor<T>
subt(DecomposedTensor<T> const &l, DecomposedTensor<T> const &r, const T factor,
     TA::Permutation const &p) {
    assert(false);
}

template <typename T>
DecomposedTensor<T> subt(DecomposedTensor<T> const &l, const T factor) {
    assert(false);
}

template <typename T>
DecomposedTensor<T>
subt(DecomposedTensor<T> const &l, const T factor, TA::Permutation const &p) {
    assert(false);
}

template <typename T>
DecomposedTensor<T> &
subt_to(DecomposedTensor<T> &l, DecomposedTensor<T> const &r) {
    assert(false);
}

template <typename T>
DecomposedTensor<T> &
subt_to(DecomposedTensor<T> &l, DecomposedTensor<T> const &r, const T factor) {
    assert(false);
}

template <typename T>
DecomposedTensor<T> &subt_to(DecomposedTensor<T> &l, const T factor) {
    assert(false);
}

} // namespace tensor
} // namespace tcc
#endif // TCC_TENSOR_DECOMPOSEDTENSORSUBTRACTION_H
