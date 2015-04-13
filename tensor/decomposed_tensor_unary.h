#pragma once
#ifndef TCC_TENSOR_DECOMPOSEDTENSORUNARY_H
#define TCC_TENSOR_DECOMPOSEDTENSORUNARY_H

#include "decomposed_tensor.h"
#include "decomposed_tensor_algebra.h"

namespace tcc {
namespace tensor {

// TODO try and find faster way of doing this.
template <typename T>
T norm(DecomposedTensor<T> const &t) {
    return algebra::combine(t).norm();
}

template <typename T>
DecomposedTensor<T>
permute(DecomposedTensor<T> const &t, TA::Permutation const &p) {
    assert(false);
}

template <typename T>
DecomposedTensor<T> clone(DecomposedTensor<T> const &t) {
    std::vector<TA::Tensor<T>> ts;
    ts.reserve(t.tensors().size());
    for (auto tensor : t.tensors()) {
        ts.push_back(tensor.clone());
    }
    return DecomposedTensor<T>(t.cut(), std::move(ts));
}

template <typename T>
DecomposedTensor<T> compress(DecomposedTensor<T> const &t, double cut) {
    assert(false);
}

template <typename T>
bool empty(DecomposedTensor<T> const &t) {
    return t.empty();
}

} // namespace tensor
} // namespace tcc

#endif // TCC_TENSOR_DECOMPOSEDTENSORUNARY_H
