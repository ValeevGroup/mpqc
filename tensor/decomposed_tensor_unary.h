#pragma once
#ifndef TCC_TENSOR_DECOMPOSEDTENSORUNARY_H
#define TCC_TENSOR_DECOMPOSEDTENSORUNARY_H

#include "decomposed_tensor.h"
#include "decomposed_tensor_algebra.h"

namespace tcc {
namespace tensor {


template <typename T>
T norm(DecomposedTensor<T> const &t) {
    if (t.empty()) {
        return T(0.0);
    }
    // Norm is bounded by ||M|| <= ||S^M|| * ||T^M||
    // if single tile then norm is exact
    auto norm_bound = 1.0;
    for (auto const &tensor : t.tensors()) {
        norm_bound *= tensor.norm();
    }
    return norm_bound;
}
template <typename T>
T squared_norm(DecomposedTensor<T> const &t) {
    return algebra::combine(t).squared_norm();
}

template <typename T>
DecomposedTensor<T>
permute(DecomposedTensor<T> const &t, TA::Permutation const &p) {
    assert(!t.empty());
    if (t.ndecomp() == 1) {
        return DecomposedTensor<T>(t.cut(), t.tensor(0).permute(p));
    } else {
        return DecomposedTensor<T>(t.cut(), t.tensor(0),
                                   t.tensor(1).permute(p));
    }
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
DecomposedTensor<T> scale(DecomposedTensor<T> const &t, T factor) {
    auto left = t.tensor(0).scale(factor);
    if (t.ndecomp() == 2) {
        assert(false);
        //        right = t.tensor(1).clone();
        //       return DecomposedTensor<T>(t.cut(), std::move(left),
        //       std::move(right));
    }

    return DecomposedTensor<T>(t.cut(), std::move(left));
}

template <typename T>
DecomposedTensor<T>
scale(DecomposedTensor<T> const &t, T factor, TA::Permutation const &p) {
    assert(false);
}

template <typename T>
DecomposedTensor<T> neg(DecomposedTensor<T> const &t) {
    auto left = t.tensor(0).neg();
    auto right = TA::Tensor<T>();
    if (t.ndecomp() == 2) {
        right = t.tensor(1).clone();
        return DecomposedTensor<T>(t.cut(), std::move(left), std::move(right));
    }

    return DecomposedTensor<T>(t.cut(), std::move(left));
}

template <typename T>
DecomposedTensor<T>
neg(DecomposedTensor<T> const &t, TA::Permutation const &p) {
    assert(false);
}

template <typename T>
DecomposedTensor<T> &
neg_to(DecomposedTensor<T> const &t, TA::Permutation const &p) {
    assert(false);
}

template <typename T>
DecomposedTensor<T> &neg_to(DecomposedTensor<T> const &t) {
    assert(false);
}

template <typename T>
DecomposedTensor<T> &scale_to(DecomposedTensor<T> &t, T factor) {
    t.tensor(0).scale_to(factor);
    return t;
}

template <typename T>
DecomposedTensor<T> &
scale_to(DecomposedTensor<T> &t, T factor, TA::Permutation const &p) {
    assert(false);
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
