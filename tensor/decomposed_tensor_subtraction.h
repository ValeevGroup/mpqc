#pragma once
#ifndef TCC_TENSOR_DECOMPOSEDTENSORSUBTRACTION_H
#define TCC_TENSOR_DECOMPOSEDTENSORSUBTRACTION_H

#include "decomposed_tensor.h"
#include "decomposed_tensor_unary.h"
#include "decomposed_tensor_nonintrusive_interface.h"
#include "../include/tiledarray.h"

namespace tcc {
namespace tensor {

template <typename T>
DecomposedTensor<T>
subt(DecomposedTensor<T> const &l, DecomposedTensor<T> const &r) {
    if (l.ndecomp() == 1) {
        if (r.ndecomp() == 1) {
            return DecomposedTensor<T>(l.cut(), l.tensor(0).subt(r.tensor(0)));
        }
    }
}

template <typename T>
DecomposedTensor<T>
subt(DecomposedTensor<T> const &l, TA::Tensor<double> const &r) {

    auto l_tensor = algebra::combine(l);

    // Get around Justus' range checks
    decltype(l_tensor) out_t(l_tensor.range());
    auto size = l_tensor.range().volume();
    std::transform(l_tensor.data(), l_tensor.data() + size, r.data(),
                   out_t.data(), [](T left, T right) { return left - right; });

    return DecomposedTensor<T>{l.cut(), std::move(out_t)};
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
    if (l.ndecomp() == 1) {
        if (r.ndecomp() == 1) {
            l.tensor(0).subt_to(r.tensor(0));
        }
    }
    return l;
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
