
#ifndef MPQC_TENSOR_DECOMPOSEDTENSORSUBTRACTION_H
#define MPQC_TENSOR_DECOMPOSEDTENSORSUBTRACTION_H

#include <tiledarray.h>
#include "mpqc/math/tensor/clr/decomposed_tensor.h"
#include "mpqc/math/tensor/clr/decomposed_tensor_addition.h"
#include "mpqc/math/tensor/clr/decomposed_tensor_nonintrusive_interface.h"
#include "mpqc/math/tensor/clr/decomposed_tensor_unary.h"

namespace mpqc {
namespace tensor {

template <typename T>
DecomposedTensor<T>
subt(DecomposedTensor<T> const &l, DecomposedTensor<T> const &r) {
    if(l.ndecomp() == 1){
        if(r.ndecomp() == 1){
            return DecomposedTensor<T>(l.cut(), l.tensor(0).subt(r.tensor(0)));
        }
    }
    return add(l,tensor::scale(r,-1.0));
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

template <typename T, typename F>
DecomposedTensor<T> &
subt_to(DecomposedTensor<T> &l, DecomposedTensor<T> const &r, const F factor) {
    if (l.ndecomp() == 1) {
        if (r.ndecomp() == 1) {
            l.tensor(0).subt_to(r.tensor(0), factor);
        }
    }
    return l;
}

template <typename T>
DecomposedTensor<T> &subt_to(DecomposedTensor<T> &l, const T factor) {
    assert(false);
}

} // namespace tensor
} // namespace mpqc
#endif // MPQC_TENSOR_DECOMPOSEDTENSORSUBTRACTION_H
