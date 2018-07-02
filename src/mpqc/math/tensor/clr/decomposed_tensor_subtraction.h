
#ifndef MPQC4_SRC_MPQC_MATH_TENSOR_CLR_DECOMPOSED_TENSOR_SUBTRACTION_H_
#define MPQC4_SRC_MPQC_MATH_TENSOR_CLR_DECOMPOSED_TENSOR_SUBTRACTION_H_

#include "mpqc/math/tensor/clr/decomposed_tensor.h"
#include "mpqc/math/tensor/clr/decomposed_tensor_addition.h"
#include "mpqc/math/tensor/clr/decomposed_tensor_nonintrusive_interface.h"
#include "mpqc/math/tensor/clr/decomposed_tensor_unary.h"
#include <tiledarray.h>

namespace mpqc {
namespace math {

template <typename T>
DecomposedTensor<T> subt(DecomposedTensor<T> const &l,
                         DecomposedTensor<T> const &r) {
  if (l.ndecomp() == 1) {
    if (r.ndecomp() == 1) {
      return DecomposedTensor<T>(l.cut(), l.tensor(0).subt(r.tensor(0)));
    }
  }
  return add(l, scale(r, -1.0));
}

template <typename T>
DecomposedTensor<T> subt(DecomposedTensor<T> const &l,
                         DecomposedTensor<T> const &r,
                         TA::Permutation const &p) {
  assert(false);
}

template <typename T>
DecomposedTensor<T> subt(DecomposedTensor<T> const &l,
                         DecomposedTensor<T> const &r, const T factor) {
  assert(false);
}

template <typename T>
DecomposedTensor<T> subt(DecomposedTensor<T> const &l,
                         DecomposedTensor<T> const &r, const T factor,
                         TA::Permutation const &p) {
  assert(false);
}

template <typename T>
DecomposedTensor<T> subt(DecomposedTensor<T> const &l, const T factor) {
  assert(false);
}

template <typename T>
DecomposedTensor<T> subt(DecomposedTensor<T> const &l, const T factor,
                         TA::Permutation const &p) {
  assert(false);
}

template <typename T>
DecomposedTensor<T> &subt_to(DecomposedTensor<T> &l,
                             DecomposedTensor<T> const &r) {
  if (l.ndecomp() == 1) {
    if (r.ndecomp() == 1) {
      l.tensor(0).subt_to(r.tensor(0));
    }
  }
  return l;
}

template <typename T, typename F>
DecomposedTensor<T> &subt_to(DecomposedTensor<T> &l,
                             DecomposedTensor<T> const &r, const F factor) {
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

}  // namespace math
}  // namespace mpqc
#endif  // MPQC4_SRC_MPQC_MATH_TENSOR_CLR_DECOMPOSED_TENSOR_SUBTRACTION_H_
