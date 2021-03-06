
#ifndef MPQC4_SRC_MPQC_MATH_TENSOR_CLR_DECOMPOSED_TENSOR_UNARY_H_
#define MPQC4_SRC_MPQC_MATH_TENSOR_CLR_DECOMPOSED_TENSOR_UNARY_H_

#include "mpqc/math/tensor/clr/decomposed_tensor.h"
#include "mpqc/math/tensor/clr/decomposed_tensor_algebra.h"

namespace mpqc {
namespace math {

/*! \brief Returns the trace of a decompose tensor
 * Currently not recommended due to implementation details
 */
template <typename T>
T trace(DecomposedTensor<T> const &t) {
  return combine(t).trace();
}

/*! \brief Returns the sum of a decomposed tensor tile
 * Currently not recommended due to implementation details
 */
template <typename T>
T sum(DecomposedTensor<T> const &t) {
  return combine(t).sum();
}

/*! \brief Returns the min of a decomposed tensor tile
 * Currently not recommended due to implementation details
 */
template <typename T>
T min(DecomposedTensor<T> const &t) {
  return combine(t).min();
}

/*! \brief Returns the max of a decomposed tensor tile
 * Currently not recommended due to implementation details
 */
template <typename T>
T max(DecomposedTensor<T> const &t) {
  return combine(t).max();
}

/*! \brief Returns the abs_min of a decomposed tensor tile
 * Currently not recommended due to implementation details
 */
template <typename T>
T abs_min(DecomposedTensor<T> const &t) {
  return combine(t).abs_min();
}

/*! \brief Returns the abs_max of a decomposed tensor tile
 * Currently not recommended due to implementation details
 */
template <typename T>
T abs_max(DecomposedTensor<T> const &t) {
  return combine(t).abs_max();
}

/*! \brief Returns the product of a decomposed tensor tile
 * Currently not recommended due to implementation details
 */
template <typename T>
T product(DecomposedTensor<T> const &t) {
  return combine(t).product();
}

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
  return combine(t).squared_norm();
}

template <typename T>
DecomposedTensor<T> permute(DecomposedTensor<T> const &t,
                            TA::Permutation const &p) {
  assert(!t.empty());
  if (t.ndecomp() == 1) {
    return DecomposedTensor<T>(t.cut(), t.tensor(0).permute(p));
  } else {
    return DecomposedTensor<T>(t.cut(), t.tensor(0), t.tensor(1).permute(p));
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

template <typename T, typename F>
DecomposedTensor<T> scale(DecomposedTensor<T> const &t, F factor) {
  auto left = t.tensor(0).scale(factor);
  if (t.ndecomp() == 2) {
    // assert(false);
    auto right = t.tensor(1).clone();
    return DecomposedTensor<T>(t.cut(), std::move(left), std::move(right));
  }

  return DecomposedTensor<T>(t.cut(), std::move(left));
}

template <typename T, typename F>
DecomposedTensor<T> scale(DecomposedTensor<T> const &t, F factor,
                          TA::Permutation const &p) {
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
DecomposedTensor<T> neg(DecomposedTensor<T> const &t,
                        TA::Permutation const &p) {
  assert(false);
}

template <typename T>
DecomposedTensor<T> &neg_to(DecomposedTensor<T> const &t,
                            TA::Permutation const &p) {
  assert(false);
}

template <typename T>
DecomposedTensor<T> &neg_to(DecomposedTensor<T> const &t) {
  assert(false);
}

template <typename T, typename F>
DecomposedTensor<T> &scale_to(DecomposedTensor<T> &t, F factor) {
  t.tensor(0).scale_to(factor);
  return t;
}

template <typename T, typename F>
DecomposedTensor<T> &scale_to(DecomposedTensor<T> &t, F factor,
                              TA::Permutation const &p) {
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

}  // namespace math
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_MATH_TENSOR_CLR_DECOMPOSED_TENSOR_UNARY_H_
