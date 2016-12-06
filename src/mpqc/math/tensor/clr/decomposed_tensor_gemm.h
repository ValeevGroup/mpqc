
#ifndef MPQC4_SRC_MPQC_MATH_TENSOR_CLR_DECOMPOSED_TENSOR_GEMM_H_
#define MPQC4_SRC_MPQC_MATH_TENSOR_CLR_DECOMPOSED_TENSOR_GEMM_H_

#include <tiledarray.h>

#include "mpqc/math/tensor/clr/decomposed_tensor.h"
#include "mpqc/math/tensor/clr/decomposed_tensor_addition.h"
#include "mpqc/math/tensor/clr/decomposed_tensor_algebra.h"
#include "mpqc/math/tensor/clr/decomposed_tensor_gemm_helper.h"

namespace mpqc {
namespace tensor {

template <typename T>
TA::TensorD gemm(DecomposedTensor<T> const &a, TA::TensorD const &b,
                 const T factor, TA::math::GemmHelper const &gh) {
  TA::TensorD ac = algebra::combine(a);
  return a.gemm(b, factor, gh);
}

template <typename T>
TA::TensorD gemm(TA::TensorD &c, DecomposedTensor<T> const &a,
                 TA::TensorD const &b, const T factor,
                 TA::math::GemmHelper const &gh) {
  TA::TensorD ac = algebra::combine(a);
  return c.gemm(ac, b, factor, gh);
}

template <typename T>
DecomposedTensor<T> gemm(DecomposedTensor<T> const &a,
                         DecomposedTensor<T> const &b, const T factor,
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
        out = detail::low_rank_gemm<2, 3, 3>{}(a, b, factor, gh);
      } else if (gh.right_rank() == 1) {
        out = detail::low_rank_gemm<2, 3, 1>{}(a, b, factor, gh);
      }
    } else if (gh.left_rank() == 2) {
      if (gh.right_rank() == 2) {
        out = detail::low_rank_gemm<2, 2, 2>{}(a, b, factor, gh);
      }
    }
  } else if (gh.result_rank() == 1) {
    if (gh.left_rank() == 3) {
      if (gh.right_rank() == 2) {
        out = detail::low_rank_gemm<1, 3, 2>{}(a, b, factor, gh);
      }
    }
    if (gh.left_rank() == 2) {
      if (gh.right_rank() == 1) {
        out = detail::low_rank_gemm<1, 2, 1>{}(a, b, factor, gh);
      }
    }
  } else {
    assert(false);
    return DecomposedTensor<T>{};
  }
  assert(!out.empty());
  return out;
}

template <typename T>
DecomposedTensor<T> &gemm(DecomposedTensor<T> &c, DecomposedTensor<T> const &a,
                          DecomposedTensor<T> const &b, const T factor,
                          TA::math::GemmHelper const &gh) {
  if (gh.result_rank() == 3) {
    if (gh.left_rank() == 3) {
      if (gh.right_rank() == 2) {  // Eri3 * D
        detail::low_rank_gemm<3, 3, 2>{}(c, a, b, factor, gh);
      }
    } else if (gh.left_rank() == 2) {
      if (gh.right_rank() == 3) {
        detail::low_rank_gemm<3, 2, 3>{}(c, a, b, factor, gh);
      }
    }
  } else if (gh.result_rank() == 2) {
    if (gh.left_rank() == 3) {
      if (gh.right_rank() == 3) {
        detail::low_rank_gemm<2, 3, 3>{}(c, a, b, factor, gh);
      }
      if (gh.right_rank() == 1) {
        detail::low_rank_gemm<2, 3, 1>{}(c, a, b, factor, gh);
      }
    } else if (gh.left_rank() == 2) {
      if (gh.right_rank() == 2) {
        detail::low_rank_gemm<2, 2, 2>{}(c, a, b, factor, gh);
      }
    }
  } else if (gh.result_rank() == 1) {
    if (gh.left_rank() == 3) {
      if (gh.right_rank() == 2) {
        detail::low_rank_gemm<1, 3, 2>{}(c, a, b, factor, gh);
      }
    }
    if (gh.left_rank() == 2) {
      if (gh.right_rank() == 1) {
        detail::low_rank_gemm<1, 2, 1>{}(c, a, b, factor, gh);
      }
    }
  } else {
    assert(false);
    return c;
  }
  return c;
}

}  // namespace tensor
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_MATH_TENSOR_CLR_DECOMPOSED_TENSOR_GEMM_H_
