#ifndef MPQC4_SRC_MPQC_MATH_TENSOR_CLR_MINIMIZE_STORAGE_H_
#define MPQC4_SRC_MPQC_MATH_TENSOR_CLR_MINIMIZE_STORAGE_H_

#include <tiledarray.h>

#include "mpqc/math/tensor/clr/decomposed_tensor.h"
#include "mpqc/math/tensor/clr/decomposed_tensor_algebra.h"
#include "tile.h"

namespace mpqc {

template <typename Policy>
inline void minimize_storage(TA::DistArray<TA::TensorD, Policy> &A) {
  A.truncate();
}

inline void minimize_storage(TA::DistArray<TA::TensorD, TA::DensePolicy> &A,
                             double truncate_threshold)
{
}

inline void minimize_storage(TA::DistArray<TA::TensorD, TA::SparsePolicy> &A,
                             double truncate_threshold) {
  if (truncate_threshold > TA::SparseShape<float>::threshold()) {
    TA::foreach_inplace(A, [=](TA::TensorD &t) {
      const auto norm = t.norm();
      const auto volume = t.range().volume();
      const auto val = norm / double(volume);
      return (val < truncate_threshold) ? 0.0 : norm;
    });
  } else {
    A.truncate();
  }
}

inline void minimize_storage(
    TA::DistArray<math::Tile<math::DecomposedTensor<double>>,
                  TA::SparsePolicy> &A,
    double clr_threshold) {
  if (clr_threshold != 0) {
    TA::foreach_inplace(
        A, [](math::Tile<math::DecomposedTensor<double>> &t_tile) {
          auto &t = t_tile.tile();
          auto input_norm = norm(t);

          auto compressed_norm = input_norm;
          if (t.cut() != 0.0) {
            if (t.ndecomp() == 1) {
              auto test = two_way_decomposition(t);
              if (!test.empty()) {
                t = test;
                compressed_norm = norm(t);
              }
            } else {
              recompress(t);
              compressed_norm = norm(t);
            }
          }

          // Both are always larger than or equal to the
          // real
          // norm.
          return std::min(input_norm, compressed_norm);
        });
  } else {
    A.truncate();
  }
}

}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_MATH_TENSOR_CLR_MINIMIZE_STORAGE_H_
