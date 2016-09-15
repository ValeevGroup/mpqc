#pragma once
#ifndef MPQC_TAROUTINES_MINIMIZE_STORAGE_H
#define MPQC_TAROUTINES_MINIMIZE_STORAGE_H

#include "../include/tiledarray.h"
#include "../common/typedefs.h"
#include "../common/namespaces.h"

#include "../tensor/decomposed_tensor.h"
#include "../tensor/decomposed_tensor_algebra.h"
#include "../tensor/mpqc_tile.h"

namespace mpqc {
namespace ta_routines {

inline void minimize_storage(TA::DistArray<TA::TensorD, SpPolicy> &A) {
    A.truncate();
}

inline void minimize_storage(TA::DistArray<TA::TensorD, SpPolicy> &A,
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

inline void
minimize_storage(TA::DistArray<tensor::Tile<tensor::DecomposedTensor<double>>,
                               SpPolicy> &A,
                 double clr_threshold) {
    if (clr_threshold != 0) {
        TA::foreach_inplace(
              A, [](tensor::Tile<tensor::DecomposedTensor<double>> &t_tile) {
                  auto &t = t_tile.tile();
                  auto input_norm = norm(t);

                  auto compressed_norm = input_norm;
                  if (t.cut() != 0.0) {
                      if (t.ndecomp() == 1) {
                          auto test = tensor::algebra::two_way_decomposition(t);
                          if (!test.empty()) {
                              t = test;
                              compressed_norm = norm(t);
                          }
                      } else {
                          tensor::algebra::recompress(t);
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

} // namespace ta_routines
} // namespace mpqc

#endif // MPQC_TAROUTINES_MINIMIZE_STORAGE_H
