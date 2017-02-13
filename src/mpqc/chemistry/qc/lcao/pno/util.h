/*
 * util.h
 *
 *  Created on: Feb 9, 2017
 *      Author: jinmeizhang
 */

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_LCAO_PNO_UTIL_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_LCAO_PNO_UTIL_H_

#include <tiledarray.h>

#include "mpqc/math/external/eigen/eigen.h"

namespace mpqc {

// reduce matrix 1/(ei + ej - Faa - Fbb)
// where a b are PNO indices
template <typename Tile, typename Policy>
void d_ab_inplace(TA::DistArray<Tile, Policy>& abij,
                  const Eigen::VectorXd& faa,
                  const double e_i, const double e_j) {

  auto convert = [&](Tile &result_tile) {

    // compute index
    const auto a0 = result_tile.range().lobound()[0];
    const auto an = result_tile.range().upbound()[0];
    const auto b0 = result_tile.range().lobound()[1];
    const auto bn = result_tile.range().upbound()[1];

    auto tile_idx = 0;
    typename Tile::value_type norm = 0.0;
    for (auto a = a0; a < an; ++a) {
      const auto e_a = faa[a];
      for (auto b = b0; b < bn; ++b, ++tile_idx) {
        const auto e_b = faa[b];
        const auto e_iajb = e_i + e_j - e_a - e_b;
        const auto old = result_tile[tile_idx];
        const auto result_abij = old / (e_iajb);
        norm += result_abij * result_abij;
        result_tile[tile_idx] = result_abij;
      }
    }
    return std::sqrt(norm);
  };

  TA::foreach_inplace(abij, convert);
  abij.world().gop.fence();
}

// compute the norm of a vector of DistArray
template <typename Tile, typename Policy>
double norm_vec_tensor(const std::vector<TA::DistArray<Tile, Policy>>& vec_array) {
  double norm_sq = 0.0;
  for (const auto& array: vec_array) {
    const double norm = TA::norm2(array);
    norm_sq = norm * norm;
  }
  return std::sqrt(norm_sq);
}
} // end of namespace mpac


#endif /* SRC_MPQC_CHEMISTRY_QC_LCAO_PNO_UTIL_H_ */
