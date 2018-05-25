//
// Created by Chong Peng on 10/6/17.
//

#ifndef SRC_MPQC_MATH_EXTERNAL_TILEDARRAY_ARRAY_MAX_N_H_
#define SRC_MPQC_MATH_EXTERNAL_TILEDARRAY_ARRAY_MAX_N_H_

#include <TiledArray/algebra/utils.h>
#include <tiledarray.h>

#include "mpqc/math/external/tiledarray/reduction.h"

namespace mpqc {

template <typename Tile, typename Policy>
std::vector<typename Tile::numeric_type> array_max_n(
    const TA::DistArray<Tile, Policy>& a, std::size_t n) {
  return a(TA::detail::dummy_annotation(a.trange().rank()))
      .reduce(detail::MaxNReduction<Tile>(n));
};

template <typename Tile, typename Policy>
std::vector<typename Tile::numeric_type> array_abs_max_n(
    const TA::DistArray<Tile, Policy>& a, std::size_t n) {
  return a(TA::detail::dummy_annotation(a.trange().rank()))
      .reduce(detail::AbsMaxNReduction<Tile>(n));
};

template <typename Tile, typename Policy>
std::vector<std::pair<typename Tile::numeric_type, std::vector<std::size_t>>>
array_max_n_index(const TA::DistArray<Tile, Policy>& a, std::size_t n) {
  return a(TA::detail::dummy_annotation(a.trange().rank()))
      .reduce(detail::MaxNIndexReduction<Tile>(n));
};

template <typename Tile, typename Policy>
std::vector<std::pair<typename Tile::numeric_type, std::vector<std::size_t>>>
array_abs_max_n_index(const TA::DistArray<Tile, Policy>& a, std::size_t n) {
  return a(TA::detail::dummy_annotation(a.trange().rank()))
      .reduce(detail::AbsMaxNIndexReduction<Tile>(n));
};

}  // namespace mpqc

#endif  // SRC_MPQC_MATH_EXTERNAL_TILEDARRAY_ARRAY_MAX_N_H_
