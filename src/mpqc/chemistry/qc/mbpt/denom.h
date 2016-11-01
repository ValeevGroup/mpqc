//
// Created by Chong Peng on 10/1/15.
//

#ifndef SRC_MPQC_CHEMISTRY_QC_MBPT_DENOM_H
#define SRC_MPQC_CHEMISTRY_QC_MBPT_DENOM_H

#include <tiledarray.h>
#include "mpqc/math/external/eigen/eigen.h"

namespace mpqc {

// reduce matrix 1/(ei + ej - ea - eb)
template <typename Tile, typename Policy>
void d_abij_inplace(TA::Array<double, 4, Tile, Policy> &abij,
                    const Eigen::VectorXd &ens, std::size_t n_occ,
                    std::size_t n_frozen) {
  auto convert = [&ens, n_occ, n_frozen](Tile &result_tile) {

    // compute index
    const auto a0 = result_tile.range().lobound()[0];
    const auto an = result_tile.range().upbound()[0];
    const auto b0 = result_tile.range().lobound()[1];
    const auto bn = result_tile.range().upbound()[1];
    const auto i0 = result_tile.range().lobound()[2];
    const auto in = result_tile.range().upbound()[2];
    const auto j0 = result_tile.range().lobound()[3];
    const auto jn = result_tile.range().upbound()[3];

    auto tile_idx = 0;
    typename Tile::value_type norm = 0.0;
    for (auto a = a0; a < an; ++a) {
      const auto e_a = ens[a + n_occ];
      for (auto b = b0; b < bn; ++b) {
        const auto e_b = ens[b + n_occ];
        for (auto i = i0; i < in; ++i) {
          const auto e_i = ens[i + n_frozen];
          for (auto j = j0; j < jn; ++j, ++tile_idx) {
            const auto e_j = ens[j + n_frozen];
            const auto e_iajb = e_i + e_j - e_a - e_b;
            const auto old = result_tile[tile_idx];
            const auto result_abij = old / (e_iajb);
            norm += result_abij * result_abij;
            result_tile[tile_idx] = result_abij;
          }
        }
      }
    }
    return std::sqrt(norm);
  };

  TA::foreach_inplace(abij, convert);
}

template <typename Tile, typename Policy>
TA::Array<double, 4, Tile, Policy> d_abij(
    TA::Array<double, 4, Tile, Policy> &abij, const Eigen::VectorXd &ens,
    std::size_t n_occ, std::size_t n_frozen) {
  auto convert = [&ens, n_occ, n_frozen](Tile &result_tile,
                                         const Tile &arg_tile) {

    result_tile = Tile(arg_tile.range());

    // compute index
    const auto a0 = result_tile.range().lobound()[0];
    const auto an = result_tile.range().upbound()[0];
    const auto b0 = result_tile.range().lobound()[1];
    const auto bn = result_tile.range().upbound()[1];
    const auto i0 = result_tile.range().lobound()[2];
    const auto in = result_tile.range().upbound()[2];
    const auto j0 = result_tile.range().lobound()[3];
    const auto jn = result_tile.range().upbound()[3];

    auto tile_idx = 0;
    typename Tile::value_type norm = 0.0;
    for (auto a = a0; a < an; ++a) {
      const auto e_a = ens[a + n_occ];
      for (auto b = b0; b < bn; ++b) {
        const auto e_b = ens[b + n_occ];
        for (auto i = i0; i < in; ++i) {
          const auto e_i = ens[i + n_frozen];
          for (auto j = j0; j < jn; ++j, ++tile_idx) {
            const auto e_j = ens[j + n_frozen];
            const auto e_iajb = e_i + e_j - e_a - e_b;
            const auto old = arg_tile[tile_idx];
            const auto result_abij = old / (e_iajb);
            norm += result_abij * result_abij;
            result_tile[tile_idx] = result_abij;
          }
        }
      }
    }
    return std::sqrt(norm);
  };

  return TA::foreach (abij, convert);
}

// create matrix d("a,b,i,j) = 1/(ei + ej - ea - eb)
template <typename Tile, typename Policy>
void create_d_abij(TA::Array<double, 4, Tile, Policy> &abij,
                   const Eigen::VectorXd &ens, std::size_t n_occ,
                   std::size_t n_frozen) {
  typedef typename TA::Array<double, 4, Tile, Policy>::range_type range_type;
  typedef typename TA::Array<double, 4, Tile, Policy>::iterator iterator;

  auto make_tile = [&ens, n_occ, n_frozen](range_type &range) {

    auto result_tile = Tile(range);

    // compute index
    const auto a0 = result_tile.range().lobound()[0];
    const auto an = result_tile.range().upbound()[0];
    const auto b0 = result_tile.range().lobound()[1];
    const auto bn = result_tile.range().upbound()[1];
    const auto i0 = result_tile.range().lobound()[2];
    const auto in = result_tile.range().upbound()[2];
    const auto j0 = result_tile.range().lobound()[3];
    const auto jn = result_tile.range().upbound()[3];

    auto tile_idx = 0;
    typename Tile::value_type tmp = 1.0;
    for (auto a = a0; a < an; ++a) {
      const auto e_a = ens[a + n_occ];
      for (auto b = b0; b < bn; ++b) {
        const auto e_b = ens[b + n_occ];
        for (auto i = i0; i < in; ++i) {
          const auto e_i = ens[i + n_frozen];
          for (auto j = j0; j < jn; ++j, ++tile_idx) {
            const auto e_j = ens[j + n_frozen];
            const auto e_iajb = e_i + e_j - e_a - e_b;
            const auto result_abij = tmp / (e_iajb);
            result_tile[tile_idx] = result_abij;
          }
        }
      }
    }
    return result_tile;
  };

  for (iterator it = abij.begin(); it != abij.end(); ++it) {
    madness::Future<Tile> tile = abij.world().taskq.add(
        make_tile, abij.trange().make_tile_range(it.ordinal()));

    *it = tile;
  }
  abij.world().gop.fence();
}

// create matrix d("a,i") = 1/(ei - ea)
template <typename Tile, typename Policy>
void create_d_ai(TA::Array<double, 2, Tile, Policy> &f_ai,
                 const Eigen::VectorXd &ens, std::size_t n_occ,
                 std::size_t n_frozen) {
  typedef typename TA::Array<double, 2, Tile, Policy>::range_type range_type;
  typedef typename TA::Array<double, 2, Tile, Policy>::iterator iterator;

  auto make_tile = [&ens, n_occ, n_frozen](range_type &range) {

    auto result_tile = Tile(range);
    const auto a0 = result_tile.range().lobound()[0];
    const auto an = result_tile.range().upbound()[0];
    const auto i0 = result_tile.range().lobound()[1];
    const auto in = result_tile.range().upbound()[1];

    auto ai = 0;
    typename Tile::value_type tmp = 1.0;
    for (auto a = a0; a < an; ++a) {
      const auto e_a = ens[a + n_occ];
      for (auto i = i0; i < in; ++i, ++ai) {
        const auto e_i = ens[i + n_frozen];
        const auto e_ia = e_i - e_a;
        const auto result_ai = tmp / (e_ia);
        result_tile[ai] = result_ai;
      }
    }
    return result_tile;
  };

  for (iterator it = f_ai.begin(); it != f_ai.end(); ++it) {
    madness::Future<Tile> tile = f_ai.world().taskq.add(
        make_tile, f_ai.trange().make_tile_range(it.ordinal()));

    *it = tile;
  }
  f_ai.world().gop.fence();
}
}

#endif  // SRC_MPQC_CHEMISTRY_QC_MBPT_DENOM_H
