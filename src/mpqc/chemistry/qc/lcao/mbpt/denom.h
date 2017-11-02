//
// Created by Chong Peng on 10/1/15.
//

#ifndef SRC_MPQC_CHEMISTRY_QC_LCAO_MBPT_DENOM_H_
#define SRC_MPQC_CHEMISTRY_QC_LCAO_MBPT_DENOM_H_

#include <tiledarray.h>

#include "mpqc/math/external/eigen/eigen.h"

namespace mpqc {
namespace lcao {

namespace detail {

// reduce matrix 1/(ei + ej - ea - eb)
template <typename Tile, typename Policy>
void d_abij_inplace(TA::Array<double, 4, Tile, Policy> &abij,
                    const EigenVector<typename Tile::numeric_type> &ens,
                    std::size_t n_occ, std::size_t n_frozen,
                    typename Tile::numeric_type shift = 0.0) {
  auto convert = [&ens, n_occ, n_frozen, shift](Tile &result_tile) {

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
    typename Tile::numeric_type norm = 0.0;
    for (auto a = a0; a < an; ++a) {
      const auto e_a = shift - ens[a + n_occ];
      for (auto b = b0; b < bn; ++b) {
        const auto e_ab = e_a - ens[b + n_occ];
        for (auto i = i0; i < in; ++i) {
          const auto e_abi = e_ab + ens[i + n_frozen];
          for (auto j = j0; j < jn; ++j, ++tile_idx) {
            const auto e_abij = e_abi + ens[j + n_frozen];
            const auto old = result_tile[tile_idx];
            const auto result_abij = old / (e_abij);
            norm += result_abij * result_abij;
            result_tile[tile_idx] = result_abij;
          }
        }
      }
    }
    return std::sqrt(norm);
  };

  TA::foreach_inplace(abij, convert);
  abij.world().gop.fence();
}

template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> d_abij(
    const TA::DistArray<Tile, Policy> &abij,
    const EigenVector<typename Tile::numeric_type> &ens, std::size_t n_occ,
    std::size_t n_frozen, typename Tile::numeric_type shift = 0.0) {
  auto convert = [&ens, n_occ, n_frozen, shift](Tile &result_tile,
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
    typename Tile::numeric_type norm = 0.0;
    for (auto a = a0; a < an; ++a) {
      const auto e_a = shift - ens[a + n_occ];
      for (auto b = b0; b < bn; ++b) {
        const auto e_ab = e_a - ens[b + n_occ];
        for (auto i = i0; i < in; ++i) {
          const auto e_abi = e_ab + ens[i + n_frozen];
          for (auto j = j0; j < jn; ++j, ++tile_idx) {
            const auto e_abij = e_abi + ens[j + n_frozen];
            const auto old = arg_tile[tile_idx];
            const auto result_abij = old / (e_abij);
            norm += result_abij * result_abij;
            result_tile[tile_idx] = result_abij;
          }
        }
      }
    }
    return std::sqrt(norm);
  };

  auto result = TA::foreach (abij, convert);
  abij.world().gop.fence();
  return result;
}

// create matrix d("a,i") = 1/(ei - ea)
template <typename Tile, typename Policy>
TA::DistArray<
    Tile, typename std::enable_if<std::is_same<Policy, TA::SparsePolicy>::value,
                                  TA::SparsePolicy>::type>
create_d_ai(madness::World &world, const TA::TiledRange &trange,
            const EigenVector<typename Tile::numeric_type> ens,
            std::size_t n_occ, std::size_t n_frozen) {
  typedef typename TA::DistArray<Tile, Policy>::range_type range_type;

  auto make_tile = [&ens, n_occ, n_frozen](range_type &range, std::size_t ord,
                                           Tile *out_tile, TA::TensorF *norms) {

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
    const auto tile_volume = result_tile.range().volume();
    const auto tile_norm = result_tile.norm();
    bool save_norm =
        tile_norm >= tile_volume * TA::SparseShape<float>::threshold();
    if (save_norm) {
      *out_tile = result_tile;
      (*norms)[ord] = tile_norm;
    }
  };

  const auto tvolume = trange.tiles_range().volume();
  std::vector<Tile> tiles(tvolume);
  TA::TensorF tile_norms(trange.tiles_range(), 0.0);
  auto pmap = TA::SparsePolicy::default_pmap(world, tvolume);

  for (auto const ord : *pmap) {
    world.taskq.add(make_tile, trange.make_tile_range(ord), ord, &tiles[ord],
                    &tile_norms);
  }

  world.gop.fence();

  TA::SparseShape<float> shape(world, tile_norms, trange);
  TA::DistArray<Tile, Policy> result(world, trange, shape, pmap);

  for (auto const ord : *pmap) {
    if (result.is_local(ord) && !result.is_zero(ord)) {
      auto &tile = tiles[ord];
      assert(!tile.empty());
      result.set(ord, tile);
    }
  }

  world.gop.fence();
  result.truncate();
  return result;
}

// create matrix d("a,i") = 1/(ei - ea)
template <typename Tile, typename Policy>
TA::DistArray<
    Tile, typename std::enable_if<std::is_same<Policy, TA::DensePolicy>::value,
                                  TA::DensePolicy>::type>
create_d_ai(madness::World &world, const TA::TiledRange &trange,
            const Eigen::VectorXd &ens, std::size_t n_occ,
            std::size_t n_frozen) {
  typedef typename TA::DistArray<Tile, Policy>::range_type range_type;

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

  TA::DistArray<Tile, Policy> result(world, trange);

  for (auto it = result.begin(); it != result.end(); ++it) {
    madness::Future<Tile> tile =
        world.taskq.add(make_tile, trange.make_tile_range(it.ordinal()));

    *it = tile;
  }

  world.gop.fence();
  return result;
}

}  // namespace detail
}  // namespace lcao
}  // namespace mpqc

#endif  // SRC_MPQC_CHEMISTRY_QC_LCAO_MBPT_DENOM_H_
