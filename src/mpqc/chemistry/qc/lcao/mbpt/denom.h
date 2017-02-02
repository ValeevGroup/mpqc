//
// Created by Chong Peng on 10/1/15.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_MBPT_DENOM_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_MBPT_DENOM_H_

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
  abij.world().gop.fence();
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

  auto result = TA::foreach (abij, convert);
  abij.world().gop.fence();
  return result;
}

template <typename Tile, typename Policy>
TA::Array<double, 4, Tile, Policy> g_lt(
    TA::Array<double, 4, Tile, Policy> &abij, const Eigen::VectorXd &ens,
    std::size_t n_occ, std::size_t n_frozen, double x) {

  auto convert = [&ens, n_occ, n_frozen, x](Tile &result_tile,
                                         const Tile &arg_tile) {

    result_tile = Tile(arg_tile.range());

    // compute index
    const auto d0 = result_tile.range().lobound()[0];
    const auto dn = result_tile.range().upbound()[0];
    const auto a0 = result_tile.range().lobound()[1];
    const auto an = result_tile.range().upbound()[1];
    const auto b0 = result_tile.range().lobound()[2];
    const auto bn = result_tile.range().upbound()[2];
    const auto i0 = result_tile.range().lobound()[3];
    const auto in = result_tile.range().upbound()[3];

    double alpha = 3.0*(ens(n_occ) - ens(n_occ - 1));

    auto tile_idx = 0;
    typename Tile::value_type norm = 0.0;
    for (auto d = d0; d < dn; ++d) {
      for (auto a = a0; a < an; ++a) {
        const auto e_a = ens[a + n_occ];
        for (auto b = b0; b < bn; ++b) {
          const auto e_b = ens[b + n_occ];
          for (auto i = i0; i < in; ++i, ++tile_idx) {
            const auto e_i = ens[i + n_frozen];
            const auto exponent = pow(x,0.5*(e_a/alpha - 1.0/6.0))*pow(x,0.5*(e_b/alpha - 1.0/6.0))*pow(x,0.5*(-e_i/alpha - 1.0/6.0));
            const double result = arg_tile[tile_idx]*exponent;
            result_tile[tile_idx] = result;
            norm += result*result;

          }
        }
      }
    }
    return std::sqrt(norm);
  };

  auto result = TA::foreach(abij, convert);
  abij.world().gop.fence();
  return result;
}

template <typename Tile, typename Policy>
TA::Array<double, 4, Tile, Policy> g_c_lt(
    TA::Array<double, 4, Tile, Policy> &abij, const Eigen::VectorXd &ens,
    std::size_t n_occ, std::size_t n_frozen, double x) {

  auto convert = [&ens, n_occ, n_frozen, x](Tile &result_tile,
                                         const Tile &arg_tile) {

    result_tile = Tile(arg_tile.range());

    // compute index
    const auto c0 = result_tile.range().lobound()[0];
    const auto cn = result_tile.range().upbound()[0];
    const auto i0 = result_tile.range().lobound()[1];
    const auto in = result_tile.range().upbound()[1];
    const auto j0 = result_tile.range().lobound()[2];
    const auto jn = result_tile.range().upbound()[2];
    const auto k0 = result_tile.range().lobound()[3];
    const auto kn = result_tile.range().upbound()[3];

    double alpha = 3.0*(ens(n_occ) - ens(n_occ - 1));

    auto tile_idx = 0;
    typename Tile::value_type norm = 0.0;
    for (auto c = c0; c < cn; ++c) {
      const auto e_c = ens[c + n_occ];
      for (auto i = i0; i < in; ++i) {
        const auto e_i = ens[i + n_frozen];
        for (auto j = j0; j < jn; ++j) {
          const auto e_j = ens[j + n_frozen];
          for (auto k = k0; k < kn; ++k, ++tile_idx) {
            const auto exponent = pow(x,0.5*(e_c/alpha - 1.0/6.0))*pow(x,0.5*(-e_i/alpha - 1.0/6.0))*pow(x,0.5*(-e_j/alpha - 1.0/6.0));
            const double result = arg_tile[tile_idx]*exponent;
            result_tile[tile_idx] = result;
            norm += result*result;

          }
        }
      }
    }
    return std::sqrt(norm);
  };

  auto result = TA::foreach(abij, convert);
  abij.world().gop.fence();
  return result;
}

template <typename Tile, typename Policy>
TA::Array<double, 4, Tile, Policy> g_abij_ltfcn(
    TA::Array<double, 4, Tile, Policy> &abij, const Eigen::VectorXd &ens,
    std::size_t n_occ, std::size_t n_frozen, double x) {
  auto convert = [&ens, n_occ, n_frozen, x](Tile &result_tile,
                                         const Tile &arg_tile) {

    double alpha = 3.0*(ens(n_occ) - ens(n_occ - 1));
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
            const auto exponent = pow(x,0.5*(e_a/alpha - 1.0/6.0))*pow(x,0.5*(e_b/alpha - 1.0/6.0))
                *pow(x,0.5*(-e_i/alpha - 1.0/6.0))*pow(x,0.5*(-e_j/alpha - 1.0/6.0));
            const double result = arg_tile[tile_idx]*exponent;
            result_tile[tile_idx] = result;
            norm += result*result;
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

template <typename Tile, typename Policy>
TA::Array<double, 4, Tile, Policy> t2_lt(const
    TA::Array<double, 4, Tile, Policy> &abij, const Eigen::VectorXd &ens,
    std::size_t n_occ, std::size_t n_frozen, double x) {

  auto convert = [&ens, n_occ, n_frozen, x](Tile &result_tile,
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

    double alpha = 3.0*(ens(n_occ) - ens(n_occ - 1));

    auto tile_idx = 0;
    typename Tile::value_type norm = 0.0;
    for (auto a = a0; a < an; ++a) {
      for (auto b = b0; b < bn; ++b) {
        const auto e_b = ens[b + n_occ];
        for (auto i = i0; i < in; ++i) {
          const auto e_i = ens[i + n_frozen];
          for (auto j = j0; j < jn; ++j, ++tile_idx) {
            const auto e_j = ens[j + n_frozen];
            const auto exponent = pow(x,0.5*(e_b/alpha - 1.0/6.0))*pow(x,0.5*(-e_i/alpha - 1.0/6.0))*pow(x,0.5*(-e_j/alpha - 1.0/6.0));
            const double result = arg_tile[tile_idx]*exponent;
            result_tile[tile_idx] = result;
            norm += result*result;

          }
        }
      }
    }
    return std::sqrt(norm);
  };

  auto result = TA::foreach(abij, convert);
  abij.world().gop.fence();
  return result;
}

template <typename Tile, typename Policy>
TA::Array<double, 2, Tile, Policy> t1_ltfcn(const
    TA::Array<double, 2, Tile, Policy> &abij, const Eigen::VectorXd &ens,
    std::size_t n_occ, std::size_t n_frozen, double x) {

  auto convert = [&ens, n_occ, n_frozen, x](Tile &result_tile,
                                         const Tile &arg_tile) {

    result_tile = Tile(arg_tile.range());

    // compute index
    const auto a0 = result_tile.range().lobound()[0];
    const auto an = result_tile.range().upbound()[0];
    const auto i0 = result_tile.range().lobound()[1];
    const auto in = result_tile.range().upbound()[1];

    double alpha = 3.0*(ens(n_occ) - ens(n_occ - 1));

    auto tile_idx = 0;
    typename Tile::value_type norm = 0.0;
    for (auto a = a0; a < an; ++a) {
      const auto e_a = ens[a + n_occ];
      for (auto i = i0; i < in; ++i, ++tile_idx) {
        const auto e_i = ens[i + n_frozen];
        const auto exponent = pow(x,0.5*(e_a/alpha - 1.0/6.0))*pow(x,0.5*(-e_i/alpha - 1.0/6.0));
        const double result = arg_tile[tile_idx]*exponent;
        result_tile[tile_idx] = result;
        norm += result*result;
      }
    }
    return std::sqrt(norm);
  };

  auto result = TA::foreach(abij, convert);
  abij.world().gop.fence();
  return result;
}

template <typename Tile, typename Policy>
TA::Array<double, 4, Tile, Policy> t22_lt(const
    TA::Array<double, 4, Tile, Policy> &abij, const Eigen::VectorXd &ens,
    std::size_t n_occ, std::size_t n_frozen, double x) {

  auto convert = [&ens, n_occ, n_frozen, x](Tile &result_tile,
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

    double alpha = 3.0*(ens(n_occ) - ens(n_occ - 1));

    auto tile_idx = 0;
    typename Tile::value_type norm = 0.0;
    for (auto a = a0; a < an; ++a) {
      const auto e_a = ens[a + n_occ];
      for (auto b = b0; b < bn; ++b) {
        const auto e_b = ens[b + n_occ];
        for (auto i = i0; i < in; ++i) {
          const auto e_i = ens[i + n_frozen];
          for (auto j = j0; j < jn; ++j, ++tile_idx) {
            const auto exponent = pow(x,0.5*(e_b/alpha - 1.0/6.0))*pow(x,0.5*(-e_i/alpha - 1.0/6.0))*pow(x,0.5*(e_a/alpha - 1.0/6.0));
            const double result = arg_tile[tile_idx]*exponent;
            result_tile[tile_idx] = result;
            norm += result*result;

          }
        }
      }
    }
    return std::sqrt(norm);
  };

  auto result = TA::foreach(abij, convert);
  abij.world().gop.fence();
  return result;
}

// create matrix d("a,i") = 1/(ei - ea)
template <typename Tile, typename Policy>
TA::DistArray<
    Tile, typename std::enable_if<std::is_same<Policy, TA::SparsePolicy>::value,
                                  TA::SparsePolicy>::type>
create_d_ai(madness::World &world, const TA::TiledRange &trange,
            const Eigen::VectorXd &ens, std::size_t n_occ,
            std::size_t n_frozen) {
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

}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_MBPT_DENOM_H_
