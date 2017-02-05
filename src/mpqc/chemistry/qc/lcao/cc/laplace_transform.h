//
// Created by Fabijan Pavosevic on 02/4/17.
//

//#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_CC_CCSD_H_
//#define MPQC4_SRC_MPQC_CHEMISTRY_QC_CC_CCSD_H_

#include <tiledarray.h>

#include "mpqc/math/external/eigen/eigen.h"

//this set of functions is used to re-scale 2-electron integrals with the Laplace transformat (exp(-orb_energy)).
namespace mpqc {

//re-scaling of g_dabi integral with the exponents of orbital energies. One occupied and two unoccupied orbitals are re-scaled (i,a,b).
template <typename Tile, typename Policy>
TA::Array<double, 4, Tile, Policy> g_dabi_laplace_transform(
    TA::Array<double, 4, Tile, Policy> &dabi, const Eigen::VectorXd &ens,
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

  auto result = TA::foreach(dabi, convert);
  dabi.world().gop.fence();
  return result;
}

//re-scaling of g_cjkl integral with the exponents of orbital energies. Two occupied and one unoccupied orbitals are re-scaled (i,j,c).
template <typename Tile, typename Policy>
TA::Array<double, 4, Tile, Policy> g_cjkl_laplace_transform(
    TA::Array<double, 4, Tile, Policy> &cjkl, const Eigen::VectorXd &ens,
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

  auto result = TA::foreach(cjkl, convert);
  cjkl.world().gop.fence();
  return result;
}

//re-scaling of g_abij integral with the exponents of orbital energies. Two occupied and two unoccupied orbitals are re-scaled (i,j,a,b).
template <typename Tile, typename Policy>
TA::Array<double, 4, Tile, Policy> g_abij_laplace_transform(
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

//re-scaling of t2 amplitudes with the exponents of orbital energies. Two occupied and one unoccupied orbitals are re-scaled (i,j,a).
template <typename Tile, typename Policy>
TA::Array<double, 4, Tile, Policy> t2_oou_laplace_transform(const
    TA::Array<double, 4, Tile, Policy> &t2, const Eigen::VectorXd &ens,
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

  auto result = TA::foreach(t2, convert);
  t2.world().gop.fence();
  return result;
}

//re-scaling of t2 amplitudes with the exponents of orbital energies. One occupied and two unoccupied orbitals are re-scaled (i,a,b).
template <typename Tile, typename Policy>
TA::Array<double, 4, Tile, Policy> t2_ouu_laplace_transform(const
    TA::Array<double, 4, Tile, Policy> &t2, const Eigen::VectorXd &ens,
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

  auto result = TA::foreach(t2, convert);
  t2.world().gop.fence();
  return result;
}

//re-scaling of t1 amplitudes with the exponents of orbital energies. One occupied and one unoccupied orbitals are re-scaled (i,a).
template <typename Tile, typename Policy>
TA::Array<double, 2, Tile, Policy> t1_laplace_transform(const
    TA::Array<double, 2, Tile, Policy> &t1, const Eigen::VectorXd &ens,
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

  auto result = TA::foreach(t1, convert);
  t1.world().gop.fence();
  return result;
}

}  // namespace mpqc

//#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_MBPT_DENOM_H_

