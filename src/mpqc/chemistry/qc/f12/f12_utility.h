//
// Created by Chong Peng on 10/21/15.
//

#ifndef MPQC_UTILITY_H
#define MPQC_UTILITY_H

#include <TiledArray/error.h>
#include <TiledArray/sparse_shape.h>
#include <TiledArray/tiled_range1.h>
#include <string>
#include <tiledarray.h>
#include <vector>

#include "mpqc/math/external/eigen/eigen.h"

namespace mpqc {
namespace f12 {

// coeffs of intermediates linear in the geminal, i.e. V and C
constexpr double C_ijij = (1. / 2 + 1. / 4) / 2;  // 3/8
constexpr double C_ijji = (1. / 2 - 1. / 4) / 2;  // 1/8
// coeffs of intermediates quadratic in the geminal, i.e. X and B
constexpr double CC_ijij = C_ijij * C_ijij + C_ijji * C_ijji;  // 10/64
constexpr double CC_ijji =
    2 * C_ijij * C_ijji;  // C_ijij * C_ijji + C_ijji * C_ijij;
                          // 6/64
// closed-shell-antisymmetrized versions of above
constexpr double C_ijij_bar = 2 * C_ijij - C_ijji;     // 5/8
constexpr double C_ijji_bar = 2 * C_ijji - C_ijij;     // -1/8
constexpr double CC_ijij_bar = 2 * CC_ijij - CC_ijji;  // 14/64
constexpr double CC_ijji_bar = 2 * CC_ijji - CC_ijij;  // 2/64

double basis_to_f12exponent(const std::string &basis_name);

std::vector<std::pair<double, double>> stg_ng_fit(std::size_t n, double zeta);

std::vector<std::pair<double, double>> gtg_params_squared(
    const std::vector<std::pair<double, double>> &pragmas);

// Gaussian Type Geminal parameters
struct GTGParams {
  double exponent;
  std::size_t n_fit;

  GTGParams() = default;
  GTGParams(double zeta, std::size_t n = 6) : exponent(zeta), n_fit(n) {}
  GTGParams(std::string basis_name, std::size_t n = 6)
      : exponent(basis_to_f12exponent(basis_name)), n_fit(n) {}

  std::vector<std::pair<double, double>> compute() {
    if (exponent == 0) {
      return std::vector<std::pair<double, double>>();
    } else {
      return stg_ng_fit(n_fit, exponent);
    }
  };
};

/// computes shape of F12 intermediates in "diagonal" approximation

/// F12 theory contains several intermediates which in the diagonal form
/// are sparse: only elements \f$ I^{ij}_{ij} \f$ and \f$ I^{ji}_{ij} \f$
/// of the intermediate \f$ I^{ij}_{kl} \f$  must be computed. Therefore
/// only shape elements [i,j,i,j] and [i,j,j,i] are needed.
/// @tparam P1_eq_P2 if \c true , restrict i <= j
/// @note make sure the occ is blocked by 1!
/// TODO must rewrite intermediate expressions to support P1_eq_P2==true, e.g.
/// <ij|ma'><ma'|r|ij> term
///      must become <ij|ma'><ma'|r|ij> + <ij|a'm><a'm|r|ij>
template <bool P1_eq_P2 = false>
TiledArray::SparseShape<float> make_ijij_ijji_shape(
    const TiledArray::TiledRange &trange) {
  // number of occ tiles
  auto n_occ = trange.data()[0].tiles_range().second;

  TiledArray::Tensor<float> tile_norms(trange.tiles_range(), 0.0);

  float max = (float) n_occ*n_occ;

  // set sparse tile
  for (auto i = 0; i < n_occ; i++) {
    const auto j_start = P1_eq_P2 ? i : 0;
    for (auto j = j_start; j < n_occ; j++) {
      const auto ijij = ((i * n_occ + j) * n_occ + i) * n_occ + j;
      const auto ijji = ((i * n_occ + j) * n_occ + j) * n_occ + i;
      tile_norms[ijij] = max;
      tile_norms[ijji] = max;
    }
  }

  TiledArray::SparseShape<float> shape(tile_norms, trange);
  return shape;
}

template <typename Tile, typename Policy>
void convert_X_ijkl(TiledArray::Array<double, 4, Tile, Policy> &ijkl,
                    const Eigen::MatrixXd &F) {
  auto convert = [&F](Tile &result_tile) {

    // compute index
    const auto i0 = result_tile.range().lobound()[0];
    const auto in = result_tile.range().upbound()[0];
    const auto j0 = result_tile.range().lobound()[1];
    const auto jn = result_tile.range().upbound()[1];
    const auto k0 = result_tile.range().lobound()[2];
    const auto kn = result_tile.range().upbound()[2];
    const auto l0 = result_tile.range().lobound()[3];
    const auto ln = result_tile.range().upbound()[3];

    auto tile_idx = 0;
    typename Tile::value_type norm = 0.0;
    for (auto i = i0; i < in; ++i) {
      const auto f_ii = F(i, i);
      for (auto j = j0; j < jn; ++j) {
        const auto f_jj = F(j, j);
        for (auto k = k0; k < kn; ++k) {
          for (auto l = l0; l < ln; ++l, ++tile_idx) {
            const auto old = result_tile[tile_idx];
            const auto result_abij = old * (f_ii + f_jj);
            norm += result_abij * result_abij;
            result_tile[tile_idx] = result_abij;
          }
        }
      }
    }
    return std::sqrt(norm);
  };

  TiledArray::foreach_inplace(ijkl, convert);
}

template <typename Tile, typename Policy>
TiledArray::Array<double, 4, Tile, Policy> convert_C_ijab(
    TiledArray::Array<double, 4, Tile, Policy> &ijab, const std::size_t n_occ,
    const std::size_t n_frozen, const Eigen::VectorXd &ens) {
  auto convert = [&ens, n_occ, n_frozen](Tile &result_tile,
                                         const Tile &arg_tile) {

    result_tile = Tile(arg_tile.range());

    // compute index
    const auto i0 = result_tile.range().lobound()[0];
    const auto in = result_tile.range().upbound()[0];
    const auto j0 = result_tile.range().lobound()[1];
    const auto jn = result_tile.range().upbound()[1];
    const auto a0 = result_tile.range().lobound()[2];
    const auto an = result_tile.range().upbound()[2];
    const auto b0 = result_tile.range().lobound()[3];
    const auto bn = result_tile.range().upbound()[3];

    auto tile_idx = 0;
    typename Tile::value_type norm = 0.0;
    for (auto i = i0; i < in; ++i) {
      auto en_i = ens[i + n_frozen];
      for (auto j = j0; j < jn; ++j) {
        auto en_j = ens[j + n_frozen];
        for (auto a = a0; a < an; ++a) {
          auto en_a = ens[a + n_occ];
          for (auto b = b0; b < bn; ++b, ++tile_idx) {
            auto en_b = ens[b + n_occ];
            const auto old = arg_tile[tile_idx];
            const auto divde = 1 / (en_i + en_j - en_a - en_b);
            const auto result_abij = old * divde;
            norm += result_abij * result_abij;
            result_tile[tile_idx] = result_abij;
          }
        }
      }
    }
    return std::sqrt(norm);
  };

  return TiledArray::foreach (ijab, convert);
}

/// tile reductor that computes pair contributions to the F12 energy

/// Given (a tile of) an F12 theory intermediate (V, B, X, and C)
/// computes its contribution to the pair energies in the diagonal F12
/// approximation.
/// The energy contribution from intermediate element \f$ I^{ij}_{kl} \f$ is
/// proportional
/// to \f$ c^{ij}_{ij} I^{ij}_{ij} + c^{ji}_{ij} I^{ij}_{ji} \f$, (NB for \f$
/// i=j \f$ both
/// terms contribute). The coefficients will
/// depend on the definition of the geminal, will differ between intermediates
/// linear and
/// quadratic in the geminal; they will also vary with the reference case, e.g.
/// closed-shell
/// spin-restricted case, open-shell, or spin-unrestricted cases.
/// -# if pairspin == none (closed-shell reference) will compute pair energies
/// for i<=j
/// -# if pairspin == aa/bb will compute pair energies for i<j
/// -# if pairspin == ab will compute pair energies for all ij
/// -# if pairspin == s0 will compute pair energies for all i<=j
/// -# if pairspin == s1 will compute pair energies for all i<j
///
template <typename Tile>
struct F12PairEnergyReductor {
  using result_type = RowMatrix<typename Tile::scalar_type>;
  using argument_type = Tile;

  enum class PairSpin {
    none,  // closed-shell only
    aa,    // same-spin
    bb,    // same-spin
    ab,    // opposite-spin
    s0,    // S=0 (singlet pair)
    s1     // S=1 (triplet pair)
  };
  double c_ijij;
  double c_ijji;
  int n1;             //!< # of geminal orbitals for electron 1
  int n2;             //!< # of geminal orbitals for electron 2
  PairSpin spincase;  //!< spin case

  F12PairEnergyReductor() = default;
  F12PairEnergyReductor(F12PairEnergyReductor const &) = default;
  F12PairEnergyReductor(double cijij, double cijji, int _n1, int _n2 = -1,
                        PairSpin sc = PairSpin::none)
      : c_ijij(cijij),
        c_ijji(cijji),
        n1(_n1),
        n2(_n2 == -1 ? _n1 : _n2),
        spincase(sc) {
    TA_USER_ASSERT(spincase != PairSpin::s0 && spincase != PairSpin::s1,
                   "spin-adapted F12PairEnergyReductor not implemented yet");
    TA_USER_ASSERT(spincase != PairSpin::aa && spincase != PairSpin::bb &&
                       spincase != PairSpin::ab,
                   "spin-orbital F12PairEnergyReductor not implemented yet");
    if (spincase != PairSpin::ab)
      TA_USER_ASSERT(n1 == n2,
                     "same-spin case but #s of orbitals for "
                     "electrons 1 and 2 differ");
  }

  result_type operator()() const { return result_type::Zero(n1, n2); }

  result_type operator()(result_type const &t) const { return t; }

  void operator()(result_type &me, result_type const &other) const {
    TA_USER_ASSERT(me.rows() == other.rows() && me.cols() == other.cols(),
                   "result vector size mismatch");
    me += other;
  }

  void operator()(result_type &me, argument_type const &tile) const {
    auto const &range = tile.range();

    TA_ASSERT(range.rank() == 4);

    auto const st = range.lobound_data();
    auto const fn = range.upbound_data();

    auto sti = st[0];
    auto fni = fn[0];
    auto stj = st[1];
    auto fnj = fn[1];
    auto stk = st[2];
    auto fnk = fn[2];
    auto stl = st[3];
    auto fnl = fn[3];

    const auto *value_ptr = tile.data();
    for (auto i = sti; i < fni; ++i) {
      for (auto j = stj; j < fnj; ++j) {
        const auto ii = std::min(i, j);
        const auto jj = std::max(i, j);
        for (auto k = stk; k < fnk; ++k) {
          for (auto l = stl; l < fnl; ++l, ++value_ptr) {
            // ijij
            if (k == i && l == j) {
              me(ii, jj) += c_ijij * (*value_ptr);
            }
            // ijji
            if (l == i && k == j) {
              me(ii, jj) += c_ijji * (*value_ptr);
            }
          }
        }
      }
    }
  }
};
}
}
#endif  // MPQC_UTILITY_H
