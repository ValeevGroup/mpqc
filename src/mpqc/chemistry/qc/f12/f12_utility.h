//
// Created by Chong Peng on 10/21/15.
//

#ifndef MPQC_UTILITY_H
#define MPQC_UTILITY_H

#include <string>
#include "../../../../../include/eigen.h"
#include "../../../../../include/tiledarray.h"
#include <vector>
#include <TiledArray/error.h>
#include <TiledArray/sparse_shape.h>
#include <TiledArray/tiled_range1.h>

namespace mpqc{
namespace f12{

// coeffs of intermediates linear in the geminal, i.e. V and C
constexpr double C_ijij = (1./2 + 1./4)/2;  // 3/8
constexpr double C_ijji = (1./2 - 1./4)/2;  // 1/8
// coeffs of intermediates quadratic in the geminal, i.e. X and B
constexpr double CC_ijij = C_ijij * C_ijij + C_ijji * C_ijji;  // 10/64
constexpr double CC_ijji = 2 * C_ijij * C_ijji;  // C_ijij * C_ijji + C_ijji * C_ijij;
                                             // 6/64
// closed-shell-antisymmetrized versions of above
constexpr double C_ijij_bar = 2 * C_ijij - C_ijji;  // 5/8
constexpr double C_ijji_bar = 2 * C_ijji - C_ijij;  // -1/8
constexpr double CC_ijij_bar = 2 * CC_ijij - CC_ijji;  // 14/64
constexpr double CC_ijji_bar = 2 * CC_ijji - CC_ijij;  // 2/64


    double basis_to_f12exponent(const std::string& basis_name);

    std::vector<std::pair<double,double>> stg_ng_fit(std::size_t n, double zeta);

    std::vector<std::pair<double,double>> gtg_params_squared(const std::vector<std::pair<double,double>>& pragmas);

    // Gaussian Type Geminal parameters
    struct GTGParams{

        double exponent;
        std::size_t n_fit;

        GTGParams() = default;
        GTGParams(double zeta, std::size_t n = 6) : exponent(zeta), n_fit(n) { }
        GTGParams(std::string basis_name, std::size_t n = 6) : exponent(basis_to_f12exponent(basis_name)), n_fit(n) {}

        std::vector<std::pair<double,double>> compute()
        {
            if(exponent==0){
                return std::vector<std::pair<double,double>>();
            }
            else{
                return stg_ng_fit(n_fit,exponent);
            }
        };
    };

    // make sure the occ is blocked by 1!
    TiledArray::SparseShape<float> make_ijij_ijji_shape(const TiledArray::TiledRange& );

template <typename Tile, typename Policy>
void convert_X_ijkl(TiledArray::Array<double, 4, Tile, Policy> &ijkl,
                                          const Eigen::MatrixXd &F)
{
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
            const auto f_ii = F(i,i);
            for (auto j = j0; j < jn; ++j) {
                const auto f_jj = F(j,j);
                for (auto k = k0; k < kn; ++k) {
                    for (auto l = l0; l < ln; ++l, ++tile_idx) {
                        const auto old = result_tile[tile_idx];
                        const auto result_abij = old*(f_ii+f_jj);
                        norm += result_abij*result_abij;
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
TiledArray::Array<double, 4, Tile, Policy> convert_C_ijab(TiledArray::Array<double, 4, Tile, Policy> &ijab,
                    const std::size_t n_occ, const Eigen::VectorXd &ens)
{
    auto convert = [&ens, n_occ](Tile &result_tile, const Tile& arg_tile) {

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
            auto en_i = ens[i];
            for (auto j = j0; j < jn; ++j) {
                auto en_j = ens[j];
                for (auto a = a0; a < an; ++a) {
                    auto en_a = ens[a+n_occ];
                    for (auto b = b0; b < bn; ++b, ++tile_idx) {
                        auto en_b = ens[b+n_occ];
                        const auto old = arg_tile[tile_idx];
                        const auto divde = 1/(en_i + en_j - en_a - en_b);
                        const auto result_abij = old*divde;
                        norm += result_abij*result_abij;
                        result_tile[tile_idx] = result_abij;
                    }
                }
            }
        }
        return std::sqrt(norm);
    };

    return TiledArray::foreach(ijab, convert);
}

/// tile functor that computes contribution to the F12 energy

/// Give (a tile of) an F12 theory intermediate (V, B, X, and C)
/// computes the energy in the diagonal approximation to the F12 energy.
/// The energy contribution from intermediate element \f$ I^{ij}_{kl} \f$ is proportional
/// to \f$ c^{ij}_{ij} I^{ij}_{ij} + c^{ji}_{ij} I^{ij}_{ji} \f$.
/// NB for \f$ i=j \f$ both terms contribute.
template <typename Tile>
struct CLF12Energy {
    using result_type =  double;
    using argument_type =  Tile;

    double c_ijij;
    double c_ijji;

    CLF12Energy() = default;
    CLF12Energy(CLF12Energy const &) = default;
    CLF12Energy(double cijij, double cijji) : c_ijij(cijij), c_ijji(cijji) {}

    result_type operator()() const { return 0.0; }

    result_type operator()(result_type const &t) const { return t; }

    void operator()(result_type &me, result_type const &other) const {
        me += other;
    }

    void operator()(result_type &me, argument_type const &tile) const {
        auto const &range = tile.range();

        TA_ASSERT(range.rank() == 4);

        auto const st = range.lobound_data();
        auto const fn = range.upbound_data();
        auto tile_idx = 0;

        auto sti = st[0];
        auto fni = fn[0];
        auto stj = st[1];
        auto fnj = fn[1];
        auto stk = st[2];
        auto fnk = fn[2];
        auto stl = st[3];
        auto fnl = fn[3];

        const auto* value_ptr = tile.data();
        for (auto i = sti; i < fni; ++i) {
          for (auto j = stj; j < fnj; ++j) {
            for (auto k = stk; k < fnk; ++k) {
              for (auto l = stl; l < fnl; ++l, ++value_ptr) {
                // ijij
                if (k == i && l == j) {
                  me += c_ijij * (*value_ptr);
                }
                // ijji
                if (l == i && k == j) {
                  me += c_ijji * (*value_ptr);
                }
              }
            }
          }
        }
    }
};

}
}
#endif //MPQC_UTILITY_H
