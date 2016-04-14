//
// Created by Chong Peng on 10/21/15.
//

#ifndef TILECLUSTERCHEM_UTILITY_H
#define TILECLUSTERCHEM_UTILITY_H

#include <string>
#include "../include/eigen.h"
#include "../include/tiledarray.h"
#include <vector>
#include <TiledArray/error.h>
#include <TiledArray/sparse_shape.h>
#include <TiledArray/tiled_range1.h>

namespace mpqc{
namespace f12{

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

}
}
#endif //TILECLUSTERCHEM_UTILITY_H
