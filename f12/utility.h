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
        int n_fit;

        GTGParams(double zeta, int n = 6) : exponent(zeta), n_fit(n) { }
        GTGParams(std::string basis_name, int n = 6) : exponent(basis_to_f12exponent(basis_name)), n_fit(n) {}

        std::vector<std::pair<double,double>> compute()
        {
            return stg_ng_fit(n_fit,exponent);
        };
    };

    // make sure the occ is blocked by 1!
    TiledArray::SparseShape<float> make_ijij_shape(const TiledArray::TiledRange& );
    TiledArray::SparseShape<float> make_ijji_shape(const TiledArray::TiledRange& );

    template <typename Tile>
    struct DiagonalSum {
        using result_type =  double;
        using argument_type =  Tile;

        DiagonalSum() = default;
        DiagonalSum(DiagonalSum const &) = default;

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
            for (auto i = st[0]; i < fn[0]; ++i) {
                for (auto j = st[1]; j < fn[1]; ++j) {
                    for (auto k = st[2]; k < fn[2]; ++k) {
                        for (auto l = st[3]; l < fn[3]; ++l, ++tile_idx) {
                            if((i == j) &&  (i == k) && (k == l)){
                                me += tile.data()[tile_idx];
                            }
                        }
                    }
                }
            }
        }
    };


template <typename Tile>
struct OffDiagonalSum {
    using result_type =  double;
    using argument_type =  Tile;

    OffDiagonalSum() = default;
    OffDiagonalSum(OffDiagonalSum const &) = default;

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
        for (auto i = st[0]; i < fn[0]; ++i) {
            for (auto j = st[1]; j < fn[1]; ++j) {
                for (auto k = st[2]; k < fn[2]; ++k) {
                    for (auto l = st[3]; l < fn[3]; ++l, ++tile_idx) {
                        if(i > j && k==i && l==j){
                            me += tile.data()[tile_idx];
                        }
                    }
                }
            }
        }
    }
};


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
            for (auto j = j0; j < jn; ++j) {
                for (auto k = k0; k < kn; ++k) {
                    for (auto l = l0; l < ln; ++l, ++tile_idx) {
                        const auto old = result_tile[tile_idx];
                        const auto f_ik = F(i,k);
                        const auto f_jl = F(j,l);
                        const auto result_abij = old*(f_ik+f_jl);
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


}
}
#endif //TILECLUSTERCHEM_UTILITY_H
