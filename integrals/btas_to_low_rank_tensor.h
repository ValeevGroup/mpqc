#pragma once
#ifndef TCC_INTEGRALS_BTASTOLOWRANKTENSOR_H
#define TCC_INTEGRALS_BTASTOLOWRANKTENSOR_H

#include "../tensor/tile_pimpl.h"
#include "../tensor/tile_variant.h"
#include "../tensor/low_rank_tile.h"
#include "../tensor/full_rank_tile.h"
#include "../tensor/tile_algebra.h"

#include "../include/btas.h"
#include "../include/eigen.h"
#include "../include/tiledarray.h"

#include "../tensor/btas_shallow_copy_wrapper.h"

#include "tile_engine.h"

namespace tcc {
namespace integrals {
namespace compute_functors {

/// Function to convert btas::Tesors into Tile_Pimpl objects

/// Will convert any 2, 3 or 4 dimensional btas::Tensor into a TilePimpl.
/// Flattening of order 3 tensors is (0, 1*2).
/// Flattening of order 4 tensors is (0*1, 2*3).
class BtasToLowRankTensor {
  public:
    using TileType = tensor::TilePimpl<double>;

    BtasToLowRankTensor() = default;
    BtasToLowRankTensor(double cut) : cut_(cut) {}

    template <std::size_t N>
    TileType operator()(tensor::ShallowTensor<N> const &bt) const {
        RowMatrixXd Tile = eigen_map(bt.tensor());
        if (Tile.lpNorm<2>() >= TiledArray::SparseShape<float>::threshold()) {
            RowMatrixXd L, R;
            bool is_full = algebra::Decompose_Matrix(Tile, L, R, cut_);
            if (!is_full) {
                tensor::TileVariant<double> tile_variant{
                    tensor::LowRankTile<double>{std::move(L), std::move(R)}};

                return tensor::TilePimpl<double>{bt.range(),
                                                 std::move(tile_variant), cut_};
            }

            tensor::TileVariant<double> tile_variant{
                tensor::FullRankTile<double>{std::move(Tile)}};

            return tensor::TilePimpl<double>{bt.range(),
                                             std::move(tile_variant)};
        } else { // For now if tile is super sparse return a 1x1 zero tile.
            Tile = RowMatrixXd::Zero(1, 1);
            tensor::TileVariant<double> tile_variant{
                tensor::FullRankTile<double>{std::move(Tile)}};

            return tensor::TilePimpl<double>{bt.range(),
                                             std::move(tile_variant)};
        }
    }

  private:
    template <std::size_t N>
    RowMatrixXd eigen_map(TileEngine<double>::TileType<N> const &btas_t) const;

    double cut_ = 1e-7;
};

template <>
RowMatrixXd BtasToLowRankTensor::eigen_map<2>(
    TileEngine<double>::TileType<2> const &btas_t) const {
    auto const &extent = btas_t.extent();
    return eMap<const RowMatrixXd>(btas_t.data(), extent[0], extent[1]);
}

template <>
RowMatrixXd BtasToLowRankTensor::eigen_map<3>(
    TileEngine<double>::TileType<3> const &btas_t) const {
    auto const &extent = btas_t.extent();
    return eMap<const RowMatrixXd>(btas_t.data(), extent[0], extent[1] * extent[2]);
}

template <>
RowMatrixXd BtasToLowRankTensor::eigen_map<4>(
    TileEngine<double>::TileType<4> const &btas_t) const {
    auto const &extent = btas_t.extent();
    return eMap<const RowMatrixXd>(btas_t.data(), extent[0] * extent[1],
                             extent[2] * extent[3]);
}

} // namespace compute_functors
} // namespace integrals
} // namespace tcc


#endif /* end of include guard: TCC_INTEGRALS_BTASTOLOWRANKTENSOR_H */
