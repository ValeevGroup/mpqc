#pragma once
#ifndef TCC_INTEGRALS_TATOLOWRANKTENSOR_H
#define TCC_INTEGRALS_TATOLOWRANKTENSOR_H

#include "../tensor/tile_pimpl.h"
#include "../tensor/tile_variant.h"
#include "../tensor/low_rank_tile.h"
#include "../tensor/full_rank_tile.h"
#include "../tensor/tile_algebra.h"

#include "../include/eigen.h"
#include "../include/tiledarray.h"


#include "tile_engine.h"

namespace tcc {
namespace integrals {
namespace compute_functors {

/// Function to convert btas::Tesors into Tile_Pimpl objects

/// Will convert any 2, 3 or 4 dimensional btas::Tensor into a TilePimpl.
/// Flattening of order 3 tensors is (0, 1*2).
/// Flattening of order 4 tensors is (0*1, 2*3).
template<std::size_t DIM>
class TaToLowRankTensor {
  public:
    using TileType = tensor::TilePimpl<double>;

    TaToLowRankTensor() = default;
    TaToLowRankTensor(double cut) : cut_(cut) {}

    TileType operator()(TiledArray::Tensor<double> const &tat) const {
        Eigen::MatrixXd Tile = eigen_map(tat);
        Eigen::MatrixXd L, R;
        bool is_full = algebra::Decompose_Matrix(Tile, L, R, cut_);
        if (!is_full) {
            tensor::TileVariant<double> tile_variant{
                tensor::LowRankTile<double>{std::move(L), std::move(R)}};

            return tensor::TilePimpl<double>{tat.range(),
                                             std::move(tile_variant), cut_};
        }

        tensor::TileVariant<double> tile_variant{
            tensor::FullRankTile<double>{std::move(Tile)}};

        return tensor::TilePimpl<double>{tat.range(), std::move(tile_variant), cut_};
    }

  private:
    using RowMajorMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                                         Eigen::RowMajor>;

    Eigen::MatrixXd eigen_map(TiledArray::Tensor<double> const &tat) const;

    double cut_ = 1e-7;
};

template <>
Eigen::MatrixXd
TaToLowRankTensor<2>::eigen_map(TiledArray::Tensor<double> const &tat) const {
    auto const &extent = tat.range().size();
    return Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic,
                                          Eigen::Dynamic, Eigen::RowMajor>,
                      Eigen::AutoAlign>(tat.data(), extent[0], extent[1]);
}
template <>
Eigen::MatrixXd
TaToLowRankTensor<3>::eigen_map(TiledArray::Tensor<double> const &tat) const {
    auto const &extent = tat.range().size();
    return Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic,
                                          Eigen::Dynamic, Eigen::RowMajor>,
                      Eigen::AutoAlign>(tat.data(), extent[0],
                                        extent[1] * extent[2]);
}
template <>
Eigen::MatrixXd
TaToLowRankTensor<4>::eigen_map(TiledArray::Tensor<double> const &tat) const {
    auto const &extent = tat.range().size();
    return Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic,
                                          Eigen::Dynamic, Eigen::RowMajor>,
                      Eigen::AutoAlign>(tat.data(), extent[0] * extent[1],
                                        extent[2] * extent[3]);
}


} // namespace compute_functors
} // namespace integrals
} // namespace tcc


#endif /* end of include guard: TCC_INTEGRALS_TATOLOWRANKTENSOR_H */
