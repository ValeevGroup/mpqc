#ifndef TileClusterChem_TILE_OPS_H
#define TileClusterChem_TILE_OPS_H

#include "low_rank_tile.h"
#include "full_rank_tile.h"
#include "tile_algebra.h"

namespace tile_ops {

template <typename T>
FullRankTile<T>
gemm(const FullRankTile<T> &left, const FullRankTile<T> &right, double alpha) {
    return FullRankTile
        <T>{algebra::cblas_gemm(left.data(), right.data(), alpha)};
}

template <typename T>
LowRankTile<T>
gemm(const LowRankTile<T> &left, const LowRankTile<T> &right, double alpha) {
    auto mid = algebra::cblas_gemm(left.matrixR(), right.matrixL(), 1.0);
    if (left.rank() > right.rank()) {
        return LowRankTile<T>{algebra::cblas_gemm(left.matrixL(), mid, alpha),
                              right.matrixR()};
    } else {
        return LowRankTile<T>{alpha * left.matrixL(),
                              algebra::cblas_gemm(mid, right.matrixR(), 1.0)};
    }
}

} // namespace tile_ops

#endif // TileClusterChem_TILE_OPS_H
