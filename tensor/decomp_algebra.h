#pragma once

#include "tensor_fwd.h"
#include "../include/eigen.h"

namespace tcc {
namespace tensor {
namespace algebra {

struct DecomposedTensor {
    DecomposedTensor(double cut) : cut_(cut) {}
    DecomposedTensor(TARange const &range, TATensor &&left, TATensor &&right,
                     double cut)
        : product_range_(range),
          cut_(cut),
          L_(std::move(left)),
          R_(std::move(right)) {}

    TARange product_range_;
    double cut_;
    TATensor L_;
    TATensor R_;
};

// Just use Eigen for now, we can speed this up later.
// For now we will always decompose {1,{2,3}} or {1,2}
// eventually there will need to be a helper object to tell us how to partiton.
inline DecomposedTensor ColPivotedQr(TATensor const &T, double cut) {
    // Get range info
    auto const &range = T.range();
    auto const &extent = range.size();

    // Decide on matricization
    const auto dim = range.dim();
    auto M = 0ul; // Nrows
    auto N = 0ul; // Ncols
    if (dim == 2) {
        M = extent[0];
        N = extent[1];
    } else if (dim == 3) {
        M = extent[0];
        N = extent[1] * extent[2];
    } else {
        assert(false);
    }


    // Compute the QR of T
    Eigen::ColPivHouseholderQR<RowMatrixXd> qr(TiledArray::eigen_map(T, M, N));

    qr.setThreshold(cut);
    const auto rank = qr.rank();

    RowMatrixXd R = RowMatrixXd(qr.matrixR()
                          .topLeftCorner(rank, N)
                          .template triangularView<Eigen::Upper>())
               * qr.colsPermutation().transpose();

    RowMatrixXd L = RowMatrixXd(qr.householderQ()).leftCols(rank);

    // Form ranges L is always M by rank
    TARange Lrange{M, rank};
    // R is either rank, N for dim 2 or rank, extent[1], extent[2] for dim = 3
    TARange Rrange;
    if (dim == 2) {
        Rrange = TARange{rank, N};
    } else if (dim == 3) {
        Rrange = TARange{rank, extent[1], extent[2]};
    }

    return DecomposedTensor{T.range(), TATensor{std::move(Lrange), L.data()},
                            TATensor{std::move(Rrange), R.data()}, cut};
}

} // namespace algebra
} // namespace tensor
} // namespace tcc
