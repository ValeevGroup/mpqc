#pragma once

#include "tensor_fwd.h"
#include "../include/eigen.h"
#include <vector>
#include <tuple>

namespace tcc {
namespace tensor {
namespace algebra {

// Just use Eigen for now, we can speed this up later.
// For now we will always decompose {1,{2,3}, 4} where we assume that
// the centers correspond to {1,2,2,1} as far as r_1 and r_2 goes in the
// integral.
inline std::tuple<TARange, std::vector<TATensor>>
exchange_decomposition(TATensor const &t, double cut) {
    // Get range info
    auto const &range = t.range();
    auto const &extent = range.size();
    assert(range.dim() == 4); // Exchange is 4 centered

    std::tuple<TARange, std::vector<TATensor>> tensors;
    std::get<0>(tensors) = range;

    // Decompose first dimension
    auto rows = extent[0] * extent[1];
    auto cols = extent[2] * extent[3];
    auto map = TA::eigen_map(t, rows, cols);

    Eigen::ColPivHouseholderQR<RowMatrixXd> qr(map);
    qr.setThreshold(cut);
    const auto rank = qr.rank();

    // If not decomposable just return tensor and exit.
    if (rank > 0.5 * rows) {
        std::get<1>(tensors).emplace_back(t.clone());
        return tensors;
    }

    RowMatrixXd eigl = RowMatrixXd(qr.householderQ()).leftCols(rank);
    RowMatrixXd eigr = RowMatrixXd(qr.matrixR()
                                .topLeftCorner(rank, cols)
                                .template triangularView<Eig::Upper>())
                * qr.colsPermutation().transpose();

    { // decompose left dims
        auto lrows = extent[0];
        auto lcols = rank * extent[1];
        eigl.resize(lrows, lcols);
        qr.compute(eigl);
        auto rank_l = qr.rank();

        TARange l_range{lrows, rank_l};
        TATensor l(l_range);
        auto l_map = TA::eigen_map(l, lrows, rank_l);
        l_map = RowMatrixXd(qr.householderQ()).leftCols(rank_l);
        std::get<1>(tensors).push_back(l);

        TARange m1_range{rank_l, extent[1], rank};
        TATensor m1(m1_range);
        auto m1_map = TA::eigen_map(m1, rank_l, lcols);
        m1_map = RowMatrixXd(qr.matrixR()
                                 .topLeftCorner(rank_l, lcols)
                                 .template triangularView<Eig::Upper>())
                 * qr.colsPermutation().transpose();
        std::get<1>(tensors).push_back(m1);
    }
    { // Decompose right dims
        auto rrows = rank * extent[2];
        auto rcols = extent[3];
        eigr.resize(rrows, rcols);
        qr.compute(eigr);
        auto rank_r = qr.rank();

        TARange m2_range{rank, extent[2], rank_r};
        TATensor m2(m2_range);
        auto m2_map = TA::eigen_map(m2, rrows, rank_r);
        m2_map = RowMatrixXd(qr.householderQ()).leftCols(rank_r);
        std::get<1>(tensors).push_back(m2);

        TARange r_range{rank_r, extent[3]};
        TATensor r(r_range);
        auto r_map = TA::eigen_map(r, rank_r, rcols);
        r_map = RowMatrixXd(qr.matrixR()
                                 .topLeftCorner(rank_r, rcols)
                                 .template triangularView<Eig::Upper>())
                 * qr.colsPermutation().transpose();
        std::get<1>(tensors).push_back(r);
    }

    return tensors;
}

inline std::tuple<TARange, std::vector<TATensor>>
coulomb_decomposition(TATensor const &t, double cut) {
    // Get range info
    auto const &range = t.range();
    auto const &extent = range.size();
    assert(range.dim() == 4); // Exchange is 4 centered

    std::tuple<TARange, std::vector<TATensor>> tensors;
    std::get<0>(tensors) = range;

    // Decompose first dimension
    auto rows = extent[0] * extent[1];
    auto cols = extent[2] * extent[3];
    auto map = TA::eigen_map(t, rows, cols);

    Eigen::ColPivHouseholderQR<RowMatrixXd> qr(map);
    qr.setThreshold(cut);
    const auto rank = qr.rank();

    // If not decomposable just return tensor and exit.
    if (rank > 0.5 * rows) {
        std::get<1>(tensors).emplace_back(t.clone());
        return tensors;
    }

    // Form the left most tensor.
    TARange l_range{rows, rank};
    TATensor l(l_range);
    auto l_map = TA::eigen_map(l, rows, rank);
    l_map = RowMatrixXd(qr.householderQ()).leftCols(rank);
    std::get<1>(tensors).push_back(l);

    // form the right tensor which is the mirror of l
    TARange r_range(rank, cols);
    TATensor r(r_range);
    auto r_map = TA::eigen_map(r, rank, cols);
    r_map = RowMatrixXd(qr.matrixR()
                            .topLeftCorner(rank, cols)
                            .template triangularView<Eig::Upper>())
            * qr.colsPermutation().transpose();

    std::get<1>(tensors).push_back(r);

    return tensors;
}

} // namespace algebra
} // namespace tensor
} // namespace tcc
