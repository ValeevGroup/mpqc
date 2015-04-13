#pragma once
#ifndef TCC_TENSOR_DECOMPOSEDTENSORALGEBRA_H
#define TCC_TENSOR_DECOMPOSEDTENSORALGEBRA_H

#include "decomposed_tensor.h"
#include "../include/eigen.h"
#include "../common/namespaces.h"
#include "../common/typedefs.h"
#include <madness/tensor/clapack.h>

#include <memory>

namespace tcc {
namespace tensor {
namespace algebra {

// calculate the qr_rank of a matrix.
inline int qr_rank(double const *data, int rows, int cols, double threshhold) {
    const auto full_rank = std::min(cols, rows);
    auto out_rank = full_rank;

    auto M = Eig::Map<const Eig::MatrixXd>(data, rows, cols);

    auto squared_sum = 0.0;
    for (auto i = (full_rank - 1); i >= 0; --i) { // rows of R

        for (auto j = (cols - 1); j >= i; --j) { // cols of R
            squared_sum += M(i, j) * M(i, j);
        }

        if (std::sqrt(squared_sum) >= threshhold) {
            return out_rank;
        }

        --out_rank; // Decriment rank and go to next row.
    }

    return out_rank;
}

// Returns true if input is low rank.
bool col_pivoted_qr(TA::Tensor<double> const &in, TA::Tensor<double> &L,
                    TA::Tensor<double> &R, double thresh) {
    auto const &extent = in.range().size();

    // for now assume rank 3 and assume (X, {ab}) or (r, {ab}) decomp
    // swapping cols and rows because dgeqp3 is a column major routine.
    int cols = extent[0];
    int rows = extent[1] * extent[2];
    auto full_rank = std::min(rows, cols);

    Eig::VectorXi J = Eig::VectorXi::Zero(cols);
    double Tau[full_rank]; // Dyanmic stack array ala C99 may not be portable
    double work;
    int LWORK = -1; // Ask for space computation
    int INFO;
    int LDA = rows;

    // Have to copy data since dgeqp3_ eats input
    const auto tensor_size = rows * cols;
    std::unique_ptr<double[]> t_copy{new double[tensor_size]};
    const auto data = in.data();
    std::copy(data, data + tensor_size, t_copy.get());

    // Call routine
    dgeqp3_(&rows, &cols, t_copy.get(), &LDA, J.data(), Tau, &work, &LWORK,
            &INFO);
    LWORK = work;
    std::unique_ptr<double[]> W{new double[LWORK]};
    dgeqp3_(&rows, &cols, t_copy.get(), &LDA, J.data(), Tau, W.get(), &LWORK,
            &INFO);

    // Determine rank and if decomposing is worth it.
    auto rank = qr_rank(t_copy.get(), rows, cols, thresh);
    if (rank > 0.5 * double(full_rank)) {
        return false;
    }

    // LAPACK assumes 1 based indexing, but we need zero.
    std::for_each(J.data(), J.data() + J.size(), [](int &val) { --val; });
    Eigen::PermutationWrapper<Eigen::VectorXi> P(J);

    // Create L tensor
    TA::Range l_range{extent[0], static_cast<std::size_t>(rank)};
    L = TA::Tensor<double>(std::move(l_range));

    // Eigen map the input
    auto A = Eig::Map<Eig::MatrixXd>(t_copy.get(), rows, cols);

    // Assign into l_tensor
    auto L_map = Eig::Map<Eig::MatrixXd>(L.data(), rank, cols);
    L_map = Eig::MatrixXd(A.topLeftCorner(rank, cols)
                                .template triangularView<Eigen::Upper>())
            * P.transpose();

    // compute q .
    dorgqr_(&rows, &rank, &rank, t_copy.get(), &rows, Tau, W.get(), &LWORK,
            &INFO);

    // From Q goes to R because of column major transpose issues.
    TA::Range q_range{static_cast<std::size_t>(rank), extent[1], extent[2]};
    R = TA::Tensor<double>(std::move(q_range));
    auto R_map = Eig::Map<Eig::MatrixXd>(R.data(), rows, rank);
    R_map = A.leftCols(rank);
    return true;
}

/// Returns an empty DecomposedTensor if the rank does not decrease
inline DecomposedTensor<double>
recompress_right(DecomposedTensor<double> const &t) {
    assert(t.ndecomp() >= 2);
    TA::Tensor<double> Rl, Rr;
    if (col_pivoted_qr(t.tensor(1), Rl, Rr, t.cut())) {
        // Make new left tensors
        constexpr auto NoT = madness::cblas::CBLAS_TRANSPOSE::NoTrans;
        auto gh = TA::math::GemmHelper(NoT, NoT, 2, 2, 2);
        auto new_left = t.tensor(0).gemm(Rl, 1.0, gh);

        // Use decomposed new left and right tensors.
        return DecomposedTensor<double>(t.cut(), std::move(new_left),
                                        std::move(Rr));

    } else {
        return DecomposedTensor<double>{};
    }
}

/// Returns an empty DecomposedTensor if the compression rank was to large.
inline DecomposedTensor<double>
two_way_decomposition(DecomposedTensor<double> const &t) {
    if (t.ndecomp() >= 2) {
        return recompress_right(t);
    }

    TA::Tensor<double> L, R;
    if (col_pivoted_qr(t.tensor(0), L, R, t.cut())) {
        return DecomposedTensor<double>(t.cut(), std::move(L), std::move(R));
    } else {
        return DecomposedTensor<double>{};
    }
}

TA::Tensor<double> combine(DecomposedTensor<double> const &t) {
    if (t.empty()) {
        return TA::Tensor<double>{};
    }

    const auto NoT = madness::cblas::CBLAS_TRANSPOSE::NoTrans;
    auto const &tensors = t.tensors();
    auto out = tensors[0].clone();
    for (auto i = 1ul; i < tensors.size(); ++i) {
        auto l_dim = out.range().dim();
        auto r_dim = tensors[i].range().dim();
        auto result_dim = l_dim + r_dim - 2; // Only one contraction index.
        auto gh = TA::math::GemmHelper(NoT, NoT, result_dim, l_dim, r_dim);
        out = out.gemm(tensors[i], 1.0, gh);
    }

    return out;
}

} // namespace algebra
} // namespace tensor
} // namespace tcc

#endif // TCC_TENSOR_DECOMPOSEDTENSORALGEBRA_H
