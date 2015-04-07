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

/// Returns an empty DecomposedTensor if the rank does note decrease
inline DecomposedTensor<double>
recompress_left(DecomposedTensor<double> const &t) {
    assert(false);
}

/// Returns an empty DecomposedTensor if the compression rank was to large.
inline DecomposedTensor<double>
two_way_decomposition(DecomposedTensor<double> const &t) {
    if (t.ndecomp() >= 2) {
        return recompress_left(t);
    }

    auto const &extent = t.range(0).size();

    // for now assume rank 3 and assume (X, {ab}) decomp
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
    const auto data = t.tensor(0).data();
    std::copy(data, data + tensor_size, t_copy.get());

    // Call routine
    dgeqp3_(&rows, &cols, t_copy.get(), &LDA, J.data(), Tau, &work, &LWORK,
            &INFO);
    LWORK = work;
    std::unique_ptr<double[]> W{new double[LWORK]};
    dgeqp3_(&rows, &cols, t_copy.get(), &LDA, J.data(), Tau, W.get(), &LWORK,
            &INFO);

    // Determine rank and if decomposing is worth it.
    auto rank = qr_rank(t_copy.get(), rows, cols, t.cut());
    if (rank > 0.5 * double(full_rank)) {
        return DecomposedTensor<double>{};
    }

    // LAPACK assumes 1 based indexing, but we need zero.
    std::for_each(J.data(), J.data() + J.size(), [](int &val) { --val; });
    Eigen::PermutationWrapper<Eigen::VectorXi> P(J);

    // Create L tensor
    TA::Range l_range{extent[0], static_cast<std::size_t>(rank)};
    TA::Tensor<double> l_tensor(std::move(l_range));

    // Eigen map the input
    auto A = Eig::Map<Eig::MatrixXd>(t_copy.get(), rows, cols);

    // Assign into l_tensor
    auto L = Eig::Map<Eig::MatrixXd>(l_tensor.data(), rank, cols);
    L = Eigen::MatrixXd(A.topLeftCorner(rank, cols)
                              .template triangularView<Eigen::Upper>())
        * P.transpose();

    // compute q .
    dorgqr_(&rows, &rank, &rank, t_copy.get(), &rows, Tau, W.get(), &LWORK,
            &INFO);

    // From Q
    TA::Range q_range {static_cast<std::size_t>(rank), extent[1], extent[2]};
    TA::Tensor<double> q_tensor(std::move(q_range));
    auto Q = Eig::Map<Eig::MatrixXd>(q_tensor.data(), rows, rank);
    Q = A.leftCols(rank);

    // Return tensors in what looks like the wrong order because the output
    // of column based QR on row major input will be LQ.
    return DecomposedTensor<double>(t.cut(), std::move(l_tensor),
                                    std::move(q_tensor));
}

} // namespace algebra
} // namespace tensor
} // namespace tcc

#endif // TCC_TENSOR_DECOMPOSEDTENSORALGEBRA_H
