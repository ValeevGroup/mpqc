#pragma once
#ifndef TCC_TENSOR_DECOMPOSEDTENSORALGEBRA_H
#define TCC_TENSOR_DECOMPOSEDTENSORALGEBRA_H

#include "decomposed_tensor.h"
#include "../include/eigen.h"
#include "../common/typedefs.h"
#include <madness/tensor/clapack.h>

namespace tcc {
namespace tensor {
namespace algebra {

inline std::size_t qr_rank(double const *data, std::size_t rows,
                           std::size_t cols, double threshhold) {
    // TODO
    assert(false);
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
    const int rows = extent[0];
    const int cols = extent[1] * extent[2];
    auto full_rank = std::min(rows, cols);

    Eig::VextorXi J = Eig::VectorXi::Zero(cols);
    double Tau[full_rank];
    double work;
    int LWORK = -1; // Ask for space computation
    int INFO;
    int LDA = rows;

    // Have to copy data since dgeqp3_ eats input
    const auto tensor_size = rows * cols;
    unique_ptr<double[]> t_copy{new double[tensor_size]};
    const auto data = t.tensor(0).data();
    std::copy(data, data + tensor_size, t_copy);

    dgeqp3_(&rows, &cols, t_copy.get(), &LDA, J.data(), Tau, &work, &LWORK,
            &INFO);
    LWORK = work;
    std::unique_ptr<double[]> W{new double[LWORK]};
    dgeqp3_(&rows, &cols, t_copy.get(), &LDA, J.data(), Tau, W.get(), &LWORK,
            &INFO);

    // TODO finish
    auto rank = qr_rank(t_copy.get(), rows, cols, t.cut());
}

} // namespace algebra
} // namespace tensor
} // namespace tcc

#endif // TCC_TENSOR_DECOMPOSEDTENSORALGEBRA_H
