#pragma once
#ifndef TCC_TENSOR_DECOMPOSEDTENSORALGEBRA_H
#define TCC_TENSOR_DECOMPOSEDTENSORALGEBRA_H

#include "decomposed_tensor.h"
#include "../include/eigen.h"
#include "../common/namespaces.h"
#include "../common/typedefs.h"
#include <madness/tensor/clapack.h>

#include <memory>
#include <limits>
#include <chrono>

extern "C" void
sgesdd_(const char *jobz, integer *m, integer *n, real4 *a, integer *lda,
        real4 *s, real4 *u, integer *ldu, real4 *vt, integer *ldvt, real4 *work,
        integer *lwork, integer *iwork, integer *info);

extern "C" void
dgesdd_(const char *jobz, integer *m, integer *n, real8 *a, integer *lda,
        real8 *s, real8 *u, integer *ldu, real8 *vt, integer *ldvt, real8 *work,
        integer *lwork, integer *iwork, integer *info);

extern "C" void
dgelqf_(integer *m, integer *n, real8 *a, integer *lda, real8 *tau, real8 *work,
        integer *lwork, integer *infoOUT);

extern "C" void
dorglq_(integer *m, integer *n, integer *k, real8 *a, integer *lda, real8 *tau,
        real8 *work, integer *lwork, integer *info);

namespace tcc {
namespace tensor {
namespace algebra {

// Compute the column pivoted qr decomposition into data, will modify input
// pointers data and J
int col_pivoted_qr(double *data, double *Tau, int rows, int cols, int *J) {
    double work_dummy;
    int LWORK = -1; // Ask for space computation
    int INFO;
    int LDA = rows;

    // Call routine
    dgeqp3_(&rows, &cols, data, &LDA, J, Tau, &work_dummy, &LWORK, &INFO);
    LWORK = work_dummy;
    std::unique_ptr<double[]> W{new double[LWORK]};
    dgeqp3_(&rows, &cols, data, &LDA, J, Tau, W.get(), &LWORK, &INFO);
    return INFO;
}

int non_pivoted_qr(double *data, double *Tau, int rows, int cols) {
    double work_dummy;
    int LWORK = -1; // Ask for space computation
    int INFO;
    int LDA = rows;

    // Call routine
    dgeqrf_(&rows, &cols, data, &LDA, Tau, &work_dummy, &LWORK, &INFO);
    LWORK = work_dummy;
    std::unique_ptr<double[]> W{new double[LWORK]};
    dgeqrf_(&rows, &cols, data, &LDA, Tau, W.get(), &LWORK, &INFO);
    return INFO;
}

int non_pivoted_lq(double *data, double *Tau, int rows, int cols) {
    double work_dummy;
    int LWORK = -1; // Ask for space computation
    int INFO;
    int LDA = rows;

    // Call routine
    dgelqf_(&rows, &cols, data, &LDA, Tau, &work_dummy, &LWORK, &INFO);
    LWORK = work_dummy;
    std::unique_ptr<double[]> W{new double[LWORK]};
    dgelqf_(&rows, &cols, data, &LDA, Tau, W.get(), &LWORK, &INFO);
    return INFO;
}

int svd(double *data, double *s, double *u, double *vt, int rows, int cols) {
    double work_dummy;
    int LWORK = -1; // Ask for space computation
    int INFO;
    int LDA = rows;
    int LDU = (rows < cols) ? rows : 1;
    int LDVT = (rows >= cols) ? cols : 1;
    const char O = 'O';
    std::unique_ptr<int[]> iwork{new int[8 * std::min(rows, cols)]};

    // Call routine
    dgesdd_(&O, &rows, &cols, data, &LDA, s, u, &LDU, vt, &LDVT, &work_dummy,
            &LWORK, iwork.get(), &INFO);
    LWORK = work_dummy;
    std::unique_ptr<double[]> W{new double[LWORK]};
    dgesdd_(&O, &rows, &cols, data, &LDA, s, u, &LDU, vt, &LDVT, W.get(),
            &LWORK, iwork.get(), &INFO);

    return INFO;
}

int form_q(double *data, double *Tau, int rows, int rank) {
    double work_dummy = 0.0;
    int LWORK = -1;
    int INFO;
    dorgqr_(&rows, &rank, &rank, data, &rows, Tau, &work_dummy, &LWORK, &INFO);
    LWORK = work_dummy;
    std::unique_ptr<double[]> work{new double[LWORK]};

    dorgqr_(&rows, &rank, &rank, data, &rows, Tau, work.get(), &LWORK, &INFO);

    return INFO;
}

int form_q_from_lq(double *data, double *Tau, int cols, int rows, int rank) {
    double work_dummy = 0.0;
    int LWORK = -1;
    int INFO;
    dorglq_(&rank, &cols, &rank, data, &rows, Tau, &work_dummy, &LWORK, &INFO);
    LWORK = work_dummy;
    std::unique_ptr<double[]> work{new double[LWORK]};

    dorglq_(&rank, &cols, &rank, data, &rows, Tau, work.get(), &LWORK, &INFO);

    return INFO;
}

inline int svd_rank(double const *s, int N, double thresh) {
    auto rank = 0;
    for (auto i = 0; i < N; ++i) {
        if (s[i] >= thresh) {
            ++rank;
        } else {
            break;
        }
    }

    return std::max(rank, 1);
}


// calculate the qr_rank of a matrix.
inline std::size_t qr_rank(double const *data, std::size_t rows,
                           std::size_t cols, double threshhold) {
    const auto full_rank = std::min(cols, rows);
    auto out_rank = full_rank;

    auto M = Eig::Map<const Eig::MatrixXd>(data, rows, cols);

    auto squared_sum = 0.0;
    for (int i = (full_rank - 1); i >= 0; --i) { // rows of R
        for (int j = (cols - 1); j >= i; --j) {  // cols of R
            squared_sum += M(i, j) * M(i, j);
        }

        if (std::sqrt(squared_sum) >= threshhold) {
            return out_rank;
        }

        --out_rank; // Decriment rank and go to next row.
    }

    return std::max(out_rank, 1ul);
}

// Returns true if input is low rank.
bool full_rank_decompose(TA::Tensor<double> const &in, TA::Tensor<double> &L,
                         TA::Tensor<double> &R, double thresh,
                         std::size_t max_out_rank
                         = std::numeric_limits<std::size_t>::max()) {
    auto const &extent = in.range().size();

    int cols = extent[0];
    int rows = 1;
    for (auto i = 1ul; i < extent.size(); ++i) {
        rows *= extent[i];
    }
    const auto size = cols * rows;
    auto full_rank = std::min(rows, cols);

    // Will hold the column pivot information and reflectors
    Eig::VectorXi J = Eig::VectorXi::Zero(cols);
    std::unique_ptr<double[]> Tau{new double[full_rank]};

    std::unique_ptr<double[]> in_data{new double[size]};
    std::copy(in.data(), in.data() + size, in_data.get());

    // Do initial qr routine
    auto qr_err
          = col_pivoted_qr(in_data.get(), Tau.get(), rows, cols, J.data());
    if (0 != qr_err) {
        std::cout << "Something went wrong with computing qr.\n";
        throw;
    }

    // Determine rank and if decomposing is worth it.
    auto rank = qr_rank(in_data.get(), rows, cols, thresh);
    if (max_out_rank == std::numeric_limits<std::size_t>::max()) {
        max_out_rank = full_rank / 2;
    }
    if (rank >= max_out_rank) {
        return false;
    }

    // LAPACK assumes 1 based indexing, but we need zero.
    std::for_each(J.data(), J.data() + J.size(), [](int &val) { --val; });
    Eigen::PermutationWrapper<Eigen::VectorXi> P(J);

    // Create L tensor
    TA::Range l_range{extent[0], static_cast<std::size_t>(rank)};
    L = TA::Tensor<double>(std::move(l_range));

    // Eigen map the input
    auto A = Eig::Map<Eig::MatrixXd>(in_data.get(), rows, cols);

    // Assign into l_tensor
    auto L_map = Eig::Map<Eig::MatrixXd>(L.data(), rank, cols);
    L_map = Eig::MatrixXd(A.topLeftCorner(rank, cols)
                                .template triangularView<Eigen::Upper>())
            * P.transpose();

    auto q_err = form_q(in_data.get(), Tau.get(), rows, rank);
    if (0 != q_err) {
        std::cout << "Something went wrong with forming q.\n";
        throw;
    }

    // From Q goes to R because of column major transpose issues.
    if (extent.size() == 3) {
        TA::Range q_range{static_cast<std::size_t>(rank), extent[1], extent[2]};
        R = TA::Tensor<double>(std::move(q_range));
    } else if (extent.size() == 2) {
        TA::Range q_range{static_cast<std::size_t>(rank), extent[1]};
        R = TA::Tensor<double>(std::move(q_range));
    } else {
        assert(false);
    }

    auto R_map = Eig::Map<Eig::MatrixXd>(R.data(), rows, rank);
    R_map = A.leftCols(rank);
    return true;
}

// Returns true if input is low rank.
void ta_tensor_col_pivoted_qr(TA::Tensor<double> &in, TA::Tensor<double> &L,
                              TA::Tensor<double> &R, double thresh) {
    auto const &extent = in.range().size();

    int cols = extent[0];
    int rows = 1;
    for (auto i = 1ul; i < extent.size(); ++i) {
        rows *= extent[i];
    }
    auto full_rank = std::min(rows, cols);

    // Will hold the column pivot information and reflectors
    Eig::VectorXi J = Eig::VectorXi::Zero(cols);
    std::unique_ptr<double[]> Tau{new double[full_rank]};

    // Do initial qr routine
    auto qr_err = col_pivoted_qr(in.data(), Tau.get(), rows, cols, J.data());
    if (0 != qr_err) {
        std::cout << "Something went wrong with computing qr.\n";
        throw;
    }

    // Determine rank and if decomposing is worth it.
    auto rank = qr_rank(in.data(), rows, cols, thresh);

    // LAPACK assumes 1 based indexing, but we need zero.
    std::for_each(J.data(), J.data() + J.size(), [](int &val) { --val; });
    Eigen::PermutationWrapper<Eigen::VectorXi> P(J);

    // Create L tensor
    TA::Range l_range{extent[0], static_cast<std::size_t>(rank)};
    L = TA::Tensor<double>(std::move(l_range));

    // Eigen map the input
    auto A = Eig::Map<Eig::MatrixXd>(in.data(), rows, cols);

    // Assign into l_tensor
    auto L_map = Eig::Map<Eig::MatrixXd>(L.data(), rank, cols);
    L_map = Eig::MatrixXd(A.topLeftCorner(rank, cols)
                                .template triangularView<Eigen::Upper>())
            * P.transpose();

    auto q_err = form_q(in.data(), Tau.get(), rows, rank);

    if (0 != q_err) {
        std::cout << "Something went wrong with forming q.\n";
        throw;
    }

    // Form Q goes to R because of column major transpose issues.
    if (extent.size() == 3) {
        TA::Range q_range{static_cast<std::size_t>(rank), extent[1], extent[2]};
        R = TA::Tensor<double>(std::move(q_range));
    } else if (extent.size() == 2) {
        TA::Range q_range{static_cast<std::size_t>(rank), extent[1]};
        R = TA::Tensor<double>(std::move(q_range));
    } else {
        assert(false);
    }

    auto R_map = Eig::Map<Eig::MatrixXd>(R.data(), rows, rank);
    R_map = A.leftCols(rank);
}

// eats in data and outputs L and R tensors.
inline void ta_tensor_qr(TA::Tensor<double> &in, TA::Tensor<double> &L,
                         TA::Tensor<double> &R) {

    auto const &extent = in.range().size();

    // Reverse map
    // We will always extract the first dimension as cols
    int cols = extent[0];

    // The rest will be smashed into rows
    const auto ndims = extent.size();
    int rows = 1;
    for (auto i = 1ul; i < ndims; ++i) {
        rows *= extent[i];
    }

    auto full_rank = std::min(rows, cols);

    // Will hold the reflectors
    std::unique_ptr<double[]> Tau{new double[full_rank]};

    // Lets start by not copying data
    /* const auto size = cols * rows; */
    /* std::unique_ptr<double[]> in_data{new double[size]}; */
    /* std::copy(in.data(), in.data() + size, in_data.get()); */

    auto qr_err = non_pivoted_qr(in.data(), Tau.get(), rows, cols);

    if (0 != qr_err) {
        std::cout << "Something went wrong with computing qr.\n";
        throw;
    }

    TA::Range l_range{cols, full_rank};
    L = TA::Tensor<double>(std::move(l_range));

    // Eigen map the input
    auto A = Eig::Map<Eig::MatrixXd>(in.data(), rows, cols);

    // Assign into l_tensor
    auto L_map = Eig::Map<Eig::MatrixXd>(L.data(), full_rank, cols);
    L_map = A.topLeftCorner(full_rank, cols)
                  .template triangularView<Eigen::Upper>();

    auto q_err = form_q(in.data(), Tau.get(), rows, full_rank);
    if (0 != q_err) {
        std::cout << "Something went wrong with forming q.\n";
        throw;
    }

    if (ndims == 2) {
        TA::Range q_range{full_rank, extent[1]};
        R = TA::Tensor<double>(std::move(q_range));
    } else if (ndims == 3) {
        TA::Range q_range{full_rank, extent[1], extent[2]};
        R = TA::Tensor<double>(std::move(q_range));
    } else {
        assert(false);
    }


    auto R_map = Eig::Map<Eig::MatrixXd>(R.data(), rows, full_rank);
    R_map = A.leftCols(full_rank);
}

// eats in data and outputs L and R tensors.
inline void ta_tensor_lq(TA::Tensor<double> &in, TA::Tensor<double> &L,
                         TA::Tensor<double> &R) {

    auto const &extent = in.range().size();

    // Reverse map
    // We will always extract the first dimension as cols
    int cols = extent[0];

    // The rest will be smashed into rows
    const auto ndims = extent.size();
    int rows = 1;
    for (auto i = 1ul; i < ndims; ++i) {
        rows *= extent[i];
    }

    auto full_rank = std::min(rows, cols);

    // Will hold the reflectors
    std::unique_ptr<double[]> Tau{new double[full_rank]};
    auto qr_err = non_pivoted_lq(in.data(), Tau.get(), rows, cols);
    if (0 != qr_err) {
        std::cout << "Something went wrong with computing qr.\n";
        throw;
    }

    if (ndims == 2) {
        TA::Range q_range{full_rank, extent[1]};
        R = TA::Tensor<double>(std::move(q_range));
    } else if (ndims == 3) {
        TA::Range q_range{full_rank, extent[1], extent[2]};
        R = TA::Tensor<double>(std::move(q_range));
    } else {
        assert(false);
    }

    auto A = Eig::Map<Eig::MatrixXd>(in.data(), rows, cols);

    auto R_map = Eig::Map<Eig::MatrixXd>(R.data(), rows, full_rank);
    R_map = A.bottomLeftCorner(rows, full_rank)
                  .template triangularView<Eigen::Lower>();

    auto q_err = form_q_from_lq(in.data(), Tau.get(), cols, rows, full_rank);
    if (0 != q_err) {
        std::cout << "Something went wrong with forming q.\n";
        throw;
    }

    TA::Range l_range{cols, full_rank};
    L = TA::Tensor<double>(std::move(l_range));

    Eig::Map<Eig::MatrixXd> L_map(L.data(), full_rank, cols);
    L_map = A.topRows(full_rank);
}

// Sacrifice input data
inline void ta_tensor_svd(TA::Tensor<double> &in, TA::Tensor<double> &L,
                          TA::Tensor<double> &R, double thresh) {
    auto const &extent = in.range().size();
    auto ndims = extent.size();
    assert(ndims == 2);

    int cols = extent[0];
    int rows = 1;
    for (auto i = 1ul; i < ndims; ++i) {
        rows *= extent[i];
    }

    Eigen::VectorXd s(std::max(rows, cols));
    std::unique_ptr<double[]> u{new double[(rows >= cols) ? 1 : rows * rows]};
    std::unique_ptr<double[]> vt{new double[(rows >= cols) ? cols * cols : 1]};

    svd(in.data(), s.data(), u.get(), vt.get(), rows, cols);

    auto rank = svd_rank(s.data(), std::max(rows, cols), thresh);

    TA::Range l_range{extent[0], static_cast<std::size_t>(rank)};
    TA::Range r_range{static_cast<std::size_t>(rank), extent[1]};
    L = TA::Tensor<double>(std::move(l_range));
    R = TA::Tensor<double>(std::move(r_range));

    if (rows >= cols) { // U is written into in
        auto Ldata = L.data();
        auto vt_ptr = vt.get();
        for (auto i = 0ul; i < extent[0]; ++i) {
            std::copy(vt_ptr, vt_ptr + rank, Ldata);
            vt_ptr += extent[0];
            Ldata += rank;
        }
        std::copy(in.data(), in.data() + R.range().volume(), R.data());

        // Roll in diagonal
        auto R_map = TA::eigen_map(R, rank, extent[1]);
        R_map = s.head(rank).asDiagonal() * R_map;
    } else { // vt is written into in
        auto Ldata = L.data();
        auto vt_ptr = in.data();
        for (auto i = 0ul; i < extent[0]; ++i) {
            std::copy(vt_ptr, vt_ptr + rank, Ldata);
            vt_ptr += extent[0];
            Ldata += rank;
        }
        std::copy(u.get(), u.get() + R.range().volume(), R.data());

        // Roll in diagonal
        auto R_map = TA::eigen_map(R, rank, extent[1]);
        R_map = s.head(rank).asDiagonal() * R_map;
    }
}

/// Currently modifies input data regardless could cause some loss of accuracy.
inline void recompress(DecomposedTensor<double> &t) {
    assert(t.ndecomp() >= 2);

    // Where matrix M = ST;
    TA::Tensor<double> Ls, Rs, Lt, Rt;
    ta_tensor_lq(t.tensor(0), Ls, Rs); // LQ Left to end up with QR;
    ta_tensor_qr(t.tensor(1), Lt, Rt); // QR right to end up with LQ;

    // Form a M matrix
    constexpr auto NoT = madness::cblas::CBLAS_TRANSPOSE::NoTrans;
    const auto gh = TA::math::GemmHelper(NoT, NoT, 2, 2, 2);
    auto M = Rs.gemm(Lt, 1.0, gh);

    // want to always do the full decomp so make input
    // max rank larger than rank of M.
    TA::Tensor<double> Lm, Rm;
    // ta_tensor_svd(M, Lm, Rm, t.cut());
    ta_tensor_col_pivoted_qr(M, Lm, Rm, t.cut());
    auto newL = Ls.gemm(Lm, 1.0, gh);

    const auto gh_r = TA::math::GemmHelper(NoT, NoT, 3, 2, 3);
    auto newR = Rm.gemm(Rt, 1.0, gh_r);
    t = DecomposedTensor<double>(t.cut(), std::move(newL), std::move(newR));
}

/// Returns an empty DecomposedTensor if the compression rank was to large.
inline DecomposedTensor<double>
two_way_decomposition(DecomposedTensor<double> const &t) {
    if (t.ndecomp() >= 2) {
        assert(false);
    }

    auto const &extent = t.tensor(0).range().size();
    const auto rows = extent[0];
    const auto cols = extent[1] * extent[2];
    const auto max_out_rank = std::min(rows, cols) / 2;
    TA::Tensor<double> L, R;
    if (full_rank_decompose(t.tensor(0), L, R, t.cut(), max_out_rank)) {
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
