/*
 * File to define linear algebra functions
 */
#ifndef Lin_Algebra_H
#define Lin_Algebra_H

#include "../include/eigen.h"
#include <TiledArray/madness.h>
#include <madness/tensor/clapack.h>

namespace algebra {

inline namespace eigen_version {
template <typename T>
RowMatrix<T> inline cblas_gemm(const RowMatrix<T> &A, const RowMatrix<T> &B,
                               double alpha) {
    return alpha * A * B;
}

template <typename T>
void inline cblas_gemm_inplace(const RowMatrix<T> &A, const RowMatrix<T> &B,
                               RowMatrix<T> &C, double alpha,
                               double beta = 1.0) {
    C = alpha * A * B + beta * C;
}

template <typename T>
bool inline Decompose_Matrix(RowMatrix<T> input, RowMatrix<T> &L,
                             RowMatrix<T> &R, double cut) {
    Eigen::ColPivHouseholderQR<RowMatrix<T>> qr(input);
    qr.setThreshold(cut);
    auto rank = qr.rank();

    auto full_rank = std::min(input.cols(), input.rows());
    if (rank >= double(full_rank) / 2.0) {
        return true;
    }

    R = RowMatrix<T>(qr.matrixR()
                         .topLeftCorner(rank, input.cols())
                         .template triangularView<Eigen::Upper>())
        * qr.colsPermutation().transpose();
    L = RowMatrix<T>(qr.householderQ()).leftCols(rank);

    return false;
}

template <typename T>
void inline ColPivotedQr(RowMatrix<T> input, RowMatrix<T> &L, RowMatrix<T> &R,
                         double cut) {
    Eigen::ColPivHouseholderQR<
        typename std::remove_reference<decltype(input)>::type> qr(input);

    qr.setThreshold(cut);
    auto rank = qr.rank();

    R = RowMatrix<T>(qr.matrixR()
                         .topLeftCorner(rank, input.cols())
                         .template triangularView<Eigen::Upper>())
        * qr.colsPermutation().transpose();

    L = RowMatrix<T>(qr.householderQ()).leftCols(rank);
}

template <typename T>
void inline CompressLeft(RowMatrix<T> &L, RowMatrix<T> &R, double cut) {

    Eigen::ColPivHouseholderQR<
        typename std::remove_reference<decltype(L)>::type> qr(L);

    qr.setThreshold(cut);
    auto rank = qr.rank();

    if (rank >= L.cols()) {
        return;
    }

    R = RowMatrix<T>(qr.matrixR()
                         .topLeftCorner(rank, qr.matrixQR().cols())
                         .template triangularView<Eigen::Upper>())
        * qr.colsPermutation().transpose() * R;

    L = RowMatrix<T>(qr.householderQ()).leftCols(rank);
}

template <typename T>
void inline CompressRight(RowMatrix<T> &L, RowMatrix<T> &R, double cut) {
    Eigen::ColPivHouseholderQR<
        typename std::remove_reference<decltype(R)>::type> qr(R);

    qr.setThreshold(cut);
    auto rank = qr.rank();

    if (rank >= R.rows()) {
        return;
    }

    R = RowMatrix<T>(qr.matrixR()
                         .topLeftCorner(rank, qr.matrixQR().cols())
                         .template triangularView<Eigen::Upper>())
        * qr.colsPermutation().transpose();

    L = L * RowMatrix<T>(qr.householderQ()).leftCols(rank);
}

} // namespace eigen_version

namespace lapack {
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> inline cblas_gemm(
    const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &A,
    const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &B,
    double alpha) {

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> C(A.rows(), B.cols());
    const int K = A.cols();
    const int M = A.rows();
    const int N = B.cols();
    if (K == 0 || M == 0 || N == 0) { // If zero leave
        for (auto i = 0; i < C.size(); ++i) {
            *(C.data() + i) = 0;
        }
        return C;
    }
    const int LDA = M, LDB = K, LDC = M;
    madness::cblas::gemm(madness::cblas::CBLAS_TRANSPOSE::NoTrans,
                         madness::cblas::CBLAS_TRANSPOSE::NoTrans, M, N, K,
                         alpha, A.data(), LDA, B.data(), LDB, 0.0, C.data(),
                         LDC);
    return C;
}

template <typename T>
RowMatrix<T> inline cblas_gemm(const RowMatrix<T> &A, const RowMatrix<T> &B,
                               double alpha) {
    return eigen_version::cblas_gemm(A, B, alpha);
}

void inline cblas_gemm_inplace(
    const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &A,
    const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &B,
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &C, double alpha,
    double beta = 1.0) {

    const int K = A.cols();
    const int M = A.rows();
    const int N = B.cols();
    if (K == 0 || M == 0 || N == 0) { // If zero leave
        C = beta * C;
        return;
    }
    const int LDA = M, LDB = K, LDC = M;
    madness::cblas::gemm(madness::cblas::CBLAS_TRANSPOSE::NoTrans,
                         madness::cblas::CBLAS_TRANSPOSE::NoTrans, M, N, K,
                         alpha, A.data(), LDA, B.data(), LDB, beta, C.data(),
                         LDC);
}

void inline cblas_gemm_inplace(
    const Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> &A,
    const Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> &B,
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> &C, double alpha,
    double beta = 1.0) {

    const int K = A.cols();
    const int M = C.rows();
    const int N = C.cols();
    if (K == 0 || M == 0 || N == 0) { // If zero leave
        C = beta * C;
        return;
    }
    const int LDA = M, LDB = K, LDC = M;
    madness::cblas::gemm(madness::cblas::CBLAS_TRANSPOSE::NoTrans,
                         madness::cblas::CBLAS_TRANSPOSE::NoTrans, M, N, K,
                         alpha, A.data(), LDA, B.data(), LDB, beta, C.data(),
                         LDC);
}

template <typename T>
void inline cblas_gemm_inplace(const RowMatrix<T> &A, const RowMatrix<T> &B,
                               RowMatrix<T> &C, double alpha,
                               double beta = 1.0) {

    eigen_version::cblas_gemm_inplace(A, B, C, alpha, beta);
}

// Only works on the output of dqeqp3
int inline qr_rank(
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> const &M,
    double thresh) {

    const auto full_rank = std::min(M.cols(), M.rows());
    auto out_rank = full_rank;

    auto squared_sum = 0.0;
    for (auto i = (full_rank - 1); i >= 0; --i) { // rows of R

        for (auto j = (M.cols() - 1); j >= i; --j) { // cols of R
            squared_sum += M(i, j) * M(i, j);
        }

        if (std::sqrt(squared_sum) >= thresh) {
            return out_rank;
        }

        --out_rank; // Decriment rank and go to next row.
    }

    return out_rank;
}

bool inline Decompose_Matrix(
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> input,
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &L,
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &R, double cut) {
    assert(input.size() >= 0);
    int M = input.rows();
    int N = input.cols();
    auto full_rank = std::min(M, N);
    Eigen::VectorXi J = Eigen::VectorXi::Zero(N);
    double Tau[full_rank];
    double work;
    int LWORK = -1; // Ask LAPACK how much space we need.
    int INFO;
    int LDA = M;

    dgeqp3_(&M, &N, input.data(), &LDA, J.data(), Tau, &work, &LWORK, &INFO);
    LWORK = work;
    std::unique_ptr<double[]> W{new double[LWORK]};
    dgeqp3_(&M, &N, input.data(), &LDA, J.data(), Tau, W.get(), &LWORK, &INFO);

    const double thresh = cut;
    auto rank = qr_rank(input, thresh);

    if (rank > 0.5 * double(full_rank)) {
        return true; // Input is full rank
    }

    // LAPACK assumes 1 based indexing, but we need zero.
    std::for_each(J.data(), J.data() + J.size(), [](int &val) { --val; });
    Eigen::PermutationWrapper<Eigen::VectorXi> P(J);
    R = Eigen::MatrixXd(input.topLeftCorner(rank, N)
                            .template triangularView<Eigen::Upper>())
        * P.transpose();

    // Form Q.
    dorgqr_(&M, &rank, &rank, input.data(), &M, Tau, W.get(), &LWORK, &INFO);
    L = input.leftCols(rank);

    return false; // Input is not full rank
}

template <typename T>
bool inline Decompose_Matrix(RowMatrix<T> input, RowMatrix<T> &L,
                             RowMatrix<T> &R, double cut) {
    return eigen_version::Decompose_Matrix(input, L, R, cut);
}

void inline ColPivotedQr(
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> input,
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &L,
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &R, double cut) {
    assert(input.size() >= 0);
    int M = input.rows();
    int N = input.cols();
    auto full_rank = std::min(M, N);
    Eigen::VectorXi J = Eigen::VectorXi::Zero(N);
    double Tau[full_rank];
    double work;
    int LWORK = -1; // Ask LAPACK how much space we need.
    int INFO;
    int LDA = M;

    dgeqp3_(&M, &N, input.data(), &LDA, J.data(), Tau, &work, &LWORK, &INFO);
    LWORK = work;
    std::unique_ptr<double[]> W{new double[LWORK]};
    dgeqp3_(&M, &N, input.data(), &LDA, J.data(), Tau, W.get(), &LWORK, &INFO);

    const double thresh = cut;
    auto rank = qr_rank(input, thresh);

    // LAPACK assumes 1 based indexing, but we need zero.
    std::for_each(J.data(), J.data() + J.size(), [](int &val) { --val; });
    Eigen::PermutationWrapper<Eigen::VectorXi> P(J);
    R = Eigen::MatrixXd(input.topLeftCorner(rank, N)
                            .template triangularView<Eigen::Upper>())
        * P.transpose();

    // Form Q.
    dorgqr_(&M, &rank, &rank, input.data(), &M, Tau, W.get(), &LWORK, &INFO);
    L = input.leftCols(rank);
}

template <typename T>
void inline ColPivotedQr(RowMatrix<T> input, RowMatrix<T> &L, RowMatrix<T> &R,
                         double cut) {
    eigen_version::ColPivotedQr(std::move(input), L, R, cut);
}

void inline CompressLeft(
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &L,
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &R, double cut,
    bool debug = false) {
    assert(L.size() >= 0);
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> input = L;
    int M = input.rows();
    int N = input.cols();
    auto full_rank = std::min(M, N);
    Eigen::VectorXi J = Eigen::VectorXi::Zero(N);
    double Tau[full_rank];
    double work;
    int LWORK = -1; // Ask LAPACK how much space we need.
    int INFO;
    int LDA = M;

    dgeqp3_(&M, &N, input.data(), &LDA, J.data(), Tau, &work, &LWORK, &INFO);
    LWORK = work;
    std::unique_ptr<double[]> W{new double[LWORK]};
    dgeqp3_(&M, &N, input.data(), &LDA, J.data(), Tau, W.get(), &LWORK, &INFO);

    const double thresh = cut;
    auto rank = qr_rank(input, thresh);

    if (!debug && rank == full_rank) {
        return;
    }

    // LAPACK assumes 1 based indexing, but we need zero.
    std::for_each(J.data(), J.data() + J.size(), [](int &val) { --val; });
    Eigen::PermutationWrapper<Eigen::VectorXi> P(J);
    R = Eigen::MatrixXd(input.topLeftCorner(rank, N)
                            .template triangularView<Eigen::Upper>())
        * P.transpose() * R;

    // Form Q.
    dorgqr_(&M, &rank, &rank, input.data(), &M, Tau, W.get(), &LWORK, &INFO);
    L = input.leftCols(rank);
}

template <typename T>
void inline CompressLeft(RowMatrix<T> &L, RowMatrix<T> &R, double cut) {
    eigen_version::CompressLeft(L, R, cut);
}

void inline CompressRight(
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &L,
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &R, double cut,
    bool debug = false) {
    assert(R.size() >= 0);
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> input = R;
    int M = input.rows();
    int N = input.cols();
    auto full_rank = std::min(M, N);
    Eigen::VectorXi J = Eigen::VectorXi::Zero(N);
    double Tau[full_rank];
    double work;
    int LWORK = -1; // Ask LAPACK how much space we need.
    int INFO;
    int LDA = M;

    dgeqp3_(&M, &N, input.data(), &LDA, J.data(), Tau, &work, &LWORK, &INFO);
    LWORK = work;
    std::unique_ptr<double[]> W{new double[LWORK]};
    dgeqp3_(&M, &N, input.data(), &LDA, J.data(), Tau, W.get(), &LWORK, &INFO);

    const double thresh = cut;
    auto rank = qr_rank(input, thresh);
    if (!debug && rank == full_rank && rank == 0) {
        return;
    }

    // LAPACK assumes 1 based indexing, but we need zero.
    std::for_each(J.data(), J.data() + J.size(), [](int &val) { --val; });
    Eigen::PermutationWrapper<Eigen::VectorXi> P(J);
    R = Eigen::MatrixXd(input.topLeftCorner(rank, N)
                            .template triangularView<Eigen::Upper>())
        * P.transpose();

    // Form Q.
    dorgqr_(&M, &rank, &rank, input.data(), &M, Tau, W.get(), &LWORK, &INFO);
    L = cblas_gemm(L, input.leftCols(rank), 1.0);
}

template <typename T>
void inline CompressRight(RowMatrix<T> &L, RowMatrix<T> &R, double cut) {
    eigen_version::CompressRight(L, R, cut);
}
} // namespace lapack

} // namespace algebra

#endif // Lin_Algebra_H
