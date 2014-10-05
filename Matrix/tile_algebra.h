/*
 * File to define linear algebra functions
 */
#ifndef Lin_Algebra_H
#define Lin_Algebra_H

#include "../include/eigen.h"
#include <TiledArray/madness.h>
#include <chrono>
#include <madness/tensor/clapack.h>

namespace algebra {

namespace eigen_version {
template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> inline cblas_gemm(
    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &A,
    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &B, double alpha) {
    return alpha * A * B;
}

template <typename T>
void inline cblas_gemm_inplace(
    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &A,
    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &B,
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &C, double alpha,
    double beta = 1.0) {
    C = alpha * A * B + beta * C;
}

template <typename T>
bool inline Decompose_Matrix(
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> input,
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &L,
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &R, double cut) {
    Eigen::ColPivHouseholderQR<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>
        qr(input);
    qr.setThreshold(cut);
    auto rank = qr.rank();

    auto full_rank = std::min(input.cols(), input.rows());
    if (rank >= double(full_rank) / 2.0) {
        return true;
    }

    R = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(
            qr.matrixR()
                .topLeftCorner(rank, input.cols())
                .template triangularView<Eigen::Upper>())
        * qr.colsPermutation().transpose();
    L = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(qr.householderQ())
            .leftCols(rank);

    return false;
}

template <typename T>
void inline ColPivotedQr(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> input,
                         Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &L,
                         Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &R,
                         double cut) {
    Eigen::ColPivHouseholderQR<typename std::remove_reference<decltype(
        input)>::type> qr(input);

    qr.setThreshold(cut);
    auto rank = qr.rank();

    R = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(
            qr.matrixR()
                .topLeftCorner(rank, input.cols())
                .template triangularView<Eigen::Upper>())
        * qr.colsPermutation().transpose();

    L = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(qr.householderQ())
            .leftCols(rank);
}

template <typename T>
void inline CompressLeft(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &L,
                         Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &R,
                         double cut) {

    Eigen::ColPivHouseholderQR<typename std::remove_reference<decltype(
        L)>::type> qr(L);

    qr.setThreshold(cut);
    auto rank = qr.rank();

    if (rank >= L.cols()) {
        return;
    }

    R = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(
            qr.matrixR()
                .topLeftCorner(rank, qr.matrixQR().cols())
                .template triangularView<Eigen::Upper>())
        * qr.colsPermutation().transpose() * R;

    L = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(qr.householderQ())
            .leftCols(rank);
}

template <typename T>
void inline CompressRight(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &L,
                          Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &R,
                          double cut) {
    Eigen::ColPivHouseholderQR<typename std::remove_reference<decltype(
        R)>::type> qr(R);

    qr.setThreshold(cut);
    auto rank = qr.rank();

    if (rank >= R.rows()) {
        return;
    }

    R = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(
            qr.matrixR()
                .topLeftCorner(rank, qr.matrixQR().cols())
                .template triangularView<Eigen::Upper>())
        * qr.colsPermutation().transpose();

    L = L
        * Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(qr.householderQ())
              .leftCols(rank);
}

} // namespace eigen_version

inline namespace lapack {
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> inline cblas_gemm(
    const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &A,
    const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &B,
    double alpha) {

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> C(A.rows(), B.cols());
    const int K = A.cols();
    const int M = A.rows();
    const int N = B.cols();
    const int LDA = M, LDB = K, LDC = M;
    madness::cblas::gemm(madness::cblas::CBLAS_TRANSPOSE::NoTrans,
                         madness::cblas::CBLAS_TRANSPOSE::NoTrans, M, N, K,
                         alpha, A.data(), LDA, B.data(), LDB, 0.0, C.data(),
                         LDC);
    return C;
}

template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> inline cblas_gemm(
    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &A,
    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &B, double alpha) {
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
    const int M = A.rows();
    const int N = B.cols();
    const int LDA = M, LDB = K, LDC = M;
    madness::cblas::gemm(madness::cblas::CBLAS_TRANSPOSE::NoTrans,
                         madness::cblas::CBLAS_TRANSPOSE::NoTrans, M, N, K,
                         alpha, A.data(), LDA, B.data(), LDB, beta, C.data(),
                         LDC);
}

template <typename T>
void inline cblas_gemm_inplace(
    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &A,
    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &B,
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &C, double alpha,
    double beta = 1.0) {
    eigen_version::cblas_gemm_inplace(A, B, C, alpha, beta);
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

    Eigen::VectorXd const &Rvalues = input.diagonal();
    const double thresh = std::max(cut * std::abs(Rvalues[0]), 1e-16);
    int rank = 0;
    for (auto i = 0; i < Rvalues.size(); ++i) {
        if (std::abs(Rvalues[i]) > thresh) {
            ++rank;
        } else {
            break;
        }
    }

    if (rank > double(full_rank) / 2.0) {
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
bool inline Decompose_Matrix(
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> input,
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &L,
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &R, double cut) {
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

    Eigen::VectorXd const &Rvalues = input.diagonal();
    if (Rvalues.size() == 0) {
        return;
    }
    const double thresh = std::max(cut * std::abs(Rvalues[0]), 1e-16);
    int rank = 0;
    for (auto i = 0; i < Rvalues.size(); ++i) {
        if (std::abs(Rvalues[i]) > thresh) {
            ++rank;
        } else {
            break;
        }
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
}

template <typename T>
void inline ColPivotedQr(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> input,
                         Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &L,
                         Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &R,
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

    Eigen::VectorXd const &Rvalues = input.diagonal();
    const double thresh = std::max(cut * std::abs(Rvalues[0]), 1e-16);
    int rank = 0;
    for (auto i = 0; i < Rvalues.size(); ++i) {
        if (std::abs(Rvalues[i]) > thresh) {
            ++rank;
        } else {
            break;
        }
    }

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
void inline CompressLeft(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &L,
                         Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &R,
                         double cut, bool debug = false) {
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

    Eigen::VectorXd const &Rvalues = input.diagonal();
    const double thresh = std::max(cut * std::abs(Rvalues[0]), 1e-16);
    int rank = 0;
    for (auto i = 0; i < Rvalues.size(); ++i) {
        if (std::abs(Rvalues[i]) > thresh) {
            ++rank;
        } else {
            break;
        }
    }

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
void inline CompressRight(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &L,
                          Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &R,
                          double cut, bool debug = false) {

    eigen_version::CompressRight(L, R, cut);
}
} // namespace lapack

} // namespace algebra

#endif // Lin_Algebra_H
