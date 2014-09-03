/*
 * File to define linear algebra functions
 */
#ifndef Lin_Algebra_H
#define Lin_Algebra_H

#include "../include/eigen.h"
#include <TiledArray/madness.h>
#include <madness/tensor/clapack.h>

template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> inline cblas_gemm(
    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &A,
    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &B) {

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> C(A.rows(), B.cols());
    const int K = A.cols();
    const int M = C.rows();
    const int N = C.cols();
    const int LDA = M, LDB = K, LDC = M;
    madness::cblas::gemm(madness::cblas::CBLAS_TRANSPOSE::NoTrans,
                         madness::cblas::CBLAS_TRANSPOSE::NoTrans, M, N, K, 1.0,
                         A.data(), LDA, B.data(), LDB, 0.0, C.data(), LDC);
    // assert(C.isApprox(A * B));
    return C;
}

/**
 * ColPivQr computes the column pivoted QR decomposition for a matrix.
 * It returns a boolean which is false if the matrix rank is less than 1/2 the
 * full rank. The function captures input by value and then uses it's space as
 * scratch.  Matrices L and R have Q and R written to them if the input matrix
 * was low rank.
 *
 * Finally the cut parameter is used to set the threshold for determining the
 * rank of the input matrix.
 *
 */
template <typename T>
bool ColPivQR(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> input,
              Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &L,
              Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &R, double cut) {
    int M = input.rows();
    int N = input.cols();
    auto full_rank = std::min(M, N);
    Eigen::VectorXi J = Eigen::VectorXi::Zero(N);
    T Tau[full_rank];
    Eigen::VectorXd Work(1);
    int LWORK = -1; // Ask LAPACK how much space we need.
    int INFO;
    int LDA = M;


    dgeqp3_(&M, &N, input.data(), &LDA, J.data(), Tau, Work.data(), &LWORK,
            &INFO);
    LWORK = Work[0];
    Work.resize(LWORK);
    dgeqp3_(&M, &N, input.data(), &LDA, J.data(), Tau, Work.data(), &LWORK,
            &INFO);

    Eigen::VectorXd const &Rvalues = input.diagonal();
    const double thresh = cut * std::abs(input(1, 1));
    int rank = 0;
    for (auto i = 0; i < Rvalues.size(); ++i) {
        if (std::abs(Rvalues[i]) >= thresh) {
            ++rank;
        } else {
            break;
        }
    }

    if (rank > double(full_rank) / 4.0) {
        return true; // Input is full rank
    }

    // LAPACK assumes 1 based indexing, but we need zero.
    std::for_each(J.data(), J.data() + J.size(), [](int &val) { --val; });
    Eigen::PermutationWrapper<Eigen::VectorXi> P(J);
    R = Eigen::MatrixXd(input.topLeftCorner(rank, N).template triangularView
                        <Eigen::Upper>()) * P.transpose();

    // Form Q.
    dorgqr_(&M, &N, &rank, input.data(), &M, Tau, Work.data(), &LWORK, &INFO);
    L = input.leftCols(rank);

    return false; // Input is not full rank
}

template <typename T>
bool
CompressQR(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> input,
           Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &L,
           Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &R, double cut) {
    int M = input.rows();
    int N = input.cols();
    auto Lfull_rank = std::min(M, N);
    Eigen::VectorXi J = Eigen::VectorXi::Zero(N);
    T Tau[Lfull_rank];
    Eigen::VectorXd Work(1);
    int LWORK = -1; // Ask LAPACK how much space we need.
    int INFO;
    int LDA = M;


    dgeqp3_(&M, &N, input.data(), &LDA, J.data(), Tau, Work.data(), &LWORK,
            &INFO);
    LWORK = Work[0];
    Work.resize(LWORK);
    dgeqp3_(&M, &N, input.data(), &LDA, J.data(), Tau, Work.data(), &LWORK,
            &INFO);

    auto Cfull_rank = std::min(decltype(R.cols())(M), R.cols());
    Eigen::VectorXd const &Rvalues = input.diagonal();
    const double thresh = cut * std::abs(input(1, 1));
    int rank = 0;
    for (auto i = 0; i < Rvalues.size(); ++i) {
        if (std::abs(Rvalues[i]) >= thresh) {
            ++rank;
        } else {
            break;
        }
    }

    // L and R should not be modified until after this point in the algo.
    if (rank > double(Cfull_rank) / 4.0) { // try with 1/4
        return true;                       // Input is full rank
    }

    // LAPACK assumes 1 based indexing, but we need zero.
    std::for_each(J.data(), J.data() + J.size(), [](int &val) { --val; });
    Eigen::PermutationWrapper<Eigen::VectorXi> P(J);
    R = Eigen::MatrixXd(input.topLeftCorner(rank, N).template triangularView
                    <Eigen::Upper>()) * P.transpose() * R;

    // Form Q.
    dorgqr_(&M, &N, &rank, input.data(), &M, Tau, Work.data(), &LWORK, &INFO);
    L = input.leftCols(rank);

    return false; // Input is not full rank
}

#endif // Lin_Algebra_H
