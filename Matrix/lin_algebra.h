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
    madness::cblas::gemm(madness::cblas::CBLAS_TRANSPOSE::NoTrans,
                         madness::cblas::CBLAS_TRANSPOSE::NoTrans, A.rows(),
                         B.cols(), B.rows(), 1.0, A.data(), 1, B.data(), 1, 1.0,
                         C.data(), 1);
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
    int rank = std::count_if(Rvalues.data(), Rvalues.data() + Rvalues.size(),
                             [=](T x) { return std::abs(x) > thresh; });

    if (rank > double(full_rank) / 2.0) {
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


#endif // Lin_Algebra_H
