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
  int LWORK = -1;
  int INFO;
  int LDA = M;
  dgeqp3_(&M, &N, input.data(), &LDA, J.data(), Tau, Work.data(), &LWORK,
          &INFO);
  LWORK = Work[0];
  Work.resize(LWORK);
  dgeqp3_(&M, &N, input.data(), &LDA, J.data(), Tau, Work.data(), &LWORK,
          &INFO);
  int rank = 0;
  double A00 = std::abs(input(1, 1));
  for (auto k = 0; k < full_rank; ++k) {
    if (std::abs(input(k, k)) < cut * A00) {
      rank = k;
      break;
    }
  }

  //if (rank > double(full_rank) / 2.0) {
  //  return false;
  //}

  for (auto i = 0; i < N; ++i) {
    --J[i];
  }
  Eigen::PermutationWrapper<Eigen::VectorXi> P(J);

  R = Eigen::MatrixXd(input.topLeftCorner(rank, N).template triangularView
                      <Eigen::Upper>()) * P.transpose();


  dorgqr_(&M, &N, &rank, input.data(), &M, Tau, Work.data(), &LWORK,
          &INFO);

  L = input.leftCols(rank);

  return true;
}


#endif // Lin_Algebra_H
