//
// Created by Chong Peng on 11/21/16.
//

#ifndef MPQC4_SRC_MPQC_MATH_LINALG_CONDITIONED_ORTHOGONALIZER_H_
#define MPQC4_SRC_MPQC_MATH_LINALG_CONDITIONED_ORTHOGONALIZER_H_

#include <tiledarray.h>
#include "mpqc/math/external/eigen/eigen.h"
#include "mpqc/math/tensor/clr/array_to_eigen.h"
#include "mpqc/math/external/eigen/eigen.h"

namespace mpqc {
namespace array_ops {

template<typename T>
std::tuple<RowMatrix<T>, RowMatrix<T>, size_t, double, double> gensqrtinv(
    const RowMatrix<T>& S, bool symmetric = false,
    double max_condition_number = 1e8) {

  using Matrix = RowMatrix<T>;

  Eigen::SelfAdjointEigenSolver<Matrix> eig_solver(S);
  auto U = eig_solver.eigenvectors();
  auto s = eig_solver.eigenvalues();
  auto s_max = s.maxCoeff();
  auto condition_number = std::min(
      s_max / std::max(s.minCoeff(), std::numeric_limits<double>::min()),
      1.0 / std::numeric_limits<double>::epsilon());
  auto threshold = s_max / max_condition_number;
  long n = s.rows();
  long n_cond = 0;
  for (long i = n - 1; i >= 0; --i) {
    if (s(i) >= threshold) {
      ++n_cond;
    } else
      i = 0;  // skip rest since eigenvalues are in ascending order
  }

  auto sigma = s.bottomRows(n_cond);
  auto result_condition_number = sigma.maxCoeff() / sigma.minCoeff();
  auto sigma_sqrt = sigma.array().sqrt().matrix().asDiagonal();
  auto sigma_invsqrt = sigma.array().sqrt().inverse().matrix().asDiagonal();

  // make canonical X/Xinv
  auto U_cond = U.block(0, n - n_cond, n, n_cond);
  Matrix X = U_cond * sigma_invsqrt;
  Matrix Xinv = U_cond * sigma_sqrt;
  // convert to symmetric, if needed
  if (symmetric) {
    X = X * U_cond.transpose();
    Xinv = Xinv * U_cond.transpose();
  }
  return std::make_tuple(X, Xinv, size_t(n_cond), condition_number,
                         result_condition_number);
}

template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> conditioned_orthogonalizer(
    TA::DistArray <Tile, Policy> S_array, double S_condition_number_threshold) {

  using Matrix = RowMatrix<typename Tile::numeric_type>;
  auto& world = S_array.world();

  Matrix S = array_ops::array_to_eigen(S_array);

  size_t obs_rank;
  double S_condition_number;
  double XtX_condition_number;
  Matrix X, Xinv;

  assert(S.rows() == S.cols());

  std::tie(X, Xinv, obs_rank, S_condition_number, XtX_condition_number) =
      gensqrtinv(S, true, S_condition_number_threshold);
  auto obs_nbf_omitted = (long)S.rows() - (long)obs_rank;

  if(world.rank() == 0){
    std::cout << "overlap condition number = " << S_condition_number;
  }

  if (obs_nbf_omitted > 0 && world.rank() == 0){
    std::cout << " (dropped " << obs_nbf_omitted << " "
              << (obs_nbf_omitted > 1 ? "fns" : "fn") << " to reduce to "
              << XtX_condition_number << ")";
    std::cout << std::endl;
  }

  if (obs_nbf_omitted > 0 && world.rank()==0)  {
    Matrix should_be_I = X.transpose() * S * X;
    Matrix I = Matrix::Identity(should_be_I.rows(), should_be_I.cols());
    std::cout << "||X^t * S * X - I||_2 = " << (should_be_I - I).norm()
              << " (should be 0)" << std::endl;
  }

//  return std::make_tuple(X, Xinv, XtX_condition_number);
  return array_ops::eigen_to_array<Tile,Policy>(world,X,S_array.trange().data()[0], S_array.trange().data()[1]);
}

}  // namespace  array_ops
}  // namespace  mpqc

#endif  // MPQC4_SRC_MPQC_MATH_LINALG_CONDITIONED_ORTHOGONALIZER_H_
