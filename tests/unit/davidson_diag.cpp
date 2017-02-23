//
// Created by Chong Peng on 2/21/17.
//

#include "mpqc/math/linalg/davidson_diag.h"
#include "mpqc/math/tensor/clr/array_to_eigen.h"

#include "catch.hpp"

using namespace mpqc;

TEST_CASE("Davidson Algorithm", "[davidson]") {
  // matrix size
  const auto n = 200;
  const auto sparse = 0.01;
  const auto n_roots = 3;
  const auto n_guess = 3;
  const auto converge = 1.0e-6;
  const auto max_iter = 20;

  // initialize matrix
  RowMatrix<double> A = RowMatrix<double>::Zero(n, n);
  for (auto i = 0; i < n; i++) {
    A(i, i) = i + 1;
  }
  A = A + sparse * RowMatrix<double>::Random(n, n);
  A = 0.5 * (A.transpose() + A);

  // eigen solve
  Eigen::SelfAdjointEigenSolver<RowMatrix<double>> es(A);
  auto e = es.eigenvalues().segment(0, n_roots);

  std::cout << e << std::endl;

  TA::TiledRange1 tr_n{0, 50, 100, 150, n};
  TA::TiledRange1 tr_guess{0, 1};

  auto A_ta = array_ops::eigen_to_array<TA::TensorD, TA::DensePolicy>(
      TA::get_default_world(), A, tr_n, tr_n);

  // build guess vector
  RowMatrix<double> guess = RowMatrix<double>::Identity(n, n_guess);

  std::vector<TA::DistArray<TA::TensorD, TA::DensePolicy>> guess_ta(n_guess);

  for (auto i = 0; i < n_guess; i++) {
    guess_ta[i] = array_ops::eigen_to_array<TA::TensorD, TA::DensePolicy>(
        TA::get_default_world(), guess.col(i), tr_n, tr_guess);
  }

  SymmDavidsonDiag<TA::DistArray<TA::TensorD, TA::DensePolicy>> dvd(n_roots,
                                                                    n_guess);

  EigenVector<double> eig = EigenVector<double>::Zero(n_roots);
  for (auto i = 0; i < max_iter; i++) {
    std::cout << "Iter: " << i << std::endl;
    const auto n_v = guess_ta.size();
    std::vector<TA::DistArray<TA::TensorD, TA::DensePolicy>> HB(n_v);

    for(auto i = 0; i < n_v; i++){
      HB[i]("i,j") = A_ta("i,k") * guess_ta[i]("k,j");
    }

    EigenVector<double> eig_new = dvd.extrapolate(HB, guess_ta);

    std::cout << eig_new << std::endl;
    std::cout << "norm= " << (eig - eig_new).norm() << "\n";

    if ((eig - eig_new).norm() < converge) {
      break;
    }

    eig = eig_new;
  }

  CHECK((e - eig).norm() < converge);
}