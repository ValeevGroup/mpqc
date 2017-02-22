//
// Created by Chong Peng on 2/21/17.
//

#include "mpqc/math/linalg/davidson_diag.h"

#include <tiledarray.h>
#include "catch.hpp"

using namespace mpqc;

TEST_CASE("Davidson Algorithm", "[davidson]"){

  // matrix size
  const auto n = 100;
  const auto sparse = 0.01;
  const auto n_roots = 3;
  const auto n_guess = 3;
  const auto converge = 1.0e-5;
  const auto max_iter = 100;

  // initialize matrix
  RowMatrix<double> A = RowMatrix<double>::Zero(n,n);
  for(auto i = 0; i < n; i++){
    A(i,i) = i+1;
  }
  A = A + sparse*RowMatrix<double>::Random(n,n);
  A = 0.5*(A.transpose() + A);

  // eigen solve
  Eigen::SelfAdjointEigenSolver<RowMatrix<double>> es(A);
  auto e = es.eigenvalues().segment(0,n_roots);

  std::cout << e << std::endl;

  TA::TiledRange1 tr_n {0,50,100};
  TA::TiledRange1 tr_guess {0,n_guess};

  auto A_ta = array_ops::eigen_to_array<TA::TensorD, TA::DensePolicy>(TA::get_default_world(), A, tr_n, tr_n);

  // build guess vector
  RowMatrix<double> guess = RowMatrix<double>::Identity(n,n_guess);
  auto guess_ta = array_ops::eigen_to_array<TA::TensorD, TA::DensePolicy>(TA::get_default_world(), guess, tr_n, tr_guess);

  SymmDavidsonDiag<TA::DistArray<TA::TensorD,TA::DensePolicy> > dvd (n_roots,n_guess);

  EigenVector<double> eig = EigenVector<double>::Zero(n_roots);
  for(auto i = 0; i < 10; i++){
    std::cout << "Iter: " << i << std::endl;
    TA::DistArray<TA::TensorD, TA::DensePolicy> HB;
    HB("i,j") = A_ta("i,k")* guess_ta("k,j");

    EigenVector<double> eig_new = dvd.extrapolate(HB, guess_ta);

    std::cout << eig_new << std::endl;

    if ( (eig - eig_new).norm() < converge ){
      break;
    }

    eig = eig_new;
  }

  CHECK((e - eig).norm() < converge);

}