//
// Created by Chong Peng on 10/6/17.
//

#include "catch.hpp"

#include "mpqc/math/external/eigen/eigen.h"
#include "mpqc/math/external/tiledarray/array_max_n.h"
#include "mpqc/math/tensor/clr/array_to_eigen.h"

using namespace mpqc;

TEST_CASE("Max N Elements of TA::DistArray", "[array-max-n]") {
  // matrix size
  std::size_t n = 500;
  std::size_t k = 5;
  const auto sparse = 0.1;

  // initialize matrix
  RowMatrix<double> A = RowMatrix<double>::Zero(n, n);
  for (auto i = 0; i < n; i++) {
    A(i, i) = i + 1;
  }
  A = A + sparse * RowMatrix<double>::Random(n, n);

  std::vector<double> ref_result = {A(n - 1, n - 1), A(n - 2, n - 2),
                                    A(n - 3, n - 3), A(n - 4, n - 4),
                                    A(n - 5, n - 5)};

  TA::TiledRange1 tr_n{0, 100, 200, 300, n};

  auto A_ta = math::eigen_to_array<TA::TensorD, TA::DensePolicy>(
      TA::get_default_world(), A, tr_n, tr_n);

  /// check array_max_n
  std::vector<double> result = array_max_n(A_ta, k);

  CHECK(result == ref_result);


  /// check array_max_n_index
  auto index_result = array_max_n_index(A_ta, k);
  for(std::size_t i = 0; i < k; i++){

    CHECK(index_result[i].first == result[i]);
    CHECK(index_result[i].second[0] == n-i-1);
    CHECK(index_result[i].second[1] == n-i-1);

  }
}

TEST_CASE("Abs Max N Elements of TA::DistArray", "[array-abs-max-n]") {
  // matrix size
  std::size_t n = 500;
  std::size_t k = 5;
  const auto sparse = 0.1;

  // initialize matrix
  RowMatrix<double> A = RowMatrix<double>::Zero(n, n);
  for (auto i = 0; i < n; i++) {
    A(i, i) =  (i%2==0) ? (i + 1) : -(i + 1);
  }
  A = A + sparse * RowMatrix<double>::Random(n, n);

  std::vector<double> ref_result = {A(n - 1, n - 1), A(n - 2, n - 2),
                                    A(n - 3, n - 3), A(n - 4, n - 4),
                                    A(n - 5, n - 5)};

  TA::TiledRange1 tr_n{0, 100, 200, 300, n};

  auto A_ta = math::eigen_to_array<TA::TensorD, TA::DensePolicy>(
      TA::get_default_world(), A, tr_n, tr_n);

  std::vector<double> result = array_abs_max_n(A_ta, k);

  CHECK(result == ref_result);

  /// check array_max_n_index
  auto index_result = array_abs_max_n_index(A_ta, k);
  for(std::size_t i = 0; i < k; i++){

    CHECK(index_result[i].first == result[i]);
    CHECK(index_result[i].second[0] == n-i-1);
    CHECK(index_result[i].second[1] == n-i-1);

  }
}