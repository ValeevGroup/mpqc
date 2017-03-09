//
// Created by Chong Peng on 2/22/17.
//

#include <tiledarray.h>

#include "catch.hpp"
#include "mpqc/math/external/eigen/eigen.h"
#include "mpqc/math/linalg/gram_schmidt.h"
#include "mpqc/math/tensor/clr/array_to_eigen.h"

using namespace mpqc;

TEST_CASE("Gram Schmidt", "[gram-schmidt]") {
  using Array = TA::DistArray<TA::TensorD, TA::DensePolicy>;

  const auto n = 400;  // vector size
  const auto v = 40;   // number of vector

  TA::TiledRange1 tr_n{0, 100, 200, 300, n};
  TA::TiledRange1 tr_v{0, 1};

  std::vector<Array> vecs(v);

  // initialize vector
  for (auto i = 0; i < v; i++) {
    EigenVector<double> vec = double(i + 1) * EigenVector<double>::Random(n);
    vecs[i] = array_ops::eigen_to_array<TA::TensorD, TA::DensePolicy>(
        TA::get_default_world(), vec, tr_n, tr_v);
  }

  gram_schmidt(vecs);

  const double tolerance = std::numeric_limits<double>::epsilon() * 100;

  for (auto i = 0; i < v; ++i) {
    for (auto j = i; j < v; ++j) {
      const auto test = dot_product(vecs[i], vecs[j]);
      //      std::cout << "i= " << i << " j= " << j << " dot= " << test <<
      //      std::endl;
      if (i == j) {
        // test if normalized
        REQUIRE(test == Approx(1.0).epsilon(tolerance));
      } else {
        // test if orthogonalized
        REQUIRE(test == Approx(0.0).epsilon(tolerance));
      }
    }
  }
}