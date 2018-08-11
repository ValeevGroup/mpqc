//
// Created by Eduard Valeyev on 7/12/18.
//

#include "catch.hpp"
#include "mpqc/math/external/eigen/random.h"
#include <tiledarray.h>

using namespace mpqc;

TEST_CASE("random_unitary", "[eigen]") {
  // only do in serial since this only uses Eigen
  if (TA::get_default_world().size() > 1)
    return;

  const auto size = 20;
  auto U = mpqc::math::random_unitary<double>(size);
  REQUIRE((U.transpose() * U - RowMatrix<double>(size,size).setIdentity()).lpNorm<Eigen::Infinity>() < 1e-15);
  REQUIRE((U * U.transpose() - RowMatrix<double>(size,size).setIdentity()).lpNorm<Eigen::Infinity>() < 1e-15);
}