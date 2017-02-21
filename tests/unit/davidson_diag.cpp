//
// Created by Chong Peng on 2/21/17.
//

#include "mpqc/math/linalg/davidson_diag.h"

#include <tiledarray.h>
#include "catch.hpp"

using namespace mpqc;

TEST_CASE("Davidson Algorithm", "[davidson]"){

  SymmDavidsonDiag<TA::DistArray<TA::TensorD,TA::DensePolicy> > (5,5);

}