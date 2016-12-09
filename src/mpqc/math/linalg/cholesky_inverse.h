#ifndef MPQC4_SRC_MPQC_MATH_LINALG_CHOLESKY_INVERSE_H_
#define MPQC4_SRC_MPQC_MATH_LINALG_CHOLESKY_INVERSE_H_

#include <tiledarray.h>
#include "mpqc/math/tensor/clr/array_to_eigen.h"

namespace mpqc {
namespace array_ops {

template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> cholesky_inverse(
    TA::DistArray<Tile, Policy> const &A) {
  auto &world = A.world();

  auto A_eig = array_to_eigen(A);
  Eigen::LLT<decltype(A_eig)> llt(A_eig);

  decltype(A_eig) L_inv_eig = decltype(A_eig)(llt.matrixL()).inverse();

  auto tr_A0 = A.trange().data()[0];
  auto tr_A1 = A.trange().data()[1];

  return array_ops::eigen_to_array<Tile, Policy>(world, L_inv_eig, tr_A0, tr_A1);
}

}  // namespace array_ops
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_MATH_LINALG_CHOLESKY_INVERSE_H_
