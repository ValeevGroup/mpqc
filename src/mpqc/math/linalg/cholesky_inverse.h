#ifndef MPQC4_SRC_MPQC_MATH_LINALG_CHOLESKY_INVERSE_H_
#define MPQC4_SRC_MPQC_MATH_LINALG_CHOLESKY_INVERSE_H_

#include "mpqc/math/tensor/clr/array_to_eigen.h"
#include "mpqc/util/core/exenv.h"

#include <tiledarray.h>

namespace mpqc {
namespace array_ops {

template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> cholesky_inverse(
    TA::DistArray<Tile, Policy> const &A) {
  auto &world = A.world();
  if (world.size() > 1) {
    ExEnv::out0() << "Warning: Using cholesky inverse in a world with more "
                     "than 1 process may result in poor scaling! World size: "
                  << world.size() << "\n";
  }

  auto A_eig = ::mpqc::array_ops::array_to_eigen(A);
  Eigen::LLT<decltype(A_eig)> llt(A_eig);

  Eigen::ComputationInfo info = llt.info();
  if (info != Eigen::ComputationInfo::Success){
    throw AlgorithmException("Eigen Cholesky Inverse Failed!",__FILE__,__LINE__);
  }

  decltype(A_eig) L_inv_eig = decltype(A_eig)(llt.matrixL()).inverse();

  auto tr_A0 = A.trange().data()[0];
  auto tr_A1 = A.trange().data()[1];

  return array_ops::eigen_to_array<Tile, Policy>(world, L_inv_eig, tr_A0,
                                                 tr_A1);
}

template <typename Tile, typename Policy>
TA::DistArray<Tile,Policy> eigen_inverse(const TA::DistArray<Tile, Policy> &A){

  auto& world = A.world();
  auto result_eig = array_ops::array_to_eigen(A);

  using Matrix = decltype(result_eig);
  // compute cholesky decomposition
  auto llt_solver = Eigen::LLT<Matrix>(result_eig);

  // check success
  Eigen::ComputationInfo info = llt_solver.info();
  if (info == Eigen::ComputationInfo::Success) {
    Matrix L = Matrix(llt_solver.matrixL());
    Matrix L_inv_eig = L.inverse();
    result_eig = L_inv_eig.transpose() * L_inv_eig;
  } else if (info == Eigen::ComputationInfo::NumericalIssue) {
    utility::print_par(world,
                       "!!!\nWarning!! NumericalIssue in Cholesky "
                           "Decomposition\n!!!\n");
  } else if (info == Eigen::ComputationInfo::NoConvergence) {
    utility::print_par(world,
                       "!!!\nWarning!! NoConvergence in Cholesky "
                           "Decomposition\n!!!\n");
  }

  if (info != Eigen::ComputationInfo::Success) {
    utility::print_par(world, "Using Eigen LU Decomposition Inverse!\n");

    Eigen::FullPivLU<Matrix> lu(result_eig);

    TA_ASSERT(lu.isInvertible());

    result_eig = lu.inverse();
  }

  auto tr_1 = A.trange().data()[0];
  auto tr_2 = A.trange().data()[1];
  auto result = array_ops::eigen_to_array<Tile, Policy>( A.world(), result_eig,
                                                   tr_1, tr_2);
  return result;
  
};

}  // namespace array_ops
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_MATH_LINALG_CHOLESKY_INVERSE_H_
