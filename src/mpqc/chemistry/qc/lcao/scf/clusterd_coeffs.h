
#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_CLUSTERD_COEFFS_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_CLUSTERD_COEFFS_H_

#include "mpqc/math/external/eigen/eigen.h"
#include <tiledarray.h>

#include "mpqc/math/external/eigen/eigen.h"

#include "mpqc/math/tensor/clr/vector_localization.h"

namespace mpqc {
namespace lcao {
namespace scf {

template <typename Tile, typename Policy>
void clustered_coeffs(
    std::vector<TA::DistArray<Tile, Policy>> const &xyz,
    TA::DistArray<Tile, Policy> &C, unsigned long occ_nclusters) {
  TA::DistArray<Tile, Policy> X, Y, Z;
  X("i,j") = C("mu,i") * xyz[0]("mu, nu") * C("nu,j");
  Y("i,j") = C("mu,i") * xyz[1]("mu, nu") * C("nu,j");
  Z("i,j") = C("mu,i") * xyz[2]("mu, nu") * C("nu,j");

  auto X_eig = array_ops::array_to_eigen(X);
  auto Y_eig = array_ops::array_to_eigen(Y);
  auto Z_eig = array_ops::array_to_eigen(Z);

  decltype(X_eig) oc_pos(X_eig.rows(), 3);
  for (auto i = 0; i < X_eig.rows(); ++i) {
    oc_pos(i, 0) = X_eig(i, i);
    oc_pos(i, 1) = Y_eig(i, i);
    oc_pos(i, 2) = Z_eig(i, i);
  }

  auto C_eig = array_ops::array_to_eigen(C);
  auto tr_occ =
      tensor::localize_vectors_with_kmeans(oc_pos, C_eig, occ_nclusters);

  C = array_ops::eigen_to_array<Tile,Policy>(C.world(), C_eig, C.trange().data()[0],
                                      tr_occ);
}

}  // namespace scf
}  // namespace lcao
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_CLUSTERD_COEFFS_H_
