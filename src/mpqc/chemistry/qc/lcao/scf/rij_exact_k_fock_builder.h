/*
 * Created by Drew Lewis on 02/25/2017
 */

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_RIJ_EXACT_K_FOCK_BUILDER_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_RIJ_EXACT_K_FOCK_BUILDER_H_

#include "mpqc/chemistry/qc/lcao/scf/builder.h"
#include "mpqc/math/linalg/cholesky_inverse.h"
#include "mpqc/util/misc/time.h"

#include <cassert>
#include <tiledarray.h>

namespace mpqc {
namespace scf {

template <typename Tile, typename Policy, typename Integral>
class RIJEXACTKBuilder : public FockBuilder<Tile, Policy> {
 public:
  using array_type = typename FockBuilder<Tile, Policy>::array_type;
  array_type L_inv_;
  Integral eri3_;
  Integral eri4_;

  double time_J_;
  double time_K_;
  double energy_J_;
  double energy_K_;

 public:
  RIJEXACTKBuilder(array_type const &M, Integral const &eri3,
                   Integral const &eri4)
      : eri3_(eri3), eri4_(eri4) {
    auto M_eig = array_ops::array_to_eigen(M);

    RowMatrixXd L_inv_eig =
        RowMatrixXd(Eigen::LLT<RowMatrixXd>(M_eig).matrixL()).inverse();

    auto tr_M = M.trange().data()[0];

    L_inv_ = array_ops::eigen_to_array<Tile, Policy>(M.world(), L_inv_eig, tr_M,
                                                     tr_M);
  }

  array_type operator()(array_type const &D, array_type const &C) override {
    auto &world = D.world();
    auto j0 = mpqc::fenced_now(world);
    array_type J;
    J("m, n") =
        eri3_("Z, m, n") *
        (L_inv_("X, Z") * (L_inv_("X, Y") * (eri3_("Y, r, s") * D("r, s"))));
    auto j1 = mpqc::fenced_now(world);
    time_J_ = mpqc::duration_in_s(j0, j1);

    // Make K
    auto k0 = mpqc::fenced_now(world);
    array_type K;
    K("mu, nu") = eri4_("mu, rho, nu, sig") * D("rho, sig");
    auto k1 = mpqc::fenced_now(world);
    time_K_ = mpqc::duration_in_s(k0, k1);

    // Make and return G
    array_type G;
    G("mu, nu") = 2 * J("mu, nu") - K("mu, nu");

    energy_J_ = D("i,j").dot(J("i,j"));
    energy_K_ = D("i,j").dot(K("i,j"));

    return G;
  }

  void register_fock(const array_type &fock,
                     FormulaRegistry<array_type> &registry) override {
    registry.insert(Formula(L"(κ|F|λ)"), fock);
  }

  inline void print_iter(std::string const &leader) override {
    ExEnv::out0() << "Time J: " << time_J_ << ", Energy J: " << energy_J_
                  << "\n";
    ExEnv::out0() << "Time K: " << time_K_ << ", Energy K: " << energy_K_
                  << "\n";
  }
};

}  // namespace scf
}  // namespace mpqc
#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_TRADITIONAL_FOUR_CENTER_FOCK_BUILDER_H_
