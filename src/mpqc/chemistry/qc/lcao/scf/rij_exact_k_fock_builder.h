/*
 * Created by Drew Lewis on 02/25/2017
 */

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_RIJ_EXACT_K_FOCK_BUILDER_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_RIJ_EXACT_K_FOCK_BUILDER_H_

#include "mpqc/chemistry/qc/lcao/scf/builder.h"
#include "mpqc/chemistry/qc/lcao/scf/traditional_four_center_fock_builder.h"
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
  using Basis = ::mpqc::lcao::gaussian::Basis;

  RIJEXACTKBuilder(array_type const &M, Integral const &eri3,
                   std::shared_ptr<const Basis> bra_basis,
                   std::shared_ptr<const Basis> ket_basis,
                   std::shared_ptr<const Basis> density_basis)
      : eri3_(eri3),
        K_builder_(std::make_shared<Direct4CBuilder>(
            M.world(), bra_basis, ket_basis, density_basis, false, true)) {
    auto M_eig = array_ops::array_to_eigen(M);

    RowMatrixXd L_inv_eig =
        RowMatrixXd(Eigen::LLT<RowMatrixXd>(M_eig).matrixL()).inverse();

    auto tr_M = M.trange().data()[0];

    L_inv_ = array_ops::eigen_to_array<Tile, Policy>(M.world(), L_inv_eig, tr_M,
                                                     tr_M);
  }

  array_type operator()(array_type const &D, array_type const &C,
                        double target_precision =
                            std::numeric_limits<double>::epsilon()) override {
    auto &world = D.world();
    auto j0 = mpqc::fenced_now(world);
    array_type J;
    J("m, n") =
        eri3_("Z, m, n") *
        (L_inv_("X, Z") * (L_inv_("X, Y") * (eri3_("Y, r, s") * D("r, s"))));
    // symmetrize to account for permutational symmetry
    J("m, n") = 0.5 * (J("m, n") + J("n, m"));
    auto j1 = mpqc::fenced_now(world);
    time_J_ = mpqc::duration_in_s(j0, j1);

    // Make K
    auto k0 = mpqc::fenced_now(world);
    array_type minus_K = (*K_builder_)(D, C, target_precision);
    auto k1 = mpqc::fenced_now(world);
    time_K_ = mpqc::duration_in_s(k0, k1);

    // Make and return G
    array_type G;
    G("mu, nu") = 2 * J("mu, nu") + minus_K("mu, nu");

    energy_J_ = D("i,j").dot(J("i,j"));
    energy_K_ = -D("i,j").dot(minus_K("i,j"));

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

 private:
  array_type L_inv_;
  Integral eri3_;
  using Direct4CBuilder = FourCenterFockBuilder<Tile, Policy>;
  std::shared_ptr<Direct4CBuilder> K_builder_;

  double time_J_;
  double time_K_;
  double energy_J_;
  double energy_K_;
};

}  // namespace scf
}  // namespace mpqc
#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_TRADITIONAL_FOUR_CENTER_FOCK_BUILDER_H_
