
#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_TRADITIONAL_FOUR_CENTER_FOCK_BUILDER_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_TRADITIONAL_FOUR_CENTER_FOCK_BUILDER_H_

#include <cassert>

#include <tiledarray.h>

#include "mpqc/chemistry/qc/lcao/scf/builder.h"

namespace mpqc {
namespace scf {

template <typename Tile, typename Policy, typename Integral>
class ReferenceFourCenterFockBuilder : public FockBuilder<Tile, Policy> {
 public:
  using array_type = typename FockBuilder<Tile, Policy>::array_type;
  Integral eri4_J_;
  Integral eri4_K_;

 public:
  ReferenceFourCenterFockBuilder(Integral const &eri4_J, Integral const &eri4_K)
      : eri4_J_(eri4_J), eri4_K_(eri4_K) {}

  array_type operator()(array_type const &D, array_type const &C) override {
    const auto make_J = eri4_J_.is_initialized();
    const auto make_K = eri4_K_.is_initialized();
    assert(make_J || make_K);

    // Make J
    array_type J;
    if (make_J) {
      J("mu, nu") = eri4_J_("mu, nu, rho, sig") * D("rho, sig");
      // symmetrize to account for petite list
      J("mu, nu") = 0.5 * (J("mu, nu") + J("nu, mu"));
    }

    // Make K
    array_type K;
    if (make_K) {
      K("mu, nu") = eri4_K_("mu, rho, nu, sig") * D("rho, sig");
      // symmetrize to account for petite list
      K("mu, nu") = 0.5 * (K("mu, nu") + K("nu, mu"));
    }

    // Make and return G
    array_type G;
    if (make_J && make_K)
      G("mu, nu") = 2 * J("mu, nu") - K("mu, nu");
    else if (make_J)
      G("mu, nu") = 2 * J("mu, nu");
    else if (make_K)
      G("mu, nu") = -K("mu, nu");

    return G;
  }

  void register_fock(const array_type &fock,
                     FormulaRegistry<array_type> &registry) override {
    registry.insert(Formula(L"(κ|F|λ)"), fock);
  }

  inline void print_iter(std::string const &leader) override {}
};

}  // namespace scf
}  // namespace mpqc
#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_TRADITIONAL_FOUR_CENTER_FOCK_BUILDER_H_
