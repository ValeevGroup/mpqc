
#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_TRADITIONAL_FOUR_CENTER_FOCK_BUILDER_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_TRADITIONAL_FOUR_CENTER_FOCK_BUILDER_H_

#include <cassert>

#include <tiledarray.h>

#include "mpqc/chemistry/qc/scf/builder.h"

namespace mpqc {
namespace scf {

template <typename Tile, typename Policy, typename Integral>
class FourCenterBuilder : public FockBuilder<Tile,Policy> {
 public:
  using array_type = typename FockBuilder<Tile,Policy>::array_type;
  Integral eri4_;

 public:
  FourCenterBuilder(Integral const &eri4) : eri4_(eri4) {}

  array_type operator()(array_type const &D, array_type const &C) override {
    // Make J
    array_type J;
    J("mu, nu") = eri4_("mu, nu, rho, sig") * D("rho, sig");

    // Make K
    array_type K;
    K("mu, nu") = eri4_("mu, rho, nu, sig") * D("rho, sig");

    // Make and return G
    array_type G;
    G("mu, nu") = 2 * J("mu, nu") - K("mu, nu");

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
