#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_DF_FOCK_BUILDER_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_DF_FOCK_BUILDER_H_

#include "mpqc/chemistry/qc/lcao/scf/builder.h"
#include "mpqc/chemistry/qc/lcao/scf/pbc/periodic_four_center_fock_builder.h"
#include "mpqc/chemistry/qc/lcao/scf/pbc/periodic_ri_j_builder.h"
#include "mpqc/chemistry/qc/lcao/scf/pbc/util.h"

namespace mpqc {
namespace scf {

template <typename Tile, typename Policy, typename Factory>
class PeriodicDFFockBuilder : public PeriodicFockBuilder<Tile, Policy> {
 public:
  using array_type = typename PeriodicFockBuilder<Tile, Policy>::array_type;
  using J_Builder = PeriodicRIJBuilder<Tile, Policy, Factory>;
  using K_Builder = PeriodicFourCenterFockBuilder<Tile, Policy>;

  PeriodicDFFockBuilder(Factory &ao_factory) : ao_factory_(ao_factory) {
    auto &world = ao_factory_.world();

    // Construct periodic RI-J builder
    auto t0_j_init = mpqc::fenced_now(world);
    j_builder_ = std::make_unique<J_Builder>(ao_factory_);
    auto t1_j_init = mpqc::fenced_now(world);
    double t_j_init = mpqc::duration_in_s(t0_j_init, t1_j_init);

    // Construct exact perioic 4-center K builder
    auto t0_k_init = mpqc::fenced_now(world);
    {
      // collect information to construct 4-center builders
      auto basis = ao_factory_.basis_registry()->retrieve(OrbitalIndex(L"λ"));
      auto dcell = ao_factory_.unitcell().dcell();
      auto R_max = ao_factory_.R_max();
      auto RJ_max = ao_factory_.RJ_max();
      auto RD_max = ao_factory_.RD_max();
      auto R_size = ao_factory_.R_size();
      auto RJ_size = ao_factory_.RJ_size();
      auto RD_size = ao_factory_.RD_size();
      auto screen = ao_factory_.screen();
      auto screen_threshold = ao_factory_.screen_threshold();

      // construct PerioidcFourCenterFockBuilder for exchange term
      k_builder_ = std::make_unique<K_Builder>(
          world, basis, basis, dcell, R_max, RJ_max, RD_max, R_size, RJ_size,
          RD_size, false, true, screen, screen_threshold);
    }
    auto t1_k_init = mpqc::fenced_now(world);
    auto t_k_init = mpqc::duration_in_s(t0_k_init, t1_k_init);

    ExEnv::out0() << "\nInit RI-J time:      " << t_j_init << " s" << std::endl;
    ExEnv::out0() << "\nInit Four-Center-K time:      " << t_k_init << " s\n"
                  << std::endl;
  }

  ~PeriodicDFFockBuilder() {}

  array_type operator()(array_type const &D, double target_precision,
                        bool) override {
    array_type G;

    // the '-' sign is embeded in K builder
    G("mu, nu") = 2.0 * compute_J(D, target_precision)("mu, nu") +
                  compute_K(D, target_precision)("mu, nu");

    return G;
  }

  void register_fock(const array_type &fock,
                     FormulaRegistry<array_type> &registry) override {
    registry.insert(Formula(L"(κ|F|λ)"), fock);
  }

 private:
  Factory &ao_factory_;
  std::unique_ptr<J_Builder> j_builder_;
  std::unique_ptr<K_Builder> k_builder_;

 private:
  array_type compute_J(const array_type &D, double target_precision) {
    return j_builder_->operator()(D, target_precision);
  }

  array_type compute_K(const array_type &D, double target_precision) {
    return k_builder_->operator()(D, target_precision, false);
  }
};

}  // namespace scf
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_DF_FOCK_BUILDER_H_
