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
    k_builder_ = std::make_unique<K_Builder>(ao_factory_, false, true);
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

    const auto J_latt_range = ao_factory_.R_max();
    const auto K_latt_range = k_builder_->fock_latt_range();
    // the '-' sign is embeded in K builder
    G = ::mpqc::pbc::detail::add(compute_J(D, target_precision),
                                 compute_K(D, target_precision), J_latt_range,
                                 K_latt_range, 2.0, 1.0);

    return G;
  }

  void register_fock(const array_type &fock,
                     FormulaRegistry<array_type> &registry) override {
    registry.insert(Formula(L"(κ|F|λ)"), fock);
  }

  Vector3i fock_latt_range() override {
    const auto J_latt_range = ao_factory_.R_max();
    const auto K_latt_range = k_builder_->fock_latt_range();
    if (J_latt_range(0) >= K_latt_range(0) &&
        J_latt_range(1) >= K_latt_range(1) &&
        J_latt_range(2) >= K_latt_range(2)) {
      return J_latt_range;
    } else if (J_latt_range(0) <= K_latt_range(0) &&
               J_latt_range(1) <= K_latt_range(1) &&
               J_latt_range(2) <= K_latt_range(2)) {
      return K_latt_range;
    } else {
      ExEnv::out0() << "\nLattice range of Coulomb: "
                    << J_latt_range.transpose()
                    << "\nLattice range of exchange: "
                    << K_latt_range.transpose() << std::endl;
      throw "Invalid lattice ranges!";
    }
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
