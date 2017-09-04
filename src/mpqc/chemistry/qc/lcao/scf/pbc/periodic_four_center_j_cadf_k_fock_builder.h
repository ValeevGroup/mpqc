#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_FOUR_CENTER_J_CADF_K_FOCK_BUILDER_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_FOUR_CENTER_J_CADF_K_FOCK_BUILDER_H_

#include "mpqc/chemistry/qc/lcao/scf/builder.h"
#include "mpqc/chemistry/qc/lcao/scf/pbc/periodic_cadf_k_builder.h"
#include "mpqc/chemistry/qc/lcao/scf/pbc/periodic_four_center_fock_builder.h"

namespace mpqc {
namespace scf {

template <typename Tile, typename Policy, typename Factory>
class PeriodicFourCenterJCADFKFockBuilder
    : public PeriodicFockBuilder<Tile, Policy> {
 public:
  using array_type = typename PeriodicFockBuilder<Tile, Policy>::array_type;
  using J_Builder = PeriodicFourCenterFockBuilder<Tile, Policy>;
  using K_Builder = PeriodicCADFKBuilder<Tile, Policy, Factory>;

  PeriodicFourCenterJCADFKFockBuilder(Factory &ao_factory,
                                      const double force_shape_threshold = 0.0)
      : ao_factory_(ao_factory) {
    auto &world = ao_factory_.world();

    // Construct periodic four-center J builder
    auto t0_j_init = mpqc::fenced_now(world);
    j_builder_ = std::make_unique<J_Builder>(ao_factory_, true, false);
    auto t1_j_init = mpqc::fenced_now(world);
    double t_j_init = mpqc::duration_in_s(t0_j_init, t1_j_init);

    // Construct periodic CADF-K builder
    auto t0_k_init = mpqc::fenced_now(world);
    k_builder_ =
        std::make_unique<K_Builder>(world, ao_factory_, force_shape_threshold);
    auto t1_k_init = mpqc::fenced_now(world);
    auto t_k_init = mpqc::duration_in_s(t0_k_init, t1_k_init);

    ExEnv::out0() << "\nInit FourCenter-J time: " << t_j_init << " s\n"
                  << "\nInit CADF-K time:       " << t_k_init << " s\n"
                  << std::endl;
  }

  ~PeriodicFourCenterJCADFKFockBuilder() {}

  array_type operator()(array_type const &D, double target_precision,
                        bool) override {
    array_type G;

    G("mu, nu") = compute_J(D, target_precision)("mu, nu") -
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
    return j_builder_->operator()(D, target_precision, false);
  }

  array_type compute_K(const array_type &D, double target_precision) {
    return k_builder_->operator()(D, target_precision);
  }
};

}  // namespace scf
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_FOUR_CENTER_J_CADF_K_FOCK_BUILDER_H_
