#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_RI_J_CADF_K_FOCK_BUILDER_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_RI_J_CADF_K_FOCK_BUILDER_H_

#include "mpqc/chemistry/qc/lcao/scf/builder.h"
#include "mpqc/chemistry/qc/lcao/scf/pbc/periodic_cadf_k_builder.h"
#include "mpqc/chemistry/qc/lcao/scf/pbc/periodic_ri_j_builder.h"
#include "mpqc/chemistry/qc/lcao/scf/pbc/util.h"

namespace mpqc {
namespace lcao {
namespace scf {

template<typename Tile, typename Policy, typename Factory>
class PeriodicRIJCADFKFockBuilder : public PeriodicFockBuilder<Tile, Policy> {
 public:
  using array_type = typename PeriodicFockBuilder<Tile, Policy>::array_type;
  using J_Builder = PeriodicRIJBuilder<Tile, Policy, Factory>;
  using K_Builder = PeriodicCADFKBuilder<Tile, Policy, Factory>;

  PeriodicRIJCADFKFockBuilder(Factory &ao_factory) : ao_factory_(ao_factory) {
    auto &world = ao_factory_.world();

    // Construct periodic RI-J builder
    auto t0_j_init = mpqc::fenced_now(world);
    j_builder_ = std::make_unique<J_Builder>(ao_factory_);
    auto t1_j_init = mpqc::fenced_now(world);
    double t_j_init = mpqc::duration_in_s(t0_j_init, t1_j_init);

    // Construct periodic CADF-K builder
    auto t0_k_init = mpqc::fenced_now(world);
    k_builder_ = std::make_unique<K_Builder>(ao_factory_);
    auto t1_k_init = mpqc::fenced_now(world);
    auto t_k_init = mpqc::duration_in_s(t0_k_init, t1_k_init);

    ExEnv::out0() << "\nInit RI-J time:      " << t_j_init << " s" << std::endl;
    ExEnv::out0() << "\nInit CADF-K time:      " << t_k_init << " s\n"
                  << std::endl;
  }

  ~PeriodicRIJCADFKFockBuilder() {}

  array_type operator()(array_type const &D, double target_precision,
                        bool) override {
    array_type G;

    const auto J = compute_J(D, target_precision);
    const auto K = compute_K(D, target_precision);
    const auto J_lattice_range = ao_factory_.R_max();
    const auto K_lattice_range = k_builder_->K_lattice_range();
    G = ::mpqc::pbc::detail::add(J, K, J_lattice_range, K_lattice_range, 2.0,
                                 -1.0);

    return G;
  }

  void register_fock(const array_type &fock,
                     FormulaRegistry<array_type> &registry) override {
    registry.insert(Formula(L"(κ|F|λ)"), fock);
  }

  Vector3i fock_lattice_range() override {
    const auto J_range = ao_factory_.R_max();
    const auto K_range = k_builder_->K_lattice_range();
    if (J_range(0) <= K_range(0) && J_range(1) <= K_range(1) &&
        J_range(2) <= K_range(2)) {
      return K_range;
    } else if (J_range(0) >= K_range(0) && J_range(1) >= K_range(1) &&
        J_range(2) >= K_range(2)) {
      return J_range;
    } else {
      throw "invalid lattice ranges";
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
    return k_builder_->operator()(D, target_precision);
  }
};

}  // namespace scf
}  // namespace lcao
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_RI_J_CADF_K_FOCK_BUILDER_H_
