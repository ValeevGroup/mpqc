#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_RI_J_CADF_K_FOCK_BUILDER_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_RI_J_CADF_K_FOCK_BUILDER_H_

#include "mpqc/chemistry/qc/lcao/scf/builder.h"
#include "mpqc/chemistry/qc/lcao/scf/pbc/periodic_cadf_k_builder.h"
#include "mpqc/chemistry/qc/lcao/scf/pbc/periodic_ri_j_builder.h"

namespace mpqc {
namespace scf {

template <typename Tile, typename Policy, typename Factory>
class PeriodicRIJCADFKFockBuilder : public PeriodicFockBuilder<Tile, Policy> {
 public:
  using array_type = typename PeriodicFockBuilder<Tile, Policy>::array_type;
  using J_Builder = PeriodicRIJBuilder<Tile, Policy, Factory>;
  using K_Builder = PeriodicCADFKBuilder<Tile, Policy, Factory>;

  PeriodicRIJCADFKFockBuilder(Factory &ao_factory,
                              const double force_shape_threshold = 0.0)
      : ao_factory_(ao_factory) {
    auto &world = ao_factory_.world();

    // Construct periodic RI-J builder
    auto t0_j_init = mpqc::fenced_now(world);
    j_builder_ = std::make_unique<J_Builder>(ao_factory_);
    auto t1_j_init = mpqc::fenced_now(world);
    double t_j_init = mpqc::duration_in_s(t0_j_init, t1_j_init);

    // Construct periodic CADF-K builder
    auto t0_k_init = mpqc::fenced_now(world);
    {
      auto obs = ao_factory_.basis_registry()->retrieve(OrbitalIndex(L"λ"));
      auto unitcell = ao_factory_.unitcell();
      auto dfbs_str = ao_factory_.df_k_basis();
      ExEnv::out0() << "\ndf_k_basis: " << dfbs_str << std::endl;
      auto dfbs = std::make_shared<lcao::gaussian::Basis>(
          lcao::gaussian::parallel_make_basis(
              world, lcao::gaussian::Basis::Factory(dfbs_str), unitcell));
      auto dcell = ao_factory_.unitcell().dcell();
      auto R_max = ao_factory_.R_max();
      auto RJ_max = ao_factory_.RJ_max();
      auto RD_max = ao_factory_.RD_max();
      auto R_size = ao_factory_.R_size();
      auto RJ_size = ao_factory_.RJ_size();
      auto RD_size = ao_factory_.RD_size();
      auto screen_threshold = ao_factory_.screen_threshold();
      auto shell_pair_threshold = ao_factory_.shell_pair_threshold();
      auto density_threshold = ao_factory_.density_threshold();
      auto target_precision = std::numeric_limits<double>::epsilon();
      auto ntiles_per_uc = obs->nclusters();
      auto natoms_per_uc = ao_factory_.unitcell().natoms();
      auto print_detail = ao_factory_.print_detail();
      k_builder_ = std::make_unique<K_Builder>(
          world, ao_factory_, obs, dfbs, dcell, R_max, RJ_max, RD_max, R_size,
          RJ_size, RD_size, ntiles_per_uc, natoms_per_uc, shell_pair_threshold,
          screen_threshold, density_threshold, target_precision, print_detail,
          force_shape_threshold);
    }
    //    k_builder_ =
    //        std::make_unique<K_Builder>(world, ao_factory_,
    //        force_shape_threshold);
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

    G("mu, nu") = 2.0 * compute_J(D, target_precision)("mu, nu") -
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
    return k_builder_->operator()(D, target_precision);
  }
};

}  // namespace scf
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_RI_J_CADF_K_FOCK_BUILDER_H_
