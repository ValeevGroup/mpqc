#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_RI_J_CADF_K_FOCK_BUILDER_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_RI_J_CADF_K_FOCK_BUILDER_H_

#include "mpqc/chemistry/qc/lcao/scf/builder.h"
#include "mpqc/chemistry/qc/lcao/scf/pbc/periodic_ri_j_builder.h"

namespace mpqc {
namespace scf {

template <typename Tile, typename Policy, typename Factory>
class PeriodicRIJCADFKFockBuilder : public PeriodicFockBuilder<Tile, Policy> {
public:
  using array_type = typename PeriodicFockBuilder<Tile, Policy>::array_type;
  using J_Builder = PeriodicRIJBuilder<Tile, Policy, Factory>;

  PeriodicRIJCADFKFockBuilder(Factory &ao_factory) : ao_factory_(ao_factory) {
    auto &world = ao_factory_.world();

    // Construct periodic RI-J builder
    auto t0_j_init = mpqc::fenced_now(world);
    j_builder_ = std::make_unique<J_Builder>(ao_factory_);
    auto t1_j_init = mpqc::fenced_now(world);
    double t_j_init = mpqc::duration_in_s(t0_j_init, t1_j_init);

    // Construct periodic CADF-K builder
    auto t0_k_init = mpqc::fenced_now(world);
    {
      // TODO
    }
    auto t1_k_init = mpqc::fenced_now(world);
    auto t_k_init = mpqc::duration_in_s(t0_k_init, t1_k_init);

    ExEnv::out0() << "\nInit RI-J time:      " << t_j_init << " s" << std::endl;
    ExEnv::out0() << "\nInit CADF-K time:      " << t_k_init << " s\n"
                  << std::endl;
  }

  ~PeriodicRIJCADFKFockBuilder() {}

  array_type operator()(array_type const &D, double target_precision, bool) override {
    array_type G;

    G("mu, nu") = 2.0 * compute_J(D, target_precision)("mu, nu");

    return G;
  }

  void register_fock(const array_type &fock,
                     FormulaRegistry<array_type> &registry) override {
    registry.insert(Formula(L"(κ|F|λ)"), fock);
  }

private:
 Factory &ao_factory_;
 std::unique_ptr<J_Builder> j_builder_;

private:
 array_type compute_J(const array_type &D, double target_precision) {
   return j_builder_->operator()(D, target_precision);
 }

 array_type compute_K(const array_type &D, double target_precision) {
   array_type K;
   return K;
 }

};

}  // namespace scf
}  // namespace mpqc

#endif // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_RI_J_CADF_K_FOCK_BUILDER_H_
