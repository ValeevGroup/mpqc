#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_MA_RI_J_FOUR_CENTER_K_FOCK_BUILDER_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_MA_RI_J_FOUR_CENTER_K_FOCK_BUILDER_H_

#include "mpqc/chemistry/qc/lcao/scf/builder.h"
#include "mpqc/chemistry/qc/lcao/scf/pbc/periodic_four_center_fock_builder.h"
#include "mpqc/chemistry/qc/lcao/scf/pbc/periodic_ma_ri_j_builder.h"

namespace mpqc {
namespace scf {

/*!
 * @brief PeriodicMARIJFourCenterKFockBuilder is an implementation of
 * PeriodicFockBuilder with multipole-accelerated RI-J and 4-center exchange.
 * For Coulomb, multipole approximation is used in Crystal Far Field and RI-J in
 * Crystal Near Field. For exchange, 4-center-K is always used.
 * @tparam Tile
 * @tparam Policy
 */
template <typename Tile, typename Policy>
class PeriodicMARIJFourCenterKFockBuilder
    : public PeriodicFockBuilder<Tile, Policy> {
 public:
  using factory_type = ::mpqc::lcao::gaussian::PeriodicAOFactory<Tile, Policy>;
  using array_type = typename factory_type::TArray;
  using J_Builder = PeriodicMARIJBuilder<Tile, Policy, factory_type>;
  using K_Builder = PeriodicFourCenterFockBuilder<Tile, Policy>;

  PeriodicMARIJFourCenterKFockBuilder(factory_type &ao_factory)
      : ao_factory_(ao_factory) {
    auto &world = ao_factory_.world();

    // Construct periodic MA-RI-J builder
    auto t0_j_init = mpqc::fenced_now(world);
    j_builder_ = std::make_unique<J_Builder>(ao_factory_);
    auto t1_j_init = mpqc::fenced_now(world);
    double t_j_init = mpqc::duration_in_s(t0_j_init, t1_j_init);

    // Construct exact perioic 4-center K builder
    auto t0_k_init = mpqc::fenced_now(world);
    k_builder_ = std::make_unique<K_Builder>(ao_factory_, false, true);
    auto t1_k_init = mpqc::fenced_now(world);
    auto t_k_init = mpqc::duration_in_s(t0_k_init, t1_k_init);

    ExEnv::out0() << "\nInit MA-RI-J time:     " << t_j_init << " s"
                  << std::endl;
    ExEnv::out0() << "\nInit Four-Center-K time:      " << t_k_init << " s\n"
                  << std::endl;
  }

  ~PeriodicMARIJFourCenterKFockBuilder() {}

  array_type operator()(array_type const &D, double target_precision,
                        bool) override {
    array_type G;

    const auto J = compute_J(D, target_precision);
    const auto K = compute_K(D, target_precision);
    const auto J_lattice_range = ao_factory_.R_max();
    const auto K_lattice_range = k_builder_->fock_lattice_range();
    // the '-' sign is embedded in K builder
    G = ::mpqc::pbc::detail::add(J, K, J_lattice_range, K_lattice_range, 2.0,
                                 1.0);

    return G;
  }

  void register_fock(const array_type &fock,
                     FormulaRegistry<array_type> &registry) override {
    registry.insert(Formula(L"(κ|F|λ)"), fock);
  }

  Vector3i fock_lattice_range() override {
    const auto J_lattice_range = ao_factory_.R_max();
    const auto K_lattice_range = k_builder_->fock_lattice_range();
    if (J_lattice_range(0) >= K_lattice_range(0) &&
        J_lattice_range(1) >= K_lattice_range(1) &&
        J_lattice_range(2) >= K_lattice_range(2)) {
      return J_lattice_range;
    } else if (J_lattice_range(0) <= K_lattice_range(0) &&
               J_lattice_range(1) <= K_lattice_range(1) &&
               J_lattice_range(2) <= K_lattice_range(2)) {
      return K_lattice_range;
    } else {
      ExEnv::out0() << "\nLattice range of Coulomb: "
                    << J_lattice_range.transpose()
                    << "\nLattice range of exchange: "
                    << K_lattice_range.transpose() << std::endl;
      throw "Invalid lattice ranges!";
    }
  }

  J_Builder &coulomb_builder() { return *j_builder_; }

 private:
  factory_type &ao_factory_;
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

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_MA_RI_J_FOUR_CENTER_K_FOCK_BUILDER_H_
