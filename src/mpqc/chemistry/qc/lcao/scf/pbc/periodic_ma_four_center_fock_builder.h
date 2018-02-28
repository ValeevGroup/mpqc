#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_MA_FOUR_CENTER_FOCK_BUILDER_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_MA_FOUR_CENTER_FOCK_BUILDER_H_

#include "mpqc/chemistry/qc/lcao/scf/builder.h"
#include "mpqc/chemistry/qc/lcao/scf/pbc/periodic_four_center_fock_builder.h"
#include "mpqc/chemistry/qc/lcao/scf/pbc/periodic_ma.h"

namespace mpqc {
namespace scf {

template <typename Tile, typename Policy>
class PeriodicMAFourCenterFockBuilder
    : public PeriodicFockBuilder<Tile, Policy> {
 public:
  using factory_type = ::mpqc::lcao::gaussian::PeriodicAOFactory<Tile, Policy>;
  using array_type = typename factory_type::TArray;
  using MA_Builder = ::mpqc::pbc::ma::PeriodicMA<factory_type>;
  using JK_Builder = PeriodicFourCenterFockBuilder<Tile, Policy>;

  PeriodicMAFourCenterFockBuilder(factory_type &factory)
      : ao_factory_(factory) {
    auto &world = ao_factory_.world();

    // Construct multipole approximation builder
    auto t0_ma_init = mpqc::fenced_now(world);
    ma_builder_ = std::make_unique<MA_Builder>(ao_factory_);
    auto t1_ma_init = mpqc::fenced_now(world);
    auto t_ma_init = mpqc::duration_in_s(t0_ma_init, t1_ma_init);

    // set RJ_max_ to be the boundary of Crystal Near Field (= CFF_boundary - 1)
    // for the rest of the calculation
    const auto &cff_boundary = ma_builder_->CFF_boundary();
    Vector3i cnf_boundary = {0, 0, 0};
    for (auto dim = 0; dim <= 2; ++dim) {
      if (cff_boundary(dim) > 0) {
        cnf_boundary(dim) = cff_boundary(dim) - 1;
      }
    }
    ao_factory_.set_rjmax(cnf_boundary);

    ExEnv::out0() << "\nThe boundary of Crystal Near Field is "
                  << ao_factory_.RJ_max().transpose() << std::endl;

    // Construct exact periodic 4-center J builder
    auto t0_j_init = mpqc::fenced_now(world);
    j_builder_ = std::make_unique<JK_Builder>(ao_factory_, true, false);
    auto t1_j_init = mpqc::fenced_now(world);
    auto t_j_init = mpqc::duration_in_s(t0_j_init, t1_j_init);

    // Construct exact periodic 4-center K builder
    auto t0_k_init = mpqc::fenced_now(world);
    k_builder_ = std::make_unique<JK_Builder>(ao_factory_, false, true);
    auto t1_k_init = mpqc::fenced_now(world);
    auto t_k_init = mpqc::duration_in_s(t0_k_init, t1_k_init);

    ExEnv::out0() << "\nInit MA time:            " << t_ma_init << " s"
                  << std::endl;
    ExEnv::out0() << "\nInit Four-Center-J time: " << t_j_init << " s"
                  << std::endl;
    ExEnv::out0() << "\nInit Four-Center-K time: " << t_k_init << " s\n"
                  << std::endl;
  }

  ~PeriodicMAFourCenterFockBuilder() {}

  array_type operator()(array_type const &D, double target_precision,
                        bool is_density_diagonal) override {

    const auto J = compute_J(D, target_precision, is_density_diagonal);
    const auto K = compute_K(D, target_precision, is_density_diagonal);
    const auto J_lattice_range = j_builder_->fock_lattice_range();
    const auto K_lattice_range = k_builder_->fock_lattice_range();

    // '2.0' and the '-' sign are embedded in J and K builders, respectively
    return ::mpqc::pbc::detail::add(J, K, J_lattice_range, K_lattice_range, 1.0, 1.0);
  }

  void register_fock(const array_type &fock,
                     FormulaRegistry<array_type> &registry) override {
    registry.insert(Formula(L"(κ|F|λ)"), fock);
  }

  Vector3i fock_lattice_range() override {
    const auto MA_lattice_range = j_builder_->fock_lattice_range();
    const auto JK_lattice_range = k_builder_->fock_lattice_range();
    if (MA_lattice_range(0) >= JK_lattice_range(0) &&
        MA_lattice_range(1) >= JK_lattice_range(1) &&
        MA_lattice_range(2) >= JK_lattice_range(2)) {
      return MA_lattice_range;
    } else if (MA_lattice_range(0) <= JK_lattice_range(0) &&
               MA_lattice_range(1) <= JK_lattice_range(1) &&
               MA_lattice_range(2) <= JK_lattice_range(2)) {
      return JK_lattice_range;
    } else {
      ExEnv::out0() << "\nLattice range of Coulomb: "
                    << MA_lattice_range.transpose()
                    << "\nLattice range of exchange: "
                    << JK_lattice_range.transpose() << std::endl;
      throw "Invalid lattice ranges!";
    }
  }

  MA_Builder &multipole_builder() {
    return *ma_builder_;
  }

 private:
  factory_type &ao_factory_;
  std::unique_ptr<MA_Builder> ma_builder_;
  std::unique_ptr<JK_Builder> j_builder_;
  std::unique_ptr<JK_Builder> k_builder_;

  array_type compute_J(const array_type &D, double target_precision, bool is_density_diagonal) {
    ma_builder_->compute_multipole_approx(D, target_precision);

    return j_builder_->operator()(D, target_precision, is_density_diagonal);
  }

  array_type compute_K(const array_type &D, double target_precision, bool is_density_diagonal) {
    return k_builder_->operator()(D, target_precision, is_density_diagonal);
  }

};

}  // namespace scf
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_MA_FOUR_CENTER_FOCK_BUILDER_H_
