#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_MA_RI_J_BUILDER_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_MA_RI_J_BUILDER_H_

#include "mpqc/chemistry/qc/lcao/scf/builder.h"
#include "mpqc/chemistry/qc/lcao/scf/pbc/periodic_ma.h"
#include "mpqc/chemistry/qc/lcao/scf/pbc/periodic_ri_j_builder.h"
#include "mpqc/chemistry/qc/lcao/scf/pbc/util.h"

#include <boost/math/special_functions/legendre.hpp>

namespace mpqc {
namespace scf {

/*!
 * @brief PeriodicMARIJBuilder is a Coulomb matrix builder that uses
 * multipole-accelerated RI-J approximation, i.e. multipole approximation in
 * Crystal Far Field and RI-J in Crystal Near Field.
 * @tparam Tile
 * @tparam Policy
 * @tparam Factory
 */
template <typename Tile, typename Policy, typename Factory>
class PeriodicMARIJBuilder {
 public:
  using array_type = typename PeriodicFockBuilder<Tile, Policy>::array_type;
  using RIJ_Builder = PeriodicRIJBuilder<Tile, Policy, Factory>;
  using MA_Builder = ::mpqc::pbc::ma::PeriodicMA<Factory>;

  PeriodicMARIJBuilder(Factory &ao_factory) : ao_factory_(ao_factory) {
    auto &world = ao_factory_.world();

    // Construct multipole approximation builder
    auto t0_ma_init = mpqc::fenced_now(world);
    ma_builder_ = std::make_unique<MA_Builder>(ao_factory_);
    auto t1_ma_init = mpqc::fenced_now(world);
    auto t_ma_init = mpqc::duration_in_s(t0_ma_init, t1_ma_init);

    // set RJ_max_ to be the boundary of Crystal Near Field (= CFF_boundary - 1)
    // for the rest of the calculation
    const auto &cff_boundary = ma_builder_->CFF_boundary();
    Vector3i cnf_boundary = ao_factory_.RJ_max();
    for (auto dim = 0; dim <= 2; ++dim) {
      if (ma_builder_->CFF_reached(dim)) {
        cnf_boundary(dim) = cff_boundary(dim) - 1;
      }
    }
    ao_factory_.set_rjmax(cnf_boundary);

    ExEnv::out0() << "\nThe boundary of Crystal Near Field is "
                  << ao_factory_.RJ_max().transpose() << std::endl;

    // Construct periodic RI-J builder
    auto t0_j_init = mpqc::fenced_now(world);
    rij_builder_ = std::make_unique<RIJ_Builder>(ao_factory_);
    auto t1_j_init = mpqc::fenced_now(world);
    double t_j_init = mpqc::duration_in_s(t0_j_init, t1_j_init);

    ExEnv::out0() << "\nInit RI-J time:      " << t_j_init << " s" << std::endl;
    ExEnv::out0() << "\nInit MA time:        " << t_ma_init << " s\n"
                  << std::endl;
  }

  ~PeriodicMARIJBuilder() {}

  array_type operator()(array_type const &D, double target_precision) {
    if (ma_builder_->CFF_reached()) {
      compute_MAJ(D, target_precision);
    }

    return compute_RIJ(D, target_precision);
  }

  MA_Builder &multipole_builder() { return *ma_builder_; }

 private:
  Factory &ao_factory_;
  std::unique_ptr<RIJ_Builder> rij_builder_;
  std::unique_ptr<MA_Builder> ma_builder_;

 private:
  array_type compute_RIJ(const array_type &D, double target_precision) {
    return rij_builder_->operator()(D, target_precision);
  }

  void compute_MAJ(const array_type &D, double target_precision) {
    ma_builder_->compute_multipole_approx(D, target_precision);
  }
};

}  // namespace scf
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_MA_RI_J_BUILDER_H_
