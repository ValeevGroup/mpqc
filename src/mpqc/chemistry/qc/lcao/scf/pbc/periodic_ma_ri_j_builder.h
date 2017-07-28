#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_MA_RI_J_BUILDER_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_MA_RI_J_BUILDER_H_

#include "mpqc/chemistry/qc/lcao/scf/builder.h"
#include "mpqc/chemistry/qc/lcao/scf/pbc/periodic_mm.h"
#include "mpqc/chemistry/qc/lcao/scf/pbc/periodic_ri_j_builder.h"
#include "mpqc/chemistry/qc/lcao/scf/pbc/util.h"

#include <boost/math/special_functions/legendre.hpp>

namespace mpqc {
namespace scf {

template <typename Tile, typename Policy, typename Factory>
class PeriodicMARIJBuilder {
 public:
  using array_type = typename PeriodicFockBuilder<Tile, Policy>::array_type;
  using RIJ_Builder = PeriodicRIJBuilder<Tile, Policy, Factory>;
  using MM = ::mpqc::pbc::mm::PeriodicMM<Factory>;

  PeriodicMARIJBuilder(Factory &ao_factory) : ao_factory_(ao_factory) {
    auto &world = ao_factory_.world();

    // Construct periodic RI-J builder
    auto t0_j_init = mpqc::fenced_now(world);
    rij_builder_ = std::make_unique<RIJ_Builder>(ao_factory_);
    auto t1_j_init = mpqc::fenced_now(world);
    double t_j_init = mpqc::duration_in_s(t0_j_init, t1_j_init);

    // Construct exact perioic 4-center K builder
    auto t0_mm_init = mpqc::fenced_now(world);
    mm_ = std::make_unique<MM>(ao_factory_);
    auto t1_mm_init = mpqc::fenced_now(world);
    auto t_mm_init = mpqc::duration_in_s(t0_mm_init, t1_mm_init);

    double tmp1 = boost::math::legendre_p(0, 0.5);
    ExEnv::out0() << "\nLegendre P(0, 0.5) = " << tmp1 << std::endl;

    ExEnv::out0() << "\nInit RI-J time:      " << t_j_init << " s" << std::endl;
    ExEnv::out0() << "\nInit MM time:      " << t_mm_init << " s\n"
                  << std::endl;
  }

  ~PeriodicMARIJBuilder() {}

  array_type operator()(array_type const &D, double target_precision) {
    array_type G;

    G("mu, nu") = compute_RIJ(D, target_precision)("mu, nu");

    auto tmp = compute_MMJ(D, target_precision);

    return G;
  }

 private:
  Factory &ao_factory_;
  std::unique_ptr<RIJ_Builder> rij_builder_;
  std::unique_ptr<MM> mm_;

 private:
  array_type compute_RIJ(const array_type &D, double target_precision) {
    return rij_builder_->operator()(D, target_precision);
  }

  array_type compute_MMJ(const array_type &D, double target_precision) {
    array_type mmj;
    return mmj;
  }
};

}  // namespace scf
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_MA_RI_J_BUILDER_H_
