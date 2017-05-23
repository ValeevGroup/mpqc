#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_CADF_K_BUILDER_H_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_CADF_K_BUILDER_H_H_

#include "mpqc/chemistry/qc/lcao/integrals/density_fitting/cadf_coeffs.h"
#include "mpqc/chemistry/qc/lcao/scf/builder.h"
#include "mpqc/chemistry/qc/lcao/scf/pbc/periodic_three_center_contraction_builder.h"

namespace mpqc {
namespace scf {

namespace detail {}  // namespace detail

template <typename Tile, typename Policy, typename Factory>
class PeriodicCADFKBuilder {
 public:
  using array_type = TA::DistArray<Tile, Policy>;
  using PTC_Builder = PeriodicThreeCenterContractionBuilder<Tile, Policy>;

  PeriodicCADFKBuilder(Factory &ao_factory) : ao_factory_(ao_factory) {
    auto &world = ao_factory_.world();

    // by-cluster orbital basis and df basis
    auto obs = ao_factory_.basis_registry()->retrieve(OrbitalIndex(L"λ"));
    auto dfbs = ao_factory_.basis_registry()->retrieve(OrbitalIndex(L"Κ"));

    auto dcell = ao_factory_.unitcell().dcell();
    auto R_max = ao_factory_.R_max();
    auto RJ_max = ao_factory_.RJ_max();
    auto RD_max = ao_factory_.RD_max();
    auto R_size = ao_factory_.R_size();
    auto RJ_size = ao_factory_.RJ_size();
    auto RD_size = ao_factory_.RD_size();

    const Vector3i ref_latt_range = {0, 0, 0};
    const auto natoms_per_uc = ao_factory_.unitcell().natoms();
    using ::mpqc::lcao::gaussian::detail::shift_basis_origin;
    Vector3d zero_shift_base(0.0, 0.0, 0.0);
    auto bs1 = shift_basis_origin(*obs, zero_shift_base, RJ_max, dcell);
    auto dfbs1 = shift_basis_origin(*dfbs, zero_shift_base, RJ_max, dcell);
    C_ = lcao::cadf_fitting_coefficients<Tile, Policy>(
        world, *obs, *bs1, *dfbs1, ref_latt_range, RJ_max, RJ_max,
        natoms_per_uc);

    ExEnv::out0() << "\nC_ = \n" << C_ << std::endl;
  }

  ~PeriodicCADFKBuilder() {}

  array_type operator()(array_type const &D, double target_precision) {
    array_type K;
    // to be implemented
    K("mu, nu") = 1.0 * D("mu, nu");
    return K;
  }

 private:
  Factory &ao_factory_;
  bool print_detail_;
  std::unique_ptr<PTC_Builder> three_center_builder_;

  array_type C_;
};

}  // namespace scf
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_CADF_K_BUILDER_H_H_
