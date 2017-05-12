#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_CADF_K_BUILDER_H_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_CADF_K_BUILDER_H_H_

#include "mpqc/chemistry/qc/lcao/scf/builder.h"
#include "mpqc/chemistry/qc/lcao/scf/pbc/periodic_three_center_contraction_builder.h"

namespace mpqc {
namespace scf {

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
};

}  // namespace scf
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_CADF_K_BUILDER_H_H_
