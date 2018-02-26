
#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PURIFICATION_DENSITY_BUILD_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PURIFICATION_DENSITY_BUILD_H_

#include <vector>
#include "mpqc/chemistry/qc/lcao/scf/density_builder.h"

namespace mpqc {
namespace scf {

template <typename Tile, typename Policy>
class PurificationDensityBuilder : public DensityBuilder<Tile, Policy> {
 public:
  using array_type = typename DensityBuilder<Tile, Policy>::array_type;

 private:
  array_type S_;
  array_type M_inv_;
  array_type I_;
  std::vector<array_type> r_xyz_ints_;

  double TcutC_;
  bool localize_;
  std::string localization_method_;
  int64_t n_coeff_clusters_;
  bool clustered_coeffs_;
  int64_t occ_;
  int64_t ncore_;

 public:
  /// @param[in] localization_method defined the localization method; valid
  ///            choices are "boys-foster" (default),
  ///            "boys-foster(valence)" (do not localize the core).
  PurificationDensityBuilder(array_type const &S, std::vector<array_type> r_xyz,
                             int64_t occ, int64_t ncore, int64_t nclusters,
                             double TcutC = 0.0, bool localize = true,
                             std::string localization_method = "boys-foster",
                             bool clustered_coeffs = false);

  std::pair<array_type, array_type> operator()(array_type const &F) override;

  inline void print_iter(std::string const &) override {}

 private:
  array_type purify(array_type const &);
  array_type orbitals(array_type const &);
};

}  // namespace scf
}  // namespace mpqc

#include "purification_density_build_impl.h"

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PURIFICATION_DENSITY_BUILD_H_
