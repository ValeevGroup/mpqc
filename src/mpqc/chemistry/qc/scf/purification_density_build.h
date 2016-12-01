
#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PURIFICATION_DENSITY_BUILD_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PURIFICATION_DENSITY_BUILD_H_

#include "mpqc/chemistry/qc/scf/density_builder.h"
#include <vector>

namespace mpqc {
namespace scf {

class PurificationDensityBuilder : public DensityBuilder {
 private:
  array_type S_;
  array_type M_inv_;
  array_type I_;
  std::vector<array_type> r_xyz_ints_;

  double TcutC_;
  bool localize_;
  int64_t n_coeff_clusters_;
  int64_t occ_;

 public:
  PurificationDensityBuilder(array_type const &S, std::vector<array_type> r_xyz,
                             int64_t occ, int64_t nclusters, double TcutC = 0.0,
                             bool localize = true);

  std::pair<array_type, array_type> operator()(array_type const &F) override;

  inline void print_iter(std::string const &) override {}

 private:
  array_type purify(array_type const &);
  array_type orbitals(array_type const &);
};

}  // namespace scf
}  // namespace mpqc
#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PURIFICATION_DENSITY_BUILD_H_
