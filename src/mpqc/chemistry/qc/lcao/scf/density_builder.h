
#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_DENSITY_BUILDER_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_DENSITY_BUILDER_H_

#include <string>
#include <utility>

#include <tiledarray.h>

namespace mpqc {
namespace scf {

namespace TA = TiledArray;

template <typename Tile, typename Policy>
class DensityBuilder {
 public:
  using array_type = TA::DistArray<Tile,Policy>;
  DensityBuilder() = default;
  DensityBuilder(DensityBuilder const &) = default;
  DensityBuilder(DensityBuilder &&) = default;
  DensityBuilder &operator=(DensityBuilder const &) = default;
  DensityBuilder &operator=(DensityBuilder &&) = default;
  virtual ~DensityBuilder() { }

  virtual std::pair<array_type, array_type> operator()(array_type const &) = 0;

  virtual void print_iter(std::string const &) = 0;
};

}  // namespace scf
}  // namespace mpqc
#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_DENSITY_BUILDER_H_
