
#ifndef MPQC_SCF_DENSITYBUILDER_H
#define MPQC_SCF_DENSITYBUILDER_H

#include <string>
#include <utility>

#include <tiledarray.h>

namespace mpqc {
namespace scf {

namespace TA = TiledArray;

class DensityBuilder {
  public:
    using array_type = TA::TSpArrayD;
    DensityBuilder() = default;
    DensityBuilder(DensityBuilder const &) = default;
    DensityBuilder(DensityBuilder &&) = default;
    DensityBuilder &operator=(DensityBuilder const &) = default;
    DensityBuilder &operator=(DensityBuilder &&) = default;
    virtual ~DensityBuilder() = default;

    virtual std::pair<array_type, array_type> operator()(array_type const &)
          = 0;

    virtual void print_iter(std::string const &) = 0;

};

} // namespace scf
} // namespace mpqc
#endif //  MPQC_SCF_DENSITYBUILDER_H
