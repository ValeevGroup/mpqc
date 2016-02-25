#pragma once
#ifndef MPQC_SCF_PURIFICATIONDENSITYBUILD_H
#define MPQC_SCF_PURIFICATIONDENSITYBUILD_H

#include "density_builder.h"
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
    PurificationDensityBuilder(array_type const &S,
                               std::vector<array_type> r_xyz, int64_t occ,
                               int64_t nclusters, double TcutC = 0.0,
                               bool localize = true);

    std::pair<array_type, array_type> operator()(array_type const &F) override;

    inline void print_iter(std::string const &) override {}
    rapidjson::Value results(rapidjson::Document &d) override;

  private:
    array_type purify(array_type const &);
    array_type orbitals(array_type const &);
};

} // namespace scf
} // namespace mpqc
#endif // MPQC_SCF_PURIFICATIONDENSITYBUILD_H
