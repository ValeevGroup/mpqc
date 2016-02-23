#pragma once
#ifndef MPQC_SCF_EIGENSOLVEDENSITYBUILDER_H
#define MPQC_SCF_EIGENSOLVEDENSITYBUILDER_H

#include "density_builder.h"
#include <vector>

namespace mpqc {
namespace scf {

class ESolveDensityBuilder : public DensityBuilder {
  private:
    array_type S_;
    array_type M_inv_;
    std::vector<array_type> r_xyz_ints_;

    double TcutC_;
    bool localize_;
    int64_t n_coeff_clusters_;
    std::string metric_decomp_type_;
    int64_t occ_;

    double condition_num_thresh_ = 1e-10;

  public:
    ESolveDensityBuilder() = default;
    ESolveDensityBuilder(ESolveDensityBuilder const &) = default;
    ESolveDensityBuilder(ESolveDensityBuilder &&) = default;
    ESolveDensityBuilder &operator=(ESolveDensityBuilder const &) = default;
    ESolveDensityBuilder &operator=(ESolveDensityBuilder &&) = default;
    ~ESolveDensityBuilder() = default;

    ESolveDensityBuilder(array_type const &S, std::vector<array_type> r_xyz,
                         int64_t occ, int64_t nclusters, double TcutC = 0.0,
                         std::string const &metric_decomp_type
                         = "cholesky inverse",
                         bool localize = true);

    std::pair<array_type, array_type> operator()(array_type const &F) override;

    inline void print_iter(std::string const &) override {}

    inline double condition_num_threshold() const {
        return condition_num_thresh_;
    }
    inline void condition_num_threshold(double thresh) {
        condition_num_thresh_ = thresh;
    }
};

} // namespace scf
} // namespace mpqc
#endif //  MPQC_SCF_EIGENSOLVEDENSITYBUILDER_H
