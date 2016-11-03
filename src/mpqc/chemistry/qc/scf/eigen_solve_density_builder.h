
#ifndef MPQC_SCF_EIGENSOLVEDENSITYBUILDER_H
#define MPQC_SCF_EIGENSOLVEDENSITYBUILDER_H

#include <mpqc/chemistry/qc/scf/density_builder.h>
#include <vector>
#include <array>

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

    double inverse_time_;
    std::vector<double> esolve_times_;
    std::vector<double> localization_times_;
    std::vector<double> clustering_times_;

    std::vector<std::array<double, 3>> coeff_storages_;
    std::vector<std::array<double, 3>> density_storages_;

  public:
    ESolveDensityBuilder() = default;
    ESolveDensityBuilder(ESolveDensityBuilder const &) = default;
    ESolveDensityBuilder(ESolveDensityBuilder &&) = default;
    ESolveDensityBuilder &operator=(ESolveDensityBuilder const &) = default;
    ESolveDensityBuilder &operator=(ESolveDensityBuilder &&) = default;
    ~ESolveDensityBuilder() noexcept = default;

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

    inline double TcutC() const { return TcutC_; }
};

} // namespace scf
} // namespace mpqc
#endif //  MPQC_SCF_EIGENSOLVEDENSITYBUILDER_H
