#pragma once
#ifndef MPQC_INTEGRALS_SCREENING_QVLSHELLINFO_H
#define MPQC_INTEGRALS_SCREENING_QVLSHELLINFO_H

#include "../../common/typedefs.h"
#include "../task_integrals_common.h"

namespace mpqc {
namespace integrals {
namespace detail {

/*! \brief Class to hold information for QVl aux basis type screening.
 */
class QVlShellInfo {
  private:
    std::vector<double> a_extents_;
    MatrixD cd_extents_;

    std::vector<Vec3D> a_centers_;
    std::vector<std::vector<Vec3D>> cd_centers_;

    std::vector<double> a_factor_;
    std::vector<double> a_am_;
    MatrixD one_over_gamma_ab_pow_;

    double erfinv_thr_;

  public:
    QVlShellInfo() = default;
    QVlShellInfo(ShellVec const &a_shells, ShellVec const &cd_shells,
                 double thresh);

    double a_extent(int64_t i) const { return a_extents_[i]; }
    Vec3D const &a_center(int64_t i) const { return a_centers_[i]; }

    double cd_extent(int64_t i, int64_t j) const { return cd_extents_(i, j); }
    Vec3D const &cd_center(int64_t i, int64_t j) const {
        return cd_centers_[i][j];
    }

    // Returns the factor that scales Qab in David's screening paper.
    double
    scaling_factor(int64_t a, int64_t c, int64_t d, double Rnorm) const {
        const auto rscale = 1 / (std::pow(Rnorm, a_am_[a] + 1));
        return rscale * a_factor_[a] * one_over_gamma_ab_pow_(c,d);
    }

  private:
    double pair_extent(Shell const &sh0, Shell const &sh1, Vec3D const &r_01);
    double a_shell_min_exp(Shell const &sh);
    Vec3D shell_weighted_center(Shell const &sh0, Shell const &sh1) const;
};

} // namespace detail
} // namespace integrals
} // namespace mpqc

#endif // MPQC_INTEGRALS_SCREENING_QVLSHELLINFO_H
