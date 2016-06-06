#pragma once
#ifndef MPQC_INTEGRALS_SCREENING_CACHEDSHELLINFO_H
#define MPQC_INTEGRALS_SCREENING_CACHEDSHELLINFO_H

#include "../../../../../../common/typedefs.h"
#include <mpqc/chemistry/qc/integrals/task_integrals_common.h>

namespace mpqc {
namespace integrals {
namespace detail {

/*! \brief Class to hold information for QQR type screening.
 */
class CachedShellInfo {
  private:
    MatrixD pair_extents_;
    std::vector<std::vector<Vec3D>> pair_centers_;
    double erfinv_thr_;

  public:
    CachedShellInfo() = default;
    CachedShellInfo(ShellVec const &shells, double thresh);

    Vec3D center(int64_t i, int64_t j) const { return pair_centers_[i][j]; }
    double extent(int64_t i, int64_t j) const { return pair_extents_(i, j); }

  private:
    double pair_extent(Shell const &sh0, Shell const &sh1, Vec3D const &r_01); 
    Vec3D shell_weighted_center(Shell const &sh0, Shell const &sh1);
};

} // namespace detail
} // namespace integrals
} // namespace mpqc

#endif // MPQC_INTEGRALS_SCREENING_CACHEDSHELLINFO_H
