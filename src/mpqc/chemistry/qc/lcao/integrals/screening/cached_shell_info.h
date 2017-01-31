
#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_SCREENING_CACHED_SHELL_INFO_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_SCREENING_CACHED_SHELL_INFO_H_


#include "mpqc/chemistry/qc/lcao/integrals/task_integrals_common.h"

namespace mpqc {
namespace lcao {
namespace gaussian {
namespace detail {

/*! \brief Class to hold information for QQR type screening.
 */
class CachedShellInfo {
  private:
    RowMatrixXd pair_extents_;
    std::vector<std::vector<Vector3d>> pair_centers_;
    double erfinv_thr_;

  public:
    CachedShellInfo() = default;
    CachedShellInfo(ShellVec const &shells, double thresh);

    Vector3d center(int64_t i, int64_t j) const { return pair_centers_[i][j]; }
    double extent(int64_t i, int64_t j) const { return pair_extents_(i, j); }

  private:
    double pair_extent(Shell const &sh0, Shell const &sh1, Vector3d const &r_01); 
    Vector3d shell_weighted_center(Shell const &sh0, Shell const &sh1);
};

}  // namespace  detail
}  // namespace  gaussian
}  // namespace  lcao
}  // namespace  mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_SCREENING_CACHED_SHELL_INFO_H_
