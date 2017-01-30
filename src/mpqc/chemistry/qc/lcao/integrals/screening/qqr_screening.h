
#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_SCREENING_QQR_SCREENING_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_SCREENING_QQR_SCREENING_H_


#include "mpqc/chemistry/qc/lcao/integrals/task_integrals_common.h"

#include "mpqc/chemistry/qc/lcao/integrals/screening/screen_base.h"
#include "mpqc/chemistry/qc/lcao/integrals/screening/schwarz_screen.h"
#include "mpqc/chemistry/qc/lcao/integrals/screening/cached_shell_info.h"

#include <boost/math/special_functions/erf.hpp>

namespace mpqc {
namespace lcao {
namespace gaussian {

/*! \brief Class to implement QQR Screening.
 *
 * Based on: Maurer, S. A.; Lambrecht, D. S.; Flaig, D. JCP. 136, 144107 (2012).
 *
 * Reverse Engineered from David Hollman's code in MPQC2/3
 *
 * QQR currently doesn't support mixed basis integrals
 *
 */
class QQR : public SchwarzScreen {
  private:
    detail::CachedShellInfo shell_info_ab_;
    static double well_seperated_thresh_;

  public:
    QQR() = default;
    QQR(SchwarzScreen ss, ShellVec const &shells)
            : SchwarzScreen(std::move(ss)),
              shell_info_ab_(shells, erfinv_ws_thresh()) {}

    static double well_sep_threshold() { return well_seperated_thresh_; }
    static void well_sep_threshold(double new_val) {
        well_seperated_thresh_ = new_val;
    }

    /// computes erfc^{-1} of the well_seperated_thresh_ element
    static double erfinv_ws_thresh() {
        return boost::math::erfc_inv(well_seperated_thresh_);
    }

    /*! \brief Four loop Inner Screen.
     *
     * Only use QQR once all 4 shells are present. This means
     * for now only Schwarz is used for loop skipping.
     */
    bool skip(int64_t, int64_t, int64_t, int64_t) override;
};

struct init_qqr_screen {
    template <typename E>
    QQR operator()(madness::World &world, ShrPool<E> &engs,
                   Basis const &bs, double threshold = 1e-10) {
        auto schwarz_screen = init_schwarz_screen{threshold}(world, engs, bs);
        return QQR(std::move(schwarz_screen), bs.flattened_shells());
    }
};

}  // namespace  gaussian
}  // namespace  lcao
}  // namespace  mpqc


#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_SCREENING_QQR_SCREENING_H_
