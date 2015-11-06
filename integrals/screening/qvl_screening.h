#pragma once
#ifndef MPQC_INTEGRALS_QVLSCREENING
#define MPQC_INTEGRALS_QVLSCREENING

#include "../../common/typedefs.h"
#include "../task_integrals_common.h"

#include "screen_base.h"
#include "schwarz_screen.h"
#include "qvl_shell_info.h"

#include <boost/math/special_functions/erf.hpp>

namespace mpqc {
namespace integrals {

/*! \brief Class to implement QVL Screening.
 *
 * Reverse Engineered from David Hollman's code in MPQC2/3
 */
class QVl : public SchwarzScreen {
  private:
    detail::QVlShellInfo shell_info_;
    static double well_seperated_thresh_;

  public:
    QVl() = default;
    QVl(SchwarzScreen ss, ShellVec const &a_shells, ShellVec const &cd_shells)
            : SchwarzScreen(std::move(ss)),
              shell_info_(a_shells, cd_shells, erfinv_ws_thresh()) {}

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
     * Only use QVl once all 4 shells are present. This means
     * for now only Schwarz is used for loop skipping.
     */
    bool skip(int64_t a, int64_t c, int64_t d) override {
        auto local_Qcd = Qcd(c, d);
        auto local_Qa = Qab(a);
        auto QaQcd = local_Qcd * local_Qa;

        auto ash = Qab_f2s(a);
        auto csh = Qcd_f2s(c);
        auto dsh = Qcd_f2s(d);

        if (QaQcd >= skip_threshold()) {
            double Rnorm = (shell_info_.a_center(ash)
                            - shell_info_.cd_center(csh, dsh)).norm();

            if (Rnorm > shell_info_.a_extent(ash)
                        + shell_info_.cd_extent(csh, dsh)) {
                QaQcd = local_Qcd
                        * shell_info_.scaling_factor(ash, csh, dsh, Rnorm);
            } 
        }

        return (QaQcd < skip_threshold()) ? true : false;
    }
};

struct init_qvl_screen {
    template <typename E>
    QVl operator()(madness::World &world, ShrPool<E> &engs,
                   basis::Basis const &dfbs, basis::Basis const &bs,
                   double threshold = 1e-10) {
        auto schwarz_screen
              = init_schwarz_screen{threshold}(world, engs, dfbs, bs);
        return QVl(std::move(schwarz_screen), dfbs.flattened_shells(),
                   bs.flattened_shells());
    }
};

} // namespace integrals
} // namespace mpqc


#endif // MPQC_INTEGRALS_QVLSCREENING
