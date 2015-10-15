#pragma once
#ifndef MPQC_INTEGRALS_QQRSCREENING
#define MPQC_INTEGRALS_QQRSCREENING

#include "../../common/typedefs.h"
#include "../task_integrals_common.h"

#include "screen_base.h"
#include "schwarz_screen.h"
#include "cached_shell_info.h"

#include <boost/math/special_functions/erf.hpp>

namespace mpqc {
namespace integrals {

/*! \brief Class to implement QQR Screening.
 *
 * Based on: Maurer, S. A.; Lambrecht, D. S.; Flaig, D. JCP. 136, 144107 (2012).
 *
 * Reverse Engineered from David Hollman's code in MPQC2/3
 *
 */
class QQR : public SchwarzScreen {
  private:
    detail::CachedShellInfo shell_info_ab_;
    detail::CachedShellInfo shell_info_cd_;
    static double well_seperated_thresh_;

  public:
    QQR() = default;
    QQR(SchwarzScreen ss, ShellVec const &sh_a, ShellVec const &sh_b,
        ShellVec const &sh_c, ShellVec const &sh_d)
            : SchwarzScreen(std::move(ss)),
              shell_info_ab_(sh_a, sh_b, erfinv_ws_thresh()),
              shell_info_cd_(sh_c, sh_d, erfinv_ws_thresh()) {}

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
     * For now only use QQR once all 4 shells are present. This means
     * for now only Schwarz is used for loop skipping.
     */
    bool
    skip(int64_t ordA, int64_t ordB, int64_t ordC, int64_t ordD, Shell const &,
         Shell const &, Shell const &, Shell const &) override;
};

struct init_qqr_screen {

    template <typename E>
    QQR operator()(detail::IdxVec const &idx, detail::ShrBases<4> const &bases,
                   ShrPool<E> const &engs) {
        auto schwarz_screen = init_schwarz_screen{}(idx, bases, engs);

        auto const &bases_ref = *bases;

        auto const &sh_a = bases_ref[0].cluster_shells()[idx[0]];
        auto const &sh_b = bases_ref[1].cluster_shells()[idx[1]];
        auto const &sh_c = bases_ref[2].cluster_shells()[idx[2]];
        auto const &sh_d = bases_ref[3].cluster_shells()[idx[3]];

        return QQR(std::move(schwarz_screen), sh_a, sh_b, sh_c, sh_d);
    }
};

} // namespace integrals
} // namespace mpqc


#endif // MPQC_INTEGRALS_QQRSCREENING
