#include "qqr_screening.h"

namespace mpqc {
namespace integrals {

bool QQR::skip(int64_t ordA, int64_t ordB, int64_t ordC, int64_t ordD,
               Shell const &, Shell const &, Shell const &,
               Shell const &) {

    auto QabQcd = four_center_estimate(ordA, ordB, ordC, ordD);

    // Don't bother with QQR if QQ is good enought
    if (QabQcd > skip_threshold()) {
        auto const &R_ab = shell_info_ab_.center(ordA, ordB);
        auto const &R_cd = shell_info_cd_.center(ordC, ordD);
        const auto R = (R_ab - R_cd).norm();
        if (R > 1) {
            const auto ext_ab = shell_info_ab_.extent(ordA, ordB);
            const auto ext_cd = shell_info_cd_.extent(ordC, ordD);

            const auto ext_total = ext_cd + ext_ab;
            const auto r_minus_extents = R - ext_total;

            if (r_minus_extents > 1) {
                QabQcd /= r_minus_extents;
            }
        }
    }

    return (QabQcd < skip_threshold()) ? true : false;
}


// Have to pick some value, going with tight threshold for now.
double QQR::well_seperated_thresh_ = 0.2;


} // namespace integrals
} // namespace mpqc
