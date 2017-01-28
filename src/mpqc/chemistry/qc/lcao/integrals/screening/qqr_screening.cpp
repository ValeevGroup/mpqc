#include "mpqc/chemistry/qc/lcao/integrals/screening/qqr_screening.h"

namespace mpqc {
namespace lcao {
namespace gaussian {

bool QQR::skip(int64_t a, int64_t b, int64_t c, int64_t d) {


    auto QabQcd = four_center_est(a,b,c,d);

    // Get the shell indices of the function indices
    auto a_sh = Qab_f2s(a);
    auto b_sh = Qab_f2s(b);
    auto c_sh = Qab_f2s(c);
    auto d_sh = Qab_f2s(d);

    // Don't bother with QQR if QQ is good enought
    if (QabQcd > skip_threshold()) {
        auto const &R_ab = shell_info_ab_.center(a_sh, b_sh);
        auto const &R_cd = shell_info_ab_.center(c_sh, d_sh);
        const auto R = (R_ab - R_cd).norm();

        if (R > 1) {
            const auto ext_ab = shell_info_ab_.extent(a_sh, b_sh);
            const auto ext_cd = shell_info_ab_.extent(c_sh, d_sh);

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


}  // namespace  gaussian
}  // namespace  lcao
}  // namespace  mpqc
