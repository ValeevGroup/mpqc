#pragma once
#ifndef MPQC_INTEGRALS_QRRSCREENING
#define MPQC_INTEGRALS_QRRSCREENING

#include "../common/typedefs.h"
#include "task_integrals_common.h"
#include "integral_screeners.h"

#include "schwarz_screen.h"

namespace mpqc {
namespace integrals {
namespace detail {

/*! \todo move this to c++ file. Especially the shell_weighted_center and 
*   pair_extent function
*/
class CachedShellInfo {
  private:
    MatrixD pair_extents_;

  public:
    CachedShellInfo(ShellVec const &shells0, ShellVec const &shells1)
            : pair_extents_(MatrixD(shells0.size(), shells1.size())) {

        for (auto s0 = 0ul; s0 < shells0.size(); ++s0) {
            auto const &sh0 = shells0[s0];

            for (auto s1 = 0ul; s1 < shells1.size(); ++s1) {
                auto const &sh1 = shells1[s1];

                auto r_01 = shell_weighted_center(sh0, sh1);
                pair_extents_(s0, s1) = pair_extent(sh0, sh1, r_01);
            }
        }
    }

  private:
    double pair_extent(Shell const &sh0, Shell const &sh1, Vec3D const &r_01) {

        double max_extent = 0.0;

        auto const &exp0 = sh0.alpha;
        auto const &exp1 = sh1.alpha;

        auto const &O0 = sh0.O;
        Vec3D center0 = {O0[0], O0[1], O0[2]};

        auto const &O1 = sh1.O;
        Vec3D center1 = {O1[0], O1[1], O1[2]};

        const auto ncontr0 = sh0.ncontr();
        const auto ncontr1 = sh1.ncontr();

        for (auto i = 0ul; i < ncontr0; ++i) {
            const auto i_exp = exp0[i];
            const auto i_scaled_center = i_exp * center0;

            for (auto j = 0ul; j < ncontr1; ++j) {
                const auto j_exp = exp1[j];

                const auto exp_sum = i_exp + j_exp;
                const auto inv_sum = 1 / exp_sum;

                Vec3D r_ij = inv_sum * (i_scaled_center + j_exp * center1);
                const auto diff = (r_ij - r_01).norm();
                // This magic number is my inverse erf for now.
                const auto ext_ij = std::sqrt(2 * inv_sum) * 1.163 + diff;

                max_extent = std::max(ext_ij, max_extent);
            }
        }

        return max_extent;
    }

    Vec3D shell_weighted_center(Shell const &sh0, Shell const &sh1) {

        Vec3D center = {0.0, 0.0, 0.0};
        double sum_of_coeff_products = 0.0;

        auto const &exp0 = sh0.alpha;
        auto const &exp1 = sh1.alpha;

        auto const &O0 = sh0.O;
        Vec3D center0 = {O0[0], O0[1], O0[2]};

        auto const &O1 = sh1.O;
        Vec3D center1 = {O1[0], O1[1], O1[2]};

        auto const &contr0 = sh0.contr;
        auto const &contr1 = sh1.contr;

        const auto ncontr0 = sh0.ncontr();
        const auto ncontr1 = sh1.ncontr();

        for (auto i = 0ul; i < ncontr0; ++i) {
            const auto i_coeff = contr0[i].coeff[0];
            const auto i_exp = exp0[i];
            const auto i_scaled_center = i_exp * center0;

            for (auto j = 0ul; j < ncontr1; ++j) {
                const auto j_coeff = contr1[j].coeff[0];
                const auto j_exp = exp1[j];

                const auto product = std::abs(i_coeff * j_coeff);
                const auto scale = product / (i_exp + j_exp);

                center += scale * (i_scaled_center + j_exp * center1);

                sum_of_coeff_products += product;
            }
        }

        center *= 1 / sum_of_coeff_products;
        return center;
    }
};

} // namespace detail

/*! \brief Class to implement QRR Screening.
 *
 * Based on: Maurer, S. A.; Lambrecht, D. S.; Flaig, D. JCP. 136, 144107 (2012).
 *
 * Reverse Engineered from David Hollman's code in MPQC2/3
 *
 */
class QRR : public SchwarzScreen {
    private:
        CachedShellInfo shell_info_ab;
        CachedShellInfo shell_info_cd;

    public:
};

} // namespace integrals
} // namespace mpqc


#endif // MPQC_INTEGRALS_QRRSCREENING
