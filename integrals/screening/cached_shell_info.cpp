#include "cached_shell_info.h"
#include "../../include/libint.h"

namespace mpqc {
namespace integrals {
namespace detail {

CachedShellInfo::CachedShellInfo(ShellVec const &shells0,
                                 ShellVec const &shells1, double thresh)
        : pair_extents_(MatrixD(shells0.size(), shells1.size())),
          pair_centers_(shells0.size(), std::vector<Vec3D>(shells1.size())),
          erfinv_thr_(thresh) {
    for (auto s0 = 0ul; s0 < shells0.size(); ++s0) {
        auto const &sh0 = shells0[s0];

        for (auto s1 = 0ul; s1 < shells1.size(); ++s1) {
            auto const &sh1 = shells1[s1];

            const auto r_01 = shell_weighted_center(sh0, sh1);
            pair_extents_(s0, s1) = pair_extent(sh0, sh1, r_01);
            pair_centers_[s0][s1] = r_01;
        }
    }
}

double CachedShellInfo::pair_extent(Shell const &sh0, Shell const &sh1,
                                    Vec3D const &r_01) {

    double max_extent = 0.0;

    auto const &exp0 = sh0.alpha;
    auto const &exp1 = sh1.alpha;

    auto const &O0 = sh0.O;
    Vec3D center0 = {O0[0], O0[1], O0[2]};

    auto const &O1 = sh1.O;
    Vec3D center1 = {O1[0], O1[1], O1[2]};

    const auto nprim0 = sh0.nprim();
    const auto nprim1 = sh1.nprim();

    for (auto i = 0ul; i < nprim0; ++i) {
        const auto i_exp = exp0[i];
        const auto i_scaled_center = i_exp * center0;

        for (auto j = 0ul; j < nprim1; ++j) {
            const auto j_exp = exp1[j];

            const auto exp_sum = i_exp + j_exp;
            const auto inv_sum = 1 / exp_sum;

            Vec3D r_ij = inv_sum * (i_scaled_center + j_exp * center1);
            const auto diff = (r_ij - r_01).norm();
            const auto ext_ij = std::sqrt(2 * inv_sum) * erfinv_thr_ + diff;

            max_extent = std::max(ext_ij, max_extent);
        }
    }

    return max_extent;
}

Vec3D CachedShellInfo::shell_weighted_center(Shell const &sh0,
                                             Shell const &sh1) {

    Vec3D center = {0.0, 0.0, 0.0};
    double sum_of_coeff_products = 0.0;

    auto const &exp0 = sh0.alpha;
    auto const &exp1 = sh1.alpha;

    auto const &O0 = sh0.O;
    Vec3D center0 = {O0[0], O0[1], O0[2]};

    auto const &O1 = sh1.O;
    Vec3D center1 = {O1[0], O1[1], O1[2]};

    const auto nprim0 = sh0.nprim();
    const auto nprim1 = sh1.nprim();

    for (auto i = 0ul; i < nprim0; ++i) {
        const auto i_coeff = sh0.contr[0].coeff[i];
        const auto i_exp = exp0[i];
        const auto i_scaled_center = i_exp * center0;

        for (auto j = 0ul; j < nprim1; ++j) {
            const auto j_coeff = sh1.contr[0].coeff[j];
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

} // namespace detail
} // namespace integrals
} // namespace mpqc
