#include "mpqc/chemistry/qc/integrals/screening/qvl_shell_info.h"

namespace mpqc {
namespace lcao {
namespace gaussian {
namespace detail {

namespace {
constexpr int64_t dfac(int64_t x) {
    return (x == -1 || x == 0) ? 1 : x * dfac(x - 2);
}

constexpr double pi = M_PI;

constexpr double abs(double x) { return (x < 0.0) ? -x : x; }

constexpr double tol = 0.000000001;

constexpr double sqrt_help(double x, double g) {
    return abs(g - x / g) < tol ? g : sqrt_help(x, (g + x / g) / 2.0);
}

constexpr double sqrt(double x) { return sqrt_help(x, 1.0); }

constexpr double QVl_top_term(int64_t l) {
    return pi * sqrt(2 * dfac(2 * l - 1));
}


constexpr std::array<double, 9> sqrt_2_dfac_term
      = {{QVl_top_term(0), QVl_top_term(1), QVl_top_term(2),
         QVl_top_term(3), QVl_top_term(4), QVl_top_term(5),
         QVl_top_term(6), QVl_top_term(7), QVl_top_term(8)}};

constexpr double power_term(int64_t l) { return double(2 * l + 3) / 4.0; }

constexpr std::array<double, 9> exp_power_term
      = {{power_term(0), power_term(1), power_term(2),
         power_term(3), power_term(4), power_term(5),
         power_term(6), power_term(7), power_term(8)}};
}

QVlShellInfo::QVlShellInfo(ShellVec const &a_shells, ShellVec const &cd_shells,
                           double thresh)
        : a_extents_(a_shells.size(), 0.0),
          cd_extents_(cd_shells.size(), cd_shells.size()),
          a_centers_(a_shells.size()),
          cd_centers_(cd_shells.size(), std::vector<Vector3d>(cd_shells.size())),
          a_factor_(a_shells.size(), 0.0),
          a_am_(a_shells.size(), 0.0),
          one_over_gamma_ab_pow_(cd_shells.size(), cd_shells.size()),
          erfinv_thr_(thresh) {

    for (auto s = 0ul; s < a_shells.size(); ++s) {
        auto const &sh = a_shells[s];
        const auto am = sh.contr[0].l;

        a_am_[s] = am;

        a_centers_[s] = {sh.O[0], sh.O[1], sh.O[2]};

        auto min_exp = a_shell_min_exp(sh);
        a_extents_[s] = sqrt(2.0 / min_exp) * erfinv_thr_;

        a_factor_[s] = sqrt_2_dfac_term[am]
                       / std::pow(min_exp, exp_power_term[am]);
    }

    for (auto s0 = 0ul; s0 < cd_shells.size(); ++s0) {
        auto const &sh0 = cd_shells[s0];

        auto min_exp0 = a_shell_min_exp(sh0);

        for (auto s1 = 0ul; s1 < cd_shells.size(); ++s1) {
            auto const &sh1 = cd_shells[s1];

            auto min_exp1 = a_shell_min_exp(sh1);

            const auto gamma14 = std::pow(min_exp0 + min_exp1, 0.25);
            one_over_gamma_ab_pow_(s0, s1) = 1/gamma14;

            const auto r_cd = shell_weighted_center(sh0, sh1);
            cd_extents_(s0, s1) = pair_extent(sh0, sh1, r_cd);
            cd_centers_[s0][s1] = std::move(r_cd);
        }
    }
}

double QVlShellInfo::a_shell_min_exp(Shell const &sh) {
    auto const &exp = sh.alpha;
    return *std::min_element(std::begin(exp), std::end(exp));
}


Vector3d
QVlShellInfo::shell_weighted_center(Shell const &sh0, Shell const &sh1) const {

    Vector3d center = {0.0, 0.0, 0.0};
    double sum_of_coeff_products = 0.0;

    auto const &exp0 = sh0.alpha;
    auto const &exp1 = sh1.alpha;

    auto const &O0 = sh0.O;
    Vector3d center0 = {O0[0], O0[1], O0[2]};

    auto const &O1 = sh1.O;
    Vector3d center1 = {O1[0], O1[1], O1[2]};

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

double QVlShellInfo::pair_extent(Shell const &sh0, Shell const &sh1,
                                 Vector3d const &r_01) {

    double max_extent = 0.0;

    auto const &exp0 = sh0.alpha;
    auto const &exp1 = sh1.alpha;

    auto const &O0 = sh0.O;
    Vector3d center0 = {O0[0], O0[1], O0[2]};

    auto const &O1 = sh1.O;
    Vector3d center1 = {O1[0], O1[1], O1[2]};

    const auto nprim0 = sh0.nprim();
    const auto nprim1 = sh1.nprim();

    for (auto i = 0ul; i < nprim0; ++i) {
        const auto i_exp = exp0[i];
        const auto i_scaled_center = i_exp * center0;

        for (auto j = 0ul; j < nprim1; ++j) {
            const auto j_exp = exp1[j];

            const auto exp_sum = i_exp + j_exp;
            const auto inv_sum = 1 / exp_sum;

            Vector3d r_ij = inv_sum * (i_scaled_center + j_exp * center1);
            const auto diff = (r_ij - r_01).norm();
            const auto ext_ij = std::sqrt(2 * inv_sum) * erfinv_thr_ + diff;

            max_extent = std::max(ext_ij, max_extent);
        }
    }

    return max_extent;
}

}  // namespace  detail
}  // namespace  gaussian
}  // namespace  lcao
}  // namespace  mpqc
