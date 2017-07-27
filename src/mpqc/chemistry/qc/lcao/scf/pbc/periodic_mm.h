#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_MM_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_MM_H_

#include "mpqc/chemistry/qc/lcao/factory/periodic_ao_factory.h"
#include "mpqc/chemistry/qc/lcao/integrals/screening/cached_shell_info.h"

namespace mpqc {
namespace pbc {
namespace mm {

/*!
 * \brief This class compute centers and extents for periodic multiple
 * approximation.
 */
template <typename Factory>
class PeriodicMM {
 public:
  using Shell = ::mpqc::lcao::gaussian::Shell;
  using ShellVec = ::mpqc::lcao::gaussian::ShellVec;
  using Basis = ::mpqc::lcao::gaussian::Basis;
  using ShellInfo = ::mpqc::lcao::gaussian::detail::CachedShellInfo;

  PeriodicMM(Factory &ao_factory_, double mm_thresh = 1.0e-8, double ws = 3)
      : ao_factory_(ao_factory_), mm_thresh_(mm_thresh), ws_(ws) {}

 private:
  Factory &ao_factory_;
  const double mm_thresh_;  /// threshold of multiple approximation error
  const double ws_;         /// well-separateness criterion

  std::shared_ptr<const Basis> obs_;
  std::shared_ptr<const Basis> dfbs_;
  RowMatrixXd pair_extents_;
  std::vector<std::vector<Vector3d>> pair_centers_;

 public:
  Vector3d uc_pair_center(const Vector3d &L) {}

 private:
  /*!
   * \brief basis_pair_center
   * \param L
   * \param bs0
   * \param bs1
   */
  void basis_pair_center(const Basis &bs0, const Basis &bs1) {
    const auto &shellvec0 = bs0.flattened_shells();
    const auto &shellvec1 = bs1.flattened_shells();

    const auto nshells0 = shellvec0.size();
    const auto nshells1 = shellvec1.size();

    pair_extents_.resize(nshells0, nshells1);
    pair_centers_ = std::vector<std::vector<Vector3d>>(
        nshells0, std::vector<Vector3d>(nshells1));

    for (auto s0 = 0ul; s0 != nshells0; ++s0) {
      const auto &shell0 = shellvec0[s0];

      for (auto s1 = 0ul; s1 != nshells1; ++s1) {
        const auto &shell1 = shellvec1[s1];

        auto &extent = pair_extents_(s0, s1);
        auto &center = pair_centers_[s0][s1];
        std::make_pair(center, extent) =
            shell_pair_center_extent(shell0, shell1);
      }
    }
  }

  /*!
   * \brief This computes the center and the extent of the product of two shells
   * \param sh0
   * \param sh1
   * \return
   */
  std::pair<Vector3d, double> shell_pair_center_extent(const Shell &sh0,
                                                       const Shell &sh1) {
    const auto &O0 = sh0.O;
    const auto &O1 = sh1.O;
    Vector3d center0 = {O0[0], O0[1], O0[2]};
    Vector3d center1 = {O1[0], O1[1], O1[2]};
    const auto &exp0 = sh0.alpha;
    const auto &exp1 = sh1.alpha;

    Vector3d R_a2b = center1 - center0;
    // compute |O1 - O0|
    double rab = R_a2b.norm();
    // unit vector from O0 to O1
    Vector3d n_a2b = 1.0 / rab * R_a2b;

    Vector3d center = 0.5 * (center0 + center1);
    double extent = 0.5 * rab;
    Vector3d l = center0;
    Vector3d r = center1;

    // determine if \c a is between \c b and \c c (a, b, and c are collinear)
    auto is_between = [](const Vector3d &a, const Vector3d &b,
                         const Vector3d &c) {
      const Vector3d bc = c - b;

      bool is_x_between = bc[0] > 0 ? (a[0] >= b[0] && a[0] <= c[0])
                                    : (a[0] <= b[0] && a[0] >= c[0]);
      bool is_y_between = bc[1] > 0 ? (a[1] >= b[1] && a[1] <= c[1])
                                    : (a[1] <= b[1] && a[1] >= c[1]);
      bool is_z_between = bc[2] > 0 ? (a[2] >= b[2] && a[2] <= c[2])
                                    : (a[2] <= b[2] && a[2] >= c[2]);

      return (is_x_between && is_y_between && is_z_between);
    };

    //  update center and extent for shell pair
    auto update = [&center, &extent, &l, &r, &n_a2b](const Vector3d &center_new,
                                                     const double extent_new) {
      Vector3d l_new = center_new - extent_new * n_a2b;
      Vector3d r_new = center_new + extent_new * n_a2b;
      l = is_between(l_new, center, l) ? l : l_new;
      r = is_between(r_new, center, r) ? r : r_new;
      center = 0.5 * (l + r);
      extent = 0.5 * ((r - l).norm());
    };

    const auto nprim0 = sh0.nprim();
    const auto nprim1 = sh1.nprim();
    for (auto i = 0ul; i != nprim0; ++i) {
      const auto exp_i = exp0[i];

      for (auto j = 0ul; j != nprim1; ++j) {
        const auto exp_j = exp1[j];

        const Vector3d center_ij =
            (exp_i * center0 + exp_j * center1) / (exp_i + exp_j);

        const double extent_ij = prim_pair_extent(exp_i, exp_j, rab);
        update(center_ij, extent_ij);
      }
    }

    return std::make_pair(center, extent);
  }

  /*!
   * \brief This computes the center of the product of two primitives
   * \param exp0
   * \param exp1
   * \param thresh
   * \return
   */
  double prim_pair_extent(const double exp0, const double exp1,
                          const double rab) {
    const auto exp_sum = exp0 + exp1;
    const auto exp_mult = exp0 * exp1;

    // compute overlap integral between two primitives, neglecting angular parts
    const auto prefactor = std::pow(4.0 * exp_mult / (exp_sum * exp_sum), 0.75);
    const auto exponent = std::exp(-1.0 * (exp_mult / exp_sum) * rab * rab);
    const auto S = prefactor * exponent;

    const auto extent2 =
        (-1.0 * std::log(mm_thresh_) + std::Log(S) + 0.5 * std::log(exp_sum)) /
        exp_sum;
    return std::sqrt(extent2);
  }
};

}  // namespace mm
}  // namespace pbc
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_MM_H_
