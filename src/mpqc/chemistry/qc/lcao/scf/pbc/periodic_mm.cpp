#include "mpqc/chemistry/qc/lcao/scf/pbc/periodic_mm.h"

namespace mpqc {
namespace pbc {
namespace detail {

BasisPairInfo::BasisPairInfo(std::shared_ptr<const Basis> bs0,
                             std::shared_ptr<const Basis> bs1,
                             const double thresh)
    : bs0_(bs0),
      bs1_(bs1),
      thresh_(thresh),
      nshells0_(bs0->flattened_shells().size()),
      nshells1_(bs1->flattened_shells().size()) {
  pair_extents_.resize(nshells0_, nshells1_);
  pair_centers_ = std::vector<std::vector<Vector3d>>(
      nshells0_, std::vector<Vector3d>(nshells1_));

  const auto &shellvec0 = bs0->flattened_shells();
  const auto &shellvec1 = bs1->flattened_shells();
  for (auto s0 = 0ul; s0 != nshells0_; ++s0) {
    const auto &shell0 = shellvec0[s0];

    for (auto s1 = 0ul; s1 != nshells1_; ++s1) {
      const auto &shell1 = shellvec1[s1];

      auto &extent = pair_extents_(s0, s1);
      auto &center = pair_centers_[s0][s1];
      std::make_pair(center, extent) = shell_pair_center_extent(shell0, shell1);
    }
  }
}

std::pair<Vector3d, double> BasisPairInfo::shell_pair_center_extent(
    const Shell &sh0, const Shell &sh1) {
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
  auto update = [&center, &extent, &l, &r, &n_a2b, &is_between](
                    const Vector3d &center_new, const double extent_new) {
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

double BasisPairInfo::prim_pair_extent(const double exp0, const double exp1,
                                       const double rab) {
  const auto exp_sum = exp0 + exp1;
  const auto exp_mult = exp0 * exp1;

  // compute overlap integral between two primitives, neglecting angular parts
  const auto prefactor = std::pow(4.0 * exp_mult / (exp_sum * exp_sum), 0.75);
  const auto exponent = std::exp(-1.0 * (exp_mult / exp_sum) * rab * rab);
  const auto S = prefactor * exponent;

  const auto extent2 =
      (-1.0 * std::log(thresh_) + std::log(S) + 0.5 * std::log(exp_sum)) /
      exp_sum;
  return std::sqrt(extent2);
}

double BasisPairInfo::extent(int64_t ij) const {
  const auto i = ij / nshells1_;
  const auto j = ij % nshells1_;
  return extent(i, j);
}

Vector3d BasisPairInfo::center(int64_t ij) const {
  const auto i = ij / nshells1_;
  const auto j = ij % nshells1_;
  return center(i, j);
}

double BasisPairInfo::max_distance_to(const Vector3d &ref_point) {
  double distance = 0.0;
  for (auto p = 0ul; p != nshells0_ * nshells1_; ++p) {
    const auto r = (center(p) - ref_point).norm();
    distance = std::max(distance, r);
  }
  return distance;
}

}  // namespace detail
}  // namespace pbc
}  // namespace mpqc
