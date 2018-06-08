#include "mpqc/chemistry/qc/lcao/scf/pbc/periodic_ma.h"

namespace mpqc {
namespace pbc {
namespace detail {

BasisPairInfo::BasisPairInfo(std::shared_ptr<const Basis> bs0,
                             std::shared_ptr<const Basis> bs1,
                             const double thresh, const double small_extent)
    : thresh_(thresh), small_extent_(small_extent), bs0_(bs0), bs1_(bs1) {
  const auto &shellvec0 = bs0->flattened_shells();
  const auto &shellvec1 = bs1->flattened_shells();

  nshells0_ = shellvec0.size();
  nshells1_ = shellvec1.size();
  npairs_ = nshells0_ * nshells1_;
  pair_extents_.resize(nshells0_, nshells1_);
  pair_centers_ = std::vector<std::vector<Vector3d>>(
      nshells0_, std::vector<Vector3d>(nshells1_));

  for (auto s0 = 0ul; s0 != nshells0_; ++s0) {
    const auto &shell0 = shellvec0[s0];
    for (auto s1 = 0ul; s1 != nshells1_; ++s1) {
      const auto &shell1 = shellvec1[s1];
      const auto pair = shell_pair_center_extent(shell0, shell1);
      pair_centers_[s0][s1] = pair.first;
      pair_extents_(s0, s1) = pair.second;
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

  // compute |O1 - O0|
  double rab = (center1 - center0).norm();

  //  update extent for shell pair w.r.t. a reference center
  auto update_extent_wrt_center =
      [](const Vector3d &center_ref, const Vector3d &center_ij,
         const double extent_ij, const double extent_old) {
        const auto extent_new = extent_ij + (center_ij - center_ref).norm();
        return std::max(extent_old, extent_new);
      };

  const auto nprim0 = sh0.nprim();
  const auto nprim1 = sh1.nprim();
  double extent_wrt_c0 = 0.0;
  double extent_wrt_c1 = 0.0;
  for (auto i = 0ul; i != nprim0; ++i) {
    const auto exp_i = exp0[i];

    for (auto j = 0ul; j != nprim1; ++j) {
      const auto exp_j = exp1[j];

      const Vector3d center_ij =
          (exp_i * center0 + exp_j * center1) / (exp_i + exp_j);
      const double extent_ij = prim_pair_extent(exp_i, exp_j, rab);
      extent_wrt_c0 = update_extent_wrt_center(center0, center_ij, extent_ij,
                                               extent_wrt_c0);
      extent_wrt_c1 = update_extent_wrt_center(center1, center_ij, extent_ij,
                                               extent_wrt_c1);
    }
  }

  return (extent_wrt_c0 <= extent_wrt_c1)
             ? std::make_pair(center0, extent_wrt_c0)
             : std::make_pair(center1, extent_wrt_c1);
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
      (std::log(S) - std::log(thresh_) + 0.5 * std::log(exp_sum)) / exp_sum;
  const auto extent = (extent2 < 0) ? small_extent_ : std::sqrt(extent2);

  return extent;
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
