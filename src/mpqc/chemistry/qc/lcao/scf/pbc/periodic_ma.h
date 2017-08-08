#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_MM_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_MM_H_

#include "mpqc/chemistry/qc/lcao/factory/periodic_ao_factory.h"
#include "mpqc/chemistry/qc/lcao/scf/pbc/util.h"

namespace mpqc {
namespace pbc {

namespace detail {

/*!
 * \brief Class to hold information of basis pairs
 */
class BasisPairInfo {
 public:
  using Shell = ::mpqc::lcao::gaussian::Shell;
  using Basis = ::mpqc::lcao::gaussian::Basis;

  BasisPairInfo() = default;
  BasisPairInfo(std::shared_ptr<const Basis> bs0,
                std::shared_ptr<const Basis> bs1, const double thresh = 1.0e-6,
                const double small_extent = 0.01);

 private:
  const double thresh_;        ///> threshold of multiple approximation error
  const double small_extent_;  ///> a small extent will be used in case no real
                               /// solution exists

  std::shared_ptr<const Basis> bs0_;
  std::shared_ptr<const Basis> bs1_;
  size_t nshells0_;
  size_t nshells1_;
  size_t npairs_;
  RowMatrixXd pair_extents_;
  std::vector<std::vector<Vector3d>> pair_centers_;

 public:
  double extent(int64_t i, int64_t j) const { return pair_extents_(i, j); }
  Vector3d center(int64_t i, int64_t j) const { return pair_centers_[i][j]; }

  double extent(int64_t ij) const;
  Vector3d center(int64_t ij) const;

  size_t npairs() const { return npairs_; }

  /*!
   * \brief This computes maximum distance between \cref_point and all charge
   * centers comprising the product density
   */
  double max_distance_to(const Vector3d &ref_point);

 private:
  /*!
   * \brief This computes the center and the extent of the product of two shells
   * \param sh0
   * \param sh1
   * \return
   */
  std::pair<Vector3d, double> shell_pair_center_extent(const Shell &sh0,
                                                       const Shell &sh1);

  /*!
   * \brief This computes the center of the product of two primitives
   * \param exp0
   * \param exp1
   * \param thresh
   * \return
   */
  double prim_pair_extent(const double exp0, const double exp1,
                          const double rab);
};

}  // namespace detail

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

  PeriodicMM(Factory &ao_factory, double mm_thresh = 1.0e-8, double ws = 3.0)
      : ao_factory_(ao_factory), mm_thresh_(mm_thresh), ws_(ws) {
    auto &world = ao_factory_.world();
    obs_ = ao_factory_.basis_registry()->retrieve(OrbitalIndex(L"λ"));
    dfbs_ = ao_factory_.basis_registry()->retrieve(OrbitalIndex(L"Κ"));

    dcell_ = ao_factory_.unitcell().dcell();
    R_max_ = ao_factory_.R_max();
    RJ_max_ = ao_factory_.RJ_max();
    RD_max_ = ao_factory_.RD_max();
    R_size_ = ao_factory_.R_size();
    RJ_size_ = ao_factory_.RJ_size();
    RD_size_ = ao_factory_.RD_size();

    using ::mpqc::lcao::gaussian::detail::shift_basis_origin;
    using ::mpqc::lcao::detail::direct_3D_idx;
    Vector3d zero_shift_base(0.0, 0.0, 0.0);
    // determine non-negligible shell pairs from ref unit cell and all nearby
    // unit cells
    auto basisR = shift_basis_origin(*obs_, zero_shift_base, R_max_, dcell_);
    auto significant_shellpair_list =
        detail::parallel_compute_shellpair_list(world, *obs_, *basisR);
    // determine necessary nearby unit cells involved in product density
    for (auto R = 0; R != R_size_; ++R) {
      const auto R_3D = direct_3D_idx(R, R_max_);

      const auto nshells = obs_->flattened_shells().size();
      const auto shell1_min = nshells * R;
      const auto shell1_max = shell1_min + nshells;

      auto is_significant = false;
      for (auto shell0 = 0; shell0 != nshells; ++shell0) {
        for (const auto &shell1 : significant_shellpair_list[shell0]) {
          if (shell1 >= shell1_min && shell1 < shell1_max) {
            is_significant = true;
            uc_near_list_.emplace_back(R_3D);
            break;
          }
        }
        if (is_significant) break;
      }
    }

    // compute centers and extents of product density between the referene
    // unit cell and its neighbours
    ref_pairs_ = compute_basis_pairs();
    // compute maximum distance between the center of reference unit cell and
    // all charge centers
    const auto ref_center = 0.5 * dcell_;
    max_distance_to_refcenter_ = ref_pairs_->max_distance_to(ref_center);
  }

 private:
  Factory &ao_factory_;
  const double mm_thresh_;  /// threshold of multiple approximation error
  const double ws_;         /// well-separateness criterion

  std::shared_ptr<Basis> obs_;
  std::shared_ptr<Basis> dfbs_;
  Vector3d dcell_;
  Vector3i R_max_;
  Vector3i RJ_max_;
  Vector3i RD_max_;
  int64_t R_size_;
  int64_t RJ_size_;
  int64_t RD_size_;

  std::vector<Vector3i> uc_near_list_;
  std::shared_ptr<detail::BasisPairInfo> ref_pairs_;
  double max_distance_to_refcenter_;

 public:
  /*!
   * \brief This determines if a unit cell \c uc_ket is in the crystal far
   * field of the bra unit cell \c uc_bra.
   */
  bool is_uc_in_CFF(const Vector3i &uc_ket,
                    const Vector3i &uc_bra = {0, 0, 0}) {
    using ::mpqc::lcao::gaussian::detail::shift_basis_origin;

    Vector3d vec_bra = uc_bra.cast<double>().cwiseProduct(dcell_);
    Vector3d vec_ket = uc_ket.cast<double>().cwiseProduct(dcell_);
    const auto vec_rel = vec_ket - vec_bra;

    const auto npairs = ref_pairs_->npairs();
    // CFF condition #1: all charge distributions are well seperated
    auto condition1 = true;
    for (auto p0 = 0ul; p0 != npairs; ++p0) {
      const auto center0 = ref_pairs_->center(p0) + vec_bra;
      const auto extent0 = ref_pairs_->extent(p0);
      for (auto p1 = 0ul; p1 != npairs; ++p1) {
        const auto center1 = ref_pairs_->center(p1) + vec_ket;
        const auto extent1 = ref_pairs_->extent(p1);

        const double rab = (center1 - center0).norm();
        if (rab < (extent0 + extent1)) {
          condition1 = false;
          break;
        }
      }
      if (!condition1) break;
    }

    // CFF condition #2: |L| >= ws * (r0_max + r1_max)
    auto condition2 = false;
    const auto L = vec_rel.norm();
    const auto uc_center_bra = vec_bra + 0.5 * dcell_;
    condition2 = (L >= ws_ * (max_distance_to_refcenter_ * 2.0)) ? true : false;

    // test:
    {
      ExEnv::out0() << "Vector bra corner = " << vec_bra.transpose()
                    << std::endl;
      ExEnv::out0() << "Vector ket corner = " << vec_ket.transpose()
                    << std::endl;
      ExEnv::out0() << "Vector bra center = " << uc_center_bra.transpose()
                    << std::endl;
      ExEnv::out0() << "Distance between bra and ket = " << L << std::endl;
      ExEnv::out0() << "r0_max = " << max_distance_to_refcenter_ << std::endl;
      auto cond1_val = (condition1) ? "true" : "false";
      auto cond2_val = (condition2) ? "true" : "false";
      ExEnv::out0() << "Is Condition 1 true? " << cond1_val << std::endl;
      ExEnv::out0() << "Is Condition 2 true? " << cond2_val << std::endl;
    }

    return (condition1 && condition2);
  }

 private:
  /*!
   * \brief This computes centers and extents of product density between unit
   * cell \c ref_uc and its neighbours
   */
  std::shared_ptr<detail::BasisPairInfo> compute_basis_pairs(
      const Vector3i &ref_uc = {0, 0, 0}) {
    using ::mpqc::lcao::gaussian::detail::shift_basis_origin;

    Vector3d uc_vec = ref_uc.cast<double>().cwiseProduct(dcell_);
    auto basis = shift_basis_origin(*obs_, uc_vec);
    auto basis_neighbour =
        shift_basis_origin(*obs_, uc_vec, uc_near_list_, dcell_);

    return std::make_shared<detail::BasisPairInfo>(basis, basis_neighbour);
  }
};

}  // namespace mm
}  // namespace pbc
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_MM_H_
