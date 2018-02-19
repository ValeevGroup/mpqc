#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_MA_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_MA_H_

#include "mpqc/chemistry/qc/lcao/factory/periodic_ao_factory.h"
#include "mpqc/chemistry/qc/lcao/scf/pbc/util.h"

namespace mpqc {
namespace pbc {

namespace detail {

/*!
 * @brief Class to hold information for basis shell pairs, including total
 * number of shell pairs, extents (radii) and centers of shell pairs, etc.
 */
class BasisPairInfo {
 public:
  using Shell = ::mpqc::lcao::gaussian::Shell;
  using Basis = ::mpqc::lcao::gaussian::Basis;

  BasisPairInfo() = default;

  /*!
   * @brief This constructs BasisPairInfo between two basis sets
   * @param bs0 first basis set
   * @param bs1 second basis set
   * @param thresh threshold of the first-order multipole expansion error
   * @param small_extent a small value that is used when no real solution exists
   * for any shell pair
   */
  BasisPairInfo(std::shared_ptr<const Basis> bs0,
                std::shared_ptr<const Basis> bs1, const double thresh = 1.0e-6,
                const double small_extent = 0.01);

 private:
  const double thresh_;        ///> threshold of multipole expansion error
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
  /*!
   * @brief This returns the extent (radius) of a pair of shells (shell \em i
   * and shell \em j)
   * @param i index of shell \em i
   * @param j index of shell \em j
   * @return the extent (radius)
   */
  double extent(int64_t i, int64_t j) const { return pair_extents_(i, j); }

  /*!
   * @brief This returns the extent (radius) of a pair of shells (shell \em i
   * and shell \em j)
   * @param ij the ordinal index of the shell pair (\em i, \em j)
   * @return the extent (radius)
   */
  double extent(int64_t ij) const;

  /*!
   * @brief This returns the weighted center of a pair of shells (shell \em i
   * and shell \em j)
   * @param i index of shell \em i
   * @param j index of shell \em j
   * @return the weighted center
   */
  Vector3d center(int64_t i, int64_t j) const { return pair_centers_[i][j]; }

  /*!
   * @brief This returns the weighted center of a pair of shells (shell \em i
   * and shell \em j)
   * @param ij the ordinal index of the shell pair (\em i, \em j)
   * @return the weighted center
   */
  Vector3d center(int64_t ij) const;

  /*!
   * @brief This returns the total number of shell pairs
   */
  size_t npairs() const { return npairs_; }

  /*!
   * \brief This computes maximum distance between \c ref_point and all charge
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
   * \brief This computes the extent (radius) of the product of two primitives
   * with exponents \c exp0 and \c exp1, respectively. \param exp0 exponent of
   * the first primitive function \param exp1 exponent of the second primitive
   * function \param rab distance between centers of two primitives \return the
   * extent (radius)
   */
  double prim_pair_extent(const double exp0, const double exp1,
                          const double rab);
};

}  // namespace detail

namespace ma {

/*!
 * \brief This class computes the contribution to Coulomb interaction from unit
 * cells in Crystal Far Field (CFF) using multipole approximation.
 */
template <typename Factory, libint2::Operator Oper = libint2::Operator::sphemultipole>
class PeriodicMA {
 public:
  using TArray = TA::DistArray<TA::TensorD, TA::SparsePolicy>;
  using Shell = ::mpqc::lcao::gaussian::Shell;
  using ShellVec = ::mpqc::lcao::gaussian::ShellVec;
  using Basis = ::mpqc::lcao::gaussian::Basis;
  template<typename T>
  using MultipoleMoment = std::array<T, libint2::operator_traits<Oper>::nopers>;
  /*!
   * @brief This constructs PeriodicMA using a \c PeriodicAOFactory object
   * @param ao_factory a \c PeriodicAOFactory object
   * @param ma_thresh threshold of multipole expansion error
   * @param ws well-separateness criterion
   */
  PeriodicMA(Factory &ao_factory, double ma_thresh = 1.0e-6, double ws = 3.0)
      : ao_factory_(ao_factory), ma_thresh_(ma_thresh), ws_(ws) {
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

    // compute centers and extents of product density between the reference
    // unit cell and its neighbours
    ref_pairs_ = construct_basis_pairs();
    // compute maximum distance between the center of mass of unit cell atoms
    // and all charge centers
    const auto &ref_com = ao_factory_.unitcell().com();
    max_distance_to_refcenter_ = ref_pairs_->max_distance_to(ref_com);

    // determine CFF boundary
    cff_boundary_ = compute_CFF_boundary(RJ_max_);

    // compute spherical multipole moments (the chargeless version, before being
    // contracted with density matrix)
    sphemm_ = ao_factory_.template compute_array<Oper>(L"<κ|O|λ>");

    // test
//    ExEnv::out0() << "\nspherical multiple:\n";
//    for (auto oper = 0u; oper != sphemm.size(); ++oper) {
//      ExEnv::out0() << "when oper = " << oper << ", sphe multipole = \n"
//                    << sphemm[oper] << std::endl;
//    }

    ExEnv::out0() << "\nThe boundary of Crystal Far Field is "
                  << cff_boundary_.transpose() << std::endl;
  }

 private:
  Factory &ao_factory_;
  const double ma_thresh_;  /// threshold of multipole expansion error
  const double ws_;         /// well-separateness criterion

  MultipoleMoment<TArray> sphemm_;

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
  Vector3i cff_boundary_;

 public:
  const Vector3i &CFF_boundary() { return cff_boundary_; }

  MultipoleMoment<double> compute_multipole_moments(const TArray &D, double target_precision) {
    std::array<double, 9> moments;
    for (auto op = 0; op != moments.size(); ++op) {
      moments[op] = -2.0 * sphemm_[op]("mu, nu") * D("mu, nu");
    }
    ExEnv::out0() << "\n****** electronic multipole moments ******\n"
                  << "monopole = "
                  << moments[0] << "\n"
                  << "dipole = ("
                  << moments[1] << ", "
                  << moments[2] << ", "
                  << moments[3] << ")\n"
                  << "quadrupole = ("
                  << moments[4] << ", "
                  << moments[5] << ", "
                  << moments[6] << ", "
                  << moments[7] << ", "
                  << moments[8] << ")\n";
  }

 private:
  /*!
   * \brief This constructs a \c BasisPairInfo object between basis sets of unit
   * cell \c ref_uc and its neighbours
   */
  std::shared_ptr<detail::BasisPairInfo> construct_basis_pairs(
      const Vector3i &ref_uc = {0, 0, 0}) {
    using ::mpqc::lcao::gaussian::detail::shift_basis_origin;

    Vector3d uc_vec = ref_uc.cast<double>().cwiseProduct(dcell_);
    auto basis = shift_basis_origin(*obs_, uc_vec);
    auto basis_neighbour = shift_basis_origin(*obs_, uc_vec, R_max_, dcell_);

    return std::make_shared<detail::BasisPairInfo>(basis, basis_neighbour,
                                                   ma_thresh_);
  }

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
    // CFF condition #1: all charge distributions are well separated
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
    const auto L = vec_rel.norm();
    bool condition2 = (L >= ws_ * (max_distance_to_refcenter_ * 2.0));

    // test:
    //    {
    //      ExEnv::out0() << "Vector bra corner = " << vec_bra.transpose()
    //                    << std::endl;
    //      ExEnv::out0() << "Vector ket corner = " << vec_ket.transpose()
    //                    << std::endl;
    //      ExEnv::out0() << "Vector bra center = " << uc_center_bra.transpose()
    //                    << std::endl;
    //      ExEnv::out0() << "Distance between bra and ket = " << L <<
    //      std::endl; ExEnv::out0() << "r0_max = " <<
    //      max_distance_to_refcenter_ << std::endl; auto cond1_val =
    //      (condition1) ? "true" : "false"; auto cond2_val = (condition2) ?
    //      "true" : "false"; ExEnv::out0() << "Is Condition 1 true? " <<
    //      cond1_val << std::endl; ExEnv::out0() << "Is Condition 2 true? " <<
    //      cond2_val << std::endl;
    //    }

    return (condition1 && condition2);
  }

  /*!
   * @brief This computes Crystal Far Field (CFF) boundary (bx, by, bz)
   * @param limit3d the range limit (lx, ly, lz) of CFF boundary so that
   * bx <= lx, by <= ly, bz <= lz. The computation of the boundary in one
   * dimension will be skipped if its limit is zero.
   * @return CFF boundary
   */
  Vector3i compute_CFF_boundary(const Vector3i &limit3d) {
    Vector3i cff_bound({0, 0, 0});

    for (auto dim = 0; dim <= 2; ++dim) {
      Vector3i uc_idx({0, 0, 0});
      bool is_in_CFF = false;
      auto idx1 = 0;
      if (limit3d(dim) > 0) {
        do {
          uc_idx(dim) = idx1;
          is_in_CFF = is_uc_in_CFF(uc_idx);
          idx1++;
        } while (idx1 <= limit3d(dim) && !is_in_CFF);

        cff_bound(dim) = uc_idx(dim);

        if (!is_in_CFF) {
          throw AlgorithmException(
              "Insufficient range limit for the boundary of Crystal Far Field",
              __FILE__, __LINE__);
        }
      }
    }

    return cff_bound;
  }
};

}  // namespace ma
}  // namespace pbc
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_MA_H_
