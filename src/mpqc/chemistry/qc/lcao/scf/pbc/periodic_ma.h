#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_MA_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_MA_H_

#include "mpqc/chemistry/qc/lcao/factory/periodic_ao_factory.h"
#include "mpqc/chemistry/qc/lcao/scf/pbc/util.h"

#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/legendre.hpp>

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
template <typename Factory,
          libint2::Operator Oper = libint2::Operator::sphemultipole>
class PeriodicMA {
 public:
  using TArray = TA::DistArray<TA::TensorD, TA::SparsePolicy>;
  using Shell = ::mpqc::lcao::gaussian::Shell;
  using ShellVec = ::mpqc::lcao::gaussian::ShellVec;
  using Basis = ::mpqc::lcao::gaussian::Basis;
  using Ord2lmMap = std::unordered_map<unsigned int, std::pair<int, int>>;
  using UnitCellList = std::vector<std::pair<Vector3i, Vector3d>>;
  using Shell2UCListMap = std::unordered_map<size_t, UnitCellList>;
  template <typename T, unsigned int nops = libint2::operator_traits<Oper>::nopers>
  using MultipoleMoment = std::array<T, nops>;
  using Shell2KernelMap = std::unordered_map<size_t, MultipoleMoment<double, (2 * MULTIPOLE_MAX_ORDER + 1) * (2 * MULTIPOLE_MAX_ORDER + 1)>>;

  /*!
   * @brief This constructs PeriodicMA using a \c PeriodicAOFactory object
   * @param ao_factory a \c PeriodicAOFactory object
   * @param ma_thresh threshold of multipole expansion error
   * @param ws well-separateness criterion
   */
  PeriodicMA(Factory &ao_factory, double ma_thresh = 1.0e-9, double bs_extent_thresh = 1.0e-6, double ws = 3.0)
      : ao_factory_(ao_factory), ma_thresh_(ma_thresh), bs_extent_thresh_(bs_extent_thresh), ws_(ws) {
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

    // determine dimensionality of crystal
    dimensionality_ = 0;
    for (auto dim = 0; dim <= 2; ++dim) {
      if (RJ_max_(dim) > 0) {
        dimensionality_++;
      }
    }
    ExEnv::out0() << "\nCrystal dimensionality : " << dimensionality_
                  << std::endl;

    // compute centers and extents of product density between the reference
    // unit cell and its neighbours
    ref_pairs_ = construct_basis_pairs();
    // compute maximum distance between the center of mass of unit cell atoms
    // and all charge centers
    ref_com_ = ao_factory_.unitcell().com();
    ExEnv::out0() << "\ncenter of mass of the reference cell is "
                  << ref_com_.transpose() << std::endl;
    max_distance_to_refcenter_ = ref_pairs_->max_distance_to(ref_com_);

    // determine CFF boundary
    cff_boundary_ = compute_CFF_boundary(RJ_max_);

    ExEnv::out0() << "\nThe boundary of Crystal Far Field is "
                  << cff_boundary_.transpose() << std::endl;

    // compute spherical multipole moments (the chargeless version, before being
    // contracted with density matrix) for electrons
    sphemm_ = ao_factory_.template compute_array<Oper>(L"<κ|O|λ>");

    // compute spherical multipole moments for nuclei (the charged version)
    O_nuc_ = compute_nuc_multipole_moments(ref_com_);

    // make a map from ordinal indices of (l, m) to (l, m) pairs
    O_ord_to_lm_map_ = make_ord_to_lm_map<MULTIPOLE_MAX_ORDER>();
    M_ord_to_lm_map_ = make_ord_to_lm_map<2 * MULTIPOLE_MAX_ORDER>();

    // test against Mathematica
//    {
//      auto M = build_interaction_kernel<2 * MULTIPOLE_MAX_ORDER>(Vector3d({-100.0, -100.0, -100.0}));
//      auto L = build_local_potential(O_nuc_, M);
//      double e = compute_energy(O_nuc_, L);
//
//      ExEnv::out0() << "\n****** TEST ******\n";
//      for (auto op = 0; op != L.size(); ++op) {
//        auto l = O_ord_to_lm_map_[op].first;
//        auto m = O_ord_to_lm_map_[op].second;
//        ExEnv::out0() << "op=" << op << ", l=" << l << ", m=" << m
//                      << ", Mlm=" << L[op] << "\n";
//      }
//      ExEnv::out0() << "MA energy for nuclei = " << e << std::endl;
//      ExEnv::out0() << "\n****** TEST COMPLETED ******\n";
//    }

  }

 private:
  Factory &ao_factory_;
  const double ma_thresh_;  /// multipole approximation is considered converged when Coulomb contribution from a spherical shell of unit cells is below this value
  const double bs_extent_thresh_;  /// threshold used in computing extent of a pair of primitives
  const double ws_;         /// well-separateness criterion

  MultipoleMoment<TArray> sphemm_;
  MultipoleMoment<double> O_nuc_;

  static constexpr unsigned int nopers_ = libint2::operator_traits<Oper>::nopers;
  static constexpr unsigned int nopers_doubled_lmax_ = (2 * MULTIPOLE_MAX_ORDER + 1) * (2 * MULTIPOLE_MAX_ORDER + 1);

  Ord2lmMap O_ord_to_lm_map_;
  Ord2lmMap M_ord_to_lm_map_;

  Shell2UCListMap cff_shell_to_unitcells_map_;
  Shell2KernelMap cff_shell_to_M_map_;
  int dimensionality_;  // dimensionality of crystal

  std::shared_ptr<Basis> obs_;
  std::shared_ptr<Basis> dfbs_;
  Vector3d dcell_;
  Vector3i R_max_;
  Vector3i RJ_max_;
  Vector3i RD_max_;
  int64_t R_size_;
  int64_t RJ_size_;
  int64_t RD_size_;
  Vector3d ref_com_;  // center of mass of the reference unit cell

  std::vector<Vector3i> uc_near_list_;
  std::shared_ptr<detail::BasisPairInfo> ref_pairs_;
  double max_distance_to_refcenter_;
  Vector3i cff_boundary_;

  double energy_cff_;
  TArray fock_cff_;
  bool energy_computed_ = false;
  bool fock_computed_ = false;

 public:
  /// @brief This returns the boundary of Crystal Far Field
  const Vector3i &CFF_boundary() { return cff_boundary_; }

  /*!
   * @brief This computes Coulomb interaction for the reference cell contributed
   * from all cells in CFF. Energy contribution is computed per spherical shell
   * from the CFF boundary to the user-given \c rjmax unless it is converged
   * @param D
   * @param target_precision
   * @return
   */
  void compute_multipole_approx(const TArray &D, double target_precision) {
    // compute electronic multipole moments for the reference unit cell
    auto O_elec_prim = compute_elec_multipole_moments(sphemm_, D);
    // store to a slightly different form
    MultipoleMoment<double> O_elec;
    for (auto op = 0; op != nopers_; ++op) {
      auto m = O_ord_to_lm_map_[op].second;
      auto sign = (m < 0 && (m % 2 == 0)) ? -1.0 : 1.0;
      O_elec[op] = sign * O_elec_prim[op];
    }

    // combine nuclear and electronic parts of multipole moments
    MultipoleMoment<double> O;
    for (auto op = 0; op != nopers_; ++op) {
      O[op] = O_elec[op] + O_nuc_[op];
    }

    MultipoleMoment<double> L;
    L.fill(0.0);

    energy_cff_ = 0.0;
    bool converged = false;
    bool out_of_rjmax = false;
    size_t cff_shell_idx = 0;  // idx of a spherical shell in CFF region
    Vector3i sphere_thickness_next = {0, 0, 0};

    do {
      auto shell_iter = cff_shell_to_unitcells_map_.find(cff_shell_idx);
      // build a list of unit cells for a spherical shell (outermost cells of
      // the sphere) if it does not exist
      if (shell_iter == cff_shell_to_unitcells_map_.end()) {
        auto unitcell_list = build_unitcells_on_a_sphere(cff_boundary_, cff_shell_idx);
        shell_iter = cff_shell_to_unitcells_map_.insert({cff_shell_idx, unitcell_list}).first;
      }

      auto M_shell_iter = cff_shell_to_M_map_.find(cff_shell_idx);
      // build interaction kernel for a spherical shell if it does not exist
      if (M_shell_iter == cff_shell_to_M_map_.end()) {
        MultipoleMoment<double, nopers_doubled_lmax_> M_shell;
        M_shell.fill(0.0);
        for (const auto &unitcell : shell_iter->second) {
          const auto &unitcell_vec = unitcell.second;

          // compute interaction kernel between multipole moments centered at the
          // reference unit cell and a remote unit cell
          auto M = build_interaction_kernel<2 * MULTIPOLE_MAX_ORDER>(ref_com_ - unitcell_vec);

          // compress M from unit cells to a spherical shell level
          for (auto op = 0; op != nopers_doubled_lmax_; ++op) {
            M_shell[op] += M[op];
          }
        }
        // save the shell-level M into a map so it can be reused
        M_shell_iter = cff_shell_to_M_map_.insert({cff_shell_idx, M_shell}).first;
      }

      // compute local potential created by all cells in a spherical shell
      auto L_shell = build_local_potential(O, M_shell_iter->second);
      // compute Coulomb energy contribution from a spherical shell
      double e_shell = compute_energy(O, L_shell);
      ExEnv::out0() << "multipole interaction energy for shell [" << cff_shell_idx << "] = " << e_shell << std::endl;

      for (auto op = 0; op != nopers_; ++op) {
        L[op] += L_shell[op];
      }
      energy_cff_ += e_shell;

      // determine if the energy is converged
      if (std::abs(e_shell) < ma_thresh_) {
        converged = true;
      }

      // determine if the next shell is out of rjmax range
      for (auto dim = 0; dim <= 2; ++dim) {
        if (RJ_max_(dim) > 0) {
          sphere_thickness_next(dim) = cff_shell_idx + 1;
        }
      }
      Vector3i corner_idx_next = cff_boundary_ + sphere_thickness_next;
      if (corner_idx_next(0) > RJ_max_(0) || corner_idx_next(1) > RJ_max_(1)
          || corner_idx_next(2) > RJ_max_(2)) {
        out_of_rjmax = true;
      }

      cff_shell_idx++;
    } while (!converged && !out_of_rjmax);

    if (!converged) {
      ExEnv::out0() << "\n!!!!!! Warning !!!!!!"
                    << "\nMultipole approximation is not converged to the given threshold!"
                    << std::endl;
    }

    ExEnv::out0() << "Coulomb energy contributed from CFF so far = " << energy_cff_ << std::endl;

    // compute Fock contribution from CFF
    fock_cff_ = sphemm_[0];
    fock_cff_("mu, nu") = -1.0 * L[0] * fock_cff_("mu, nu");
    if (MULTIPOLE_MAX_ORDER >= 1) {
      for (auto op = 1; op != nopers_; ++op) {
        auto l = O_ord_to_lm_map_[op].first;
        auto m = O_ord_to_lm_map_[op].second;
        auto signl = (l % 2 == 0) ? 1.0 : -1.0;
        auto signm = (m < 0 && (m % 2 != 0)) ? -1.0 : 1.0;
        auto prefactor = -1.0 * signl * signm * L[op];
        fock_cff_("mu, nu") += prefactor * sphemm_[op]("mu, nu");
      }
    }
    fock_computed_ = true;
    energy_computed_ = true;
  }

  /// compute electronic multipole moments for the reference unit cell
  MultipoleMoment<double> compute_elec_multipole_moments(const TArray &D) {
    return compute_elec_multipole_moments(sphemm_, D);
  }

  /// return whether Fock contributed from CFF is computed
  bool fock_computed() {return fock_computed_; }

  /// return whether Coulomb interaction energy contributed from CFF is computed
  bool energy_computed() {return energy_computed_; }

  const TArray &get_fock() { return fock_cff_; }
  
  double get_energy() { return energy_cff_; }

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
                                                   bs_extent_thresh_);
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
          ExEnv::out0() << "\n!!!!!! Warning !!!!!!"
                        << "\nThe range limit is not enough to reach Crystal Far Field. Use larger `rjmax` in the input."
                        << std::endl;
        }
      }
    }

    return cff_bound;
  }

  /*!
   * @brief This makes a map from the ordinal index of a (\em l, \em m) pair to
   * the corresponding \em l and \em m
   * @tparam lmax max value of l
   * @return \c unordered_map with
   * key: ordinal index
   * mapped value: (\em l, \em m) pair
   */
  template <unsigned int lmax>
  Ord2lmMap make_ord_to_lm_map() {
    Ord2lmMap result;
    result.reserve((lmax + 1) * (lmax + 1));

    auto ord_idx = 0u;
    for (auto l = 0; l <= lmax; ++l) {
      for (auto m = -l; m <= l; ++m) {
        result[ord_idx] = std::make_pair(l, m);
        ord_idx++;
      }
    }

    return result;
  }

  /*!
   * @brief This builds the real interaction kernel \latexonly M$_{l, m}$
   * \endlatexonly between two distant multipole moments centered at P and Q,
   * respectively
   *
   * \latexonly
   * \begin{eqnarray*}
   * M_{l,m}(\mathbf{r}) =
   * \begin{cases}
   *  \frac{(l-m)!}{|\mathbf{r}|^{l+1}} P_l^m(cos\theta) cos(m \phi),
   *    \text{ if m $\geq$ 0} \\
   *  \frac{(l-m)!}{|\mathbf{r}|^{l+1}} P_l^m(cos\theta) sin(m \phi),
   *    \text{ if m $<$ 0}
   * \end{cases}
   * \end{eqnarray*}
   * \endlatexonly
   *
   * @tparam lmax max value of l
   * @param c12 non-zero vector from center Q to center P
   * @return an array of doubles with array size = number of components in
   * \c Oper
   */
  template<unsigned int lmax>
  MultipoleMoment<double, (lmax + 1) * (lmax + 1)> build_interaction_kernel(const Vector3d &c12) {
    const auto c12_norm = c12.norm();
    MPQC_ASSERT(c12_norm > 0.0);

    using namespace boost::math;
    MultipoleMoment<double, (lmax + 1) * (lmax + 1)> result;

    const auto cos_theta = c12(2) / c12_norm;
    const auto c12_xy_norm = std::sqrt(c12(0) * c12(0) + c12(1) * c12(1));

    // When Rx^2 + Ry^2 = 0, φ cannot computed. We use the fact that in this
    // case cosθ = 1, and then associated Legendre polynomials
    // P_{l, m}(cosθ) = 0 if m != 0
    // P_{l, m}(cosθ) = 1 if m == 0
    if (c12_xy_norm < std::numeric_limits<double>::epsilon()) {
      for (auto l = 0, ord_idx = 0; l <= lmax; ++l) {
        auto frac = factorial<double>(l) / std::pow(c12_norm, l + 1);
        for (auto m = -l; m <= l; ++m, ++ord_idx) {
          result[ord_idx] = (m == 0) ? frac : 0.0;
        }
      }
      return result;
    } else {
      // build a map that returns cos(m * phi) for each m using recurrence relation
      auto cos_m_phi_map = cos_m_phi_map_recursive<lmax>(c12);
      // build a map that returns sin(m * phi) for each m using recurrence relation
      auto sin_m_phi_map = sin_m_phi_map_recursive<lmax>(c12);
      // build a map that returns associated Legendre polynomial P_{l, m}(cosθ)
      // for each ordinal index of (l, m) pair, using recurrence relation
      auto legendre_map = associated_legendre_map_recursive<lmax>(cos_theta);

      // fill the result
      for (auto l = 0, ord_idx = 0; l <= lmax; ++l) {
        const auto inv_denom = 1.0 / std::pow(c12_norm, l + 1);
        for (auto m = -l; m <= l; ++m, ++ord_idx) {
          const auto num = factorial<double>(l - m);
          const auto phi_part = (m >= 0) ? cos_m_phi_map[m] : sin_m_phi_map[m];
          result[ord_idx] = num * inv_denom * legendre_map[ord_idx] * phi_part;
        }
      }

      return result;
    }
  }

  /*!
   * @brief This builds the real multipole expansion \latexonly O$_{l, m}$
   * \endlatexonly for a vector \em r with respect to the origin
   *
   * \latexonly
   * \begin{eqnarray*}
   * O_{l,m}(\mathbf{r}) =
   * \begin{cases}
   *  \frac{|\mathbf{r}|^{l}}{(l+m)!} P_l^m(cos\theta) cos(m \phi),
   *    \text{ if m $\geq$ 0} \\
   *  \frac{|\mathbf{r}|^{l}}{(l+m)!} P_l^m(cos\theta) sin(m \phi),
   *    \text{ if m $<$ 0}
   * \end{cases}
   * \end{eqnarray*}
   * \endlatexonly
   *
   * @tparam lmax max value of l
   * @param r non-zero vector from the origin to \em r
   * @return an array of doubles with array size = number of components in
   * \c Oper
   */
  template<unsigned int lmax>
  MultipoleMoment<double, (lmax + 1) * (lmax + 1)> build_multipole_expansion(const Vector3d &r) {
    const auto r_norm = r.norm();

    MultipoleMoment<double, (lmax + 1) * (lmax + 1)> result;

    if (r_norm < std::numeric_limits<double>::epsilon()) {
      result.fill(0.0);
      return result;
    } else {
      using namespace boost::math;

      const auto cos_theta = r(2) / r_norm;
      const auto r_xy_norm = std::sqrt(r(0) * r(0) + r(1) * r(1));

      // When Rx^2 + Ry^2 = 0, φ cannot computed. We use the fact that in this
      // case cosθ = 1, and then associated Legendre polynomials
      // P_{l, m}(cosθ) = 0 if m != 0
      // P_{l, m}(cosθ) = 1 if m == 0
      if (r_xy_norm < std::numeric_limits<double>::epsilon()) {
        for (auto l = 0, ord_idx = 0; l <= lmax; ++l) {
          auto frac = std::pow(r_norm, l) / factorial<double>(l);
          for (auto m = -l; m <= l; ++m, ++ord_idx) {
            result[ord_idx] = (m == 0) ? frac : 0.0;
          }
        }
        return result;
      } else {
        // build a map that returns cos(m * phi) for each m using recurrence relation
        auto cos_m_phi_map = cos_m_phi_map_recursive<lmax>(r);
        // build a map that returns sin(m * phi) for each m using recurrence relation
        auto sin_m_phi_map = sin_m_phi_map_recursive<lmax>(r);
        // build a map that returns associated Legendre polynomial P_{l, m}(cosθ)
        // for each ordinal index of (l, m) pair, using recurrence relation
        auto legendre_map = associated_legendre_map_recursive<lmax>(cos_theta);

        // fill the result
        for (auto l = 0, ord_idx = 0; l <= lmax; ++l) {
          const auto num = std::pow(r_norm, l);
          for (auto m = -l; m <= l; ++m, ++ord_idx) {
            const auto inv_denom = 1.0 / factorial<double>(l + m);
            const auto phi_part = (m >= 0) ? cos_m_phi_map[m] : sin_m_phi_map[m];
            result[ord_idx] = num * inv_denom * legendre_map[ord_idx] * phi_part;
          }
        }

        return result;
      }
    }
  }

  /*!
   * @brief This builds a map that returns sin(m * phi) for each m using recurrence relation
   * @tparam lmax
   * @param r
   * @return
   */
  template<unsigned int lmax>
  std::unordered_map<int, double> sin_m_phi_map_recursive(const Vector3d &r) {
    const auto r_xy_norm = std::sqrt(r(0) * r(0) + r(1) * r(1));
    MPQC_ASSERT(r_xy_norm >= std::numeric_limits<double>::epsilon());

    const auto inv_x2py2 = 1.0 / r_xy_norm;
    const auto cos_phi = r(0) * inv_x2py2;
    const auto sin_phi = r(1) * inv_x2py2;

    std::unordered_map<int, double> sin_m_phi_map;
    sin_m_phi_map.reserve(2 * lmax + 1);

    sin_m_phi_map[0] = 0.0;
    sin_m_phi_map[0] = 0.0;
    if (lmax >= 1) {
      sin_m_phi_map[1] = sin_phi;
      sin_m_phi_map[-1] = -sin_phi;
    }
    if (lmax >= 2) {
      for (auto m = 2; m <= lmax; ++m) {
        sin_m_phi_map[m] = 2.0 * sin_m_phi_map[m - 1] * cos_phi - sin_m_phi_map[m - 2];
        sin_m_phi_map[-m] = -sin_m_phi_map[m];
      }
    }

    return sin_m_phi_map;
  }

  /*!
   * @brief This builds a map that returns cos(m * phi) for each m using recurrence relation
   * @tparam lmax
   * @param r
   * @return
   */
  template<unsigned int lmax>
  std::unordered_map<int, double> cos_m_phi_map_recursive(const Vector3d &r) {
    const auto r_xy_norm = std::sqrt(r(0) * r(0) + r(1) * r(1));
    MPQC_ASSERT(r_xy_norm >= std::numeric_limits<double>::epsilon());

    const auto cos_phi = r(0) / r_xy_norm;

    std::unordered_map<int, double> cos_m_phi_map;
    cos_m_phi_map.reserve(2 * lmax + 1);

    cos_m_phi_map[0] = 1.0;
    if (lmax >= 1) {
      cos_m_phi_map[1] = cos_phi;
      cos_m_phi_map[-1] = cos_phi;
    }
    if (lmax >= 2) {
      for (auto m = 2; m <= lmax; ++m) {
        cos_m_phi_map[m] = 2.0 * cos_m_phi_map[m - 1] * cos_phi - cos_m_phi_map[m - 2];
        cos_m_phi_map[-m] = cos_m_phi_map[m];
      }
    }

    return cos_m_phi_map;
  }

  /*!
   * This builds a map that returns associated Legendre polynomial P_{l, m}(cosθ)
   * for each ordinal index of (l, m) pair, using recurrence relation
   * @tparam lmax
   * @param cos_theta
   * @return
   */
  template<unsigned int lmax>
  std::unordered_map<int, double> associated_legendre_map_recursive(double cos_theta) {
    using namespace boost::math;

    double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);

    std::unordered_map<int, double> legendre_map;
    legendre_map.reserve((lmax + 1) * (lmax + 1));

    legendre_map[0] = 1.0;  // (l = 0, m = 0)
    if (lmax >= 1) {
      // compute P_{l, m} for all l >= 1 and m = 0
      legendre_map[2] = cos_theta;  // (l = 1, m = 0)
      if (lmax >= 2) {
        for (auto l = 2; l <= lmax; ++l) {
          const auto ord_l_0 = l * l + l;  // (l, 0)
          const auto ord_lm1_0 = l * l - l;  // (l - 1, 0)
          const auto ord_lm2_0 = l * l - 3 * l + 2;  // (l - 2, 0)
          legendre_map[ord_l_0] = legendre_next(l - 1, 0, cos_theta, legendre_map[ord_lm1_0], legendre_map[ord_lm2_0]);
        }
      }

      // compute P_{l, m} for all l >= 1 and m != 0
      for (auto m = 1; m <= lmax; ++m) {
        const auto sign = (m % 2 == 0) ? 1 : -1;
        const auto m2 = m * m;

        const auto ord_m_m = m2 + m + m;  // (l = m, m)
        const auto ord_mm1_mm1 = m2 - 1;  // (l = m - 1, m - 1)
        legendre_map[ord_m_m] = -1.0 * (2.0 * m - 1.0) * sin_theta * legendre_map[ord_mm1_mm1];
        legendre_map[m2] = sign * legendre_map[ord_m_m] / factorial<double>(m + m);  // (l = m, -m)

        if (lmax >= 2 && lmax >= m + 1) {
          const auto ord_mp1_m = m2 + 4 * m + 2;  // (l = m + 1, m)
          legendre_map[ord_mp1_m] = (2.0 * m + 1.0) * cos_theta * legendre_map[ord_m_m];
          legendre_map[ord_mp1_m - 2 * m] = sign * legendre_map[ord_mp1_m] / factorial<double>(2 * m + 1);  // (l = m + 1, -m)
        }

        if (lmax >= 3 && lmax >= m + 2) {
          for (auto l = m + 2; l <= lmax; ++l) {
            const auto l2 = l * l;
            const auto ord_l_m = l2 + l + m;  // (l, m)
            const auto ord_lm1_m = l2 - l + m;  // (l - 1, m)
            const auto ord_lm2_m = l2 - 3 * l + 2 + m;  // (l - 2, m)
            legendre_map[ord_l_m] = legendre_next(l - 1, m, cos_theta, legendre_map[ord_lm1_m], legendre_map[ord_lm2_m]);
            legendre_map[ord_l_m - 2 * m] = sign * factorial<double>(l - m) * legendre_map[ord_l_m] / factorial<double>(l + m);  // (l, -m)
          }
        }
      }
    }

    return legendre_map;
  }

  /*!
   * @brief This forms real multipole moments with cos(m φ) and sin(m φ) parts,
   * respectively.
   * @tparam nopers total number of componments in multipole expansion operator
   * @param M multipole moments with mixed cosine and sine parts
   * @param ord_to_lm_map a map from ordinal index of (l, m) to a specific
   * (l, m) pair
   * @return a pair of multipole moments with cos(m φ) and sin(m φ) parts,
   * respectively.
   */
  template <unsigned int nopers>
  std::pair<MultipoleMoment<double, nopers>, MultipoleMoment<double, nopers>>
  form_symm_and_antisymm_moments(const MultipoleMoment<double, nopers> &M, Ord2lmMap &ord_to_lm_map) {
    assert(M.size() == nopers);
    assert(ord_to_lm_map.size() == nopers);
    MultipoleMoment<double, nopers> Mplus, Mminus;
    for (auto op = 0u; op != nopers; ++op) {
      auto m = ord_to_lm_map[op].second;
      auto sign_of_m = (m > 0) ? 1 : ((m < 0) ? -1 : 0);
      auto nega1_m = (m % 2 == 0) ? 1.0 : -1.0;  // (-1)^m
      switch (sign_of_m) {
        case 1:
          Mplus[op] = M[op];
          Mminus[op] = -nega1_m * M[op - 2 * m];
          break;
        case 0:
          Mplus[op] = M[op];
          Mminus[op] = 0.0;
          break;
        case -1:
          Mplus[op] = nega1_m * M[op - 2 * m];
          Mminus[op] = M[op];
          break;
        default:
          throw ProgrammingError("Invalid sign of m.", __FILE__, __LINE__);
      }
    }
    return std::make_pair(Mplus, Mminus);
  }

  /*!
   * @brief This builds the local multipole expansion of the potential at the origin
   * created by a remote unit cell
   * @param O multipole moments of the charge distribution of the remote unit cell
   * @param M multipole interaction kernel between the origin and the remote unit cell
   * @return
   */
  MultipoleMoment<double> build_local_potential(const MultipoleMoment<double>& O,
                                                const MultipoleMoment<double, nopers_doubled_lmax_> &M) {
    MultipoleMoment<double> result;

    // form O+, O-
    auto O_pm = form_symm_and_antisymm_moments<nopers_>(O, O_ord_to_lm_map_);
    auto &O_plus = O_pm.first;
    auto &O_minus = O_pm.second;

    // form M+, M-
    auto M_pm = form_symm_and_antisymm_moments<nopers_doubled_lmax_>(M, M_ord_to_lm_map_);
    auto &M_plus = M_pm.first;
    auto &M_minus = M_pm.second;

    // form multipole expansion of local potential created by the remote unit cell
    for (auto op1 = 0; op1 != nopers_; ++op1) {
      auto l1 = O_ord_to_lm_map_[op1].first;
      auto m1 = O_ord_to_lm_map_[op1].second;
      auto accu = 0.0;
      for (auto op2 = 0; op2 != nopers_; ++op2) {
        auto l1pl2 = l1 + O_ord_to_lm_map_[op2].first;
        auto m1pm2 = m1 + O_ord_to_lm_map_[op2].second;
        auto ord_M = l1pl2 * l1pl2 + l1pl2 + m1pm2;
        if (m1 >= 0) {
          accu += M_plus[ord_M] * O_plus[op2];
          accu += M_minus[ord_M] * O_minus[op2];
        } else {
          accu += M_minus[ord_M] * O_plus[op2];
          accu -= M_plus[ord_M] * O_minus[op2];
        }
      }
      result[op1] = accu;
    }

    return result;
  }

  /*!
   * @brief This computes multipole interaction energy for give multipole moments and local potential
   * @param O multipole moments
   * @param L local potential created by distant charges
   * @return
   */
  double compute_energy(const MultipoleMoment<double> &O,
                        const MultipoleMoment<double> &L) {
    double result = 0.0;
    for (auto op = 0; op != nopers_; ++op) {
      auto l = O_ord_to_lm_map_[op].first;
      auto m = O_ord_to_lm_map_[op].second;
      auto sign = (l % 2 == 0) ? 1.0 : -1.0;
      auto delta = (m == 0) ? 1.0 : 2.0;
      // scale by 0.5 because energy per cell is half the interaction energy
      result += 0.5 * sign * delta * O[op] * L[op];
    }

    return result;
  }

  /*!
   * @brief This builds a list a unit cells in a give spherical shell (outermost
   * cells of the sphere)
   * @param sphere_start
   * @param shell_idx
   * @return
   */
  UnitCellList build_unitcells_on_a_sphere(const Vector3i &sphere_start,
                                           size_t shell_idx) {
    UnitCellList result;

    // compute thickness of the sphere to CFF boundary
    Vector3i sphere_thickness = {0, 0, 0};
    for (auto dim = 0; dim <= 2; ++dim) {
      if (RJ_max_(dim) > 0) {
        sphere_thickness(dim) = shell_idx;
      }
    }

    // index (x, y, z) of the unit cell at the corner with positive x, y, z
    Vector3i corner_idx = sphere_start + sphere_thickness;

    // a lambda that returns a (idx, vector) pair for a given idx of a unit cell
    auto make_pair = [](const Vector3i &idx, const Vector3d &dcell) {
      return std::make_pair(idx, idx.cast<double>().cwiseProduct(dcell));
    };

    // build the unit cell list depending on the crystal dimensionality
    if (dimensionality_ == 0) {
      result.emplace_back(std::make_pair(Vector3i::Zero(), Vector3d::Zero()));
    } else if (dimensionality_ == 1) {
      auto dim = (RJ_max_(0) > 0) ? 0 : ((RJ_max_(1) > 0) ? 1 : 2);
      MPQC_ASSERT(corner_idx(dim) > 0);
      Vector3i neg1_idx = -1 * corner_idx;
      result.emplace_back(make_pair(corner_idx, dcell_));
      result.emplace_back(make_pair(neg1_idx, dcell_));
    } else if (dimensionality_ == 2) {
      auto dim_a = (RJ_max_(0) == 0) ? 1 : 0;
      auto dim_b = (RJ_max_(2) == 0) ? 1 : 2;
      const auto a_max = corner_idx(dim_a);
      const auto b_max = corner_idx(dim_b);
      MPQC_ASSERT(a_max > 0 && b_max > 0);
      Vector3i pos2_idx = Vector3i::Zero();
      Vector3i neg2_idx = Vector3i::Zero();
      pos2_idx(dim_a) = a_max;
      neg2_idx(dim_a) = -a_max;
      for (auto b = -b_max; b <= b_max; ++b) {
        pos2_idx(dim_b) = b;
        neg2_idx(dim_b) = b;
        result.emplace_back(make_pair(pos2_idx, dcell_));
        result.emplace_back(make_pair(neg2_idx, dcell_));
      }

      pos2_idx.setZero();
      neg2_idx.setZero();
      pos2_idx(dim_b) = b_max;
      neg2_idx(dim_b) = -b_max;
      for (auto a = -a_max + 1; a <= a_max - 1; ++a) {
        pos2_idx(dim_a) = a;
        neg2_idx(dim_a) = a;
        result.emplace_back(make_pair(pos2_idx, dcell_));
        result.emplace_back(make_pair(neg2_idx, dcell_));
      }
    } else if (dimensionality_ == 3) {
      const auto x_max = corner_idx[0];
      const auto y_max = corner_idx[1];
      const auto z_max = corner_idx[2];
      MPQC_ASSERT(x_max > 0 && y_max > 0 && z_max > 0);
      Vector3i pos3_idx = Vector3i::Zero();
      Vector3i neg3_idx = Vector3i::Zero();
      pos3_idx(0) = x_max;
      neg3_idx(0) = -x_max;
      for (auto y = -y_max; y <= y_max; ++y) {
        pos3_idx(1) = y;
        neg3_idx(1) = y;
        for (auto z = -z_max; z <= z_max; ++z) {
          pos3_idx(2) = z;
          neg3_idx(2) = z;
          result.emplace_back(make_pair(pos3_idx, dcell_));
          result.emplace_back(make_pair(neg3_idx, dcell_));
        }
      }

      pos3_idx.setZero();
      neg3_idx.setZero();
      pos3_idx(1) = y_max;
      neg3_idx(1) = -y_max;
      for (auto x = -x_max + 1; x <= x_max - 1; ++x) {
        pos3_idx(0) = x;
        neg3_idx(0) = x;
        for (auto z = -z_max; z <= z_max; ++z) {
          pos3_idx(2) = z;
          neg3_idx(2) = z;
          result.emplace_back(make_pair(pos3_idx, dcell_));
          result.emplace_back(make_pair(neg3_idx, dcell_));
        }
      }

      pos3_idx.setZero();
      neg3_idx.setZero();
      pos3_idx(2) = z_max;
      neg3_idx(2) = -z_max;
      for (auto x = -x_max + 1; x <= x_max - 1; ++x) {
        pos3_idx(0) = x;
        neg3_idx(0) = x;
        for (auto y = -y_max + 1; y <= y_max - 1; ++y) {
          pos3_idx(1) = y;
          neg3_idx(1) = y;
          result.emplace_back(make_pair(pos3_idx, dcell_));
          result.emplace_back(make_pair(neg3_idx, dcell_));
        }
      }
    } else {
      throw ProgrammingError("Invalid crystal dimensionality.", __FILE__, __LINE__);
    }

    return result;
  }

  /*!
   * @brief This computes electronic multipole moments for the reference unit cell
   * @param sphemm
   * @param D
   * @return
   */
  MultipoleMoment<double> compute_elec_multipole_moments(
      const MultipoleMoment<TArray> &sphemm, const TArray &D) {
    MultipoleMoment<double> result;
    using ::mpqc::pbc::detail::dot_product;
    for (auto op = 0; op != nopers_; ++op) {
      result[op] = -2.0 * dot_product(sphemm[op], D, R_max_, RD_max_);  // 2 = alpha + beta, -1 = electron charge
    }

    return result;
  }

  /*!
   * @brief This computes nuclear multipole moments for the reference unit cell
   * @param center expansion center
   * @return
   */
  MultipoleMoment<double> compute_nuc_multipole_moments(const Vector3d &center) {
    const auto &atoms = ao_factory_.unitcell().atoms();
    const auto natoms = atoms.size();

    std::unordered_map<unsigned int, MultipoleMoment<double>> O_atoms;

    for (auto i = 0u; i != natoms; ++i) {
      const Vector3d &r = atoms[i].center() - center;
      O_atoms[i] = build_multipole_expansion<MULTIPOLE_MAX_ORDER>(r);
    }

    MultipoleMoment<double> result;
    for (auto op = 0; op != nopers_; ++op) {
      double accu = 0.0;
      for (auto i = 0; i != natoms; ++i) {
        accu += O_atoms[i][op] * atoms[i].charge();
      }
      result[op] = accu;
    }

    return result;
  }

};

}  // namespace ma
}  // namespace pbc
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_PBC_PERIODIC_MA_H_
