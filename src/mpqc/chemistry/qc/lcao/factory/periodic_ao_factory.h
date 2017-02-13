#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_LCAO_FACTORY_PERIODIC_AO_FACTORY_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_LCAO_FACTORY_PERIODIC_AO_FACTORY_H_

#include "mpqc/chemistry/qc/lcao/factory/ao_factory.h"
#include "mpqc/chemistry/qc/lcao/factory/factory_utility.h"

#include <iosfwd>
#include <vector>

#include "mpqc/chemistry/molecule/unit_cell.h"
#include "mpqc/chemistry/qc/lcao/integrals/integrals.h"
#include "mpqc/chemistry/units/units.h"
#include "mpqc/math/external/eigen/eigen.h"
#include "mpqc/math/tensor/clr/array_to_eigen.h"
#include "mpqc/util/keyval/keyval.h"
#include "mpqc/util/misc/time.h"

#include <mpqc/chemistry/qc/lcao/integrals/screening/schwarz_screen.h>
#include <unsupported/Eigen/MatrixFunctions>

typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic,
                      Eigen::RowMajor>
    Matrixz;
typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> Vectorz;

// constant
const std::complex<double> I(0.0, 1.0);

namespace mpqc {
namespace lcao {

namespace detail {
/*!
 * \brief This extends 1D tiled range by repeating it multiple times
 *
 * \param tr0 the original TiledRange1 object
 * \param size the number of times for repeating
 * \return the extended TiledRange1 object
 */
TA::TiledRange1 extend_trange1(TA::TiledRange1 tr0, int64_t size);

/*!
 * \brief This sorts eigenvalues and eigenvectors
 * in ascending order of the real parts of eigenvalues
 *
 * \param eigVal the vector of complex eigenvalues
 * \param eigVec the complex matrix consisting of complex eigenvectors
 */
void sort_eigen(Vectorz &eigVal, Matrixz &eigVec);

/*!
 * \brief This takes the ordinal index of a lattice
 * and returns the corresponding direct-space lattice vector
 *
 * \param ord_index the ordinal index of the lattice
 * \param latt_max the range of included lattices
 * \param dcell the direct unit cell params
 * \return the direct-space lattice vector
 */
Vector3d direct_vector(int64_t ord_idx, Vector3i latt_max, Vector3d dcell);

/*!
 * \brief This takes the ordinal index of a k-space (reciprocal-space) lattice
 * and returns the corresponding k-space lattice vector
 *
 * \param ord_idx the ordinal index of the reciprocal lattice
 * \param nk the range of included k points
 * \param dcell the direct unit cell params
 * \return the k-space lattice vector
 */
Vector3d k_vector(int64_t ord_idx, Vector3i nk, Vector3d dcell);

/*!
 * \brief This takes the 3D index of a direct lattice
 * and returns the corresponding ordinal index
 *
 * \param x the direct lattice index on x axis
 * \param y the direct lattice index on y axis
 * \param z the direct lattice index on z axis
 * \param latt_max the range of included lattices
 * \return the ordinal index in direct space
 */
int64_t direct_ord_idx(int64_t x, int64_t y, int64_t z, Vector3i latt_max);

/*!
 * \brief This takes the 3D index of a reciprocal lattice
 * and returns the corresponding ordinal index
 *
 * \param x the reciprocal lattice index on x' axis
 * \param y the reciprocal lattice index on y' axis
 * \param z the reciprocal lattice index on z' axis
 * \param nk the range of included k points
 * \return the ordinal index in k space
 */
int64_t k_ord_idx(int64_t x, int64_t y, int64_t z, Vector3i nk);

/*!
 * \brief This shifts the position of a Molecule object
 *
 * \note All atom positions will be shifted
 * \param mol the Molecule object
 * \param shift the 3D vector of the shift
 * \return the shared pointer of the shifted Molecule object
 */
std::shared_ptr<Molecule> shift_mol_origin(const Molecule &mol, Vector3d shift);

}  // namespace mpqc::lcao::detail

namespace gaussian {

template <typename Tile, typename Policy>
class PeriodicAOFactory;

/*!
 * \brief This constructs a PeriodicAOFactory object
 *
 * \tparam Tile the type of TA Tensor
 * \tparam Policy the type of TA policy
 * \param kv KeyVal object
 * \return the shared pointer of PeriodicAOFactory object
 */
template <typename Tile, typename Policy>
std::shared_ptr<PeriodicAOFactory<Tile, Policy>> construct_periodic_ao_factory(
    const KeyVal &kv) {
  std::shared_ptr<PeriodicAOFactory<Tile, Policy>> pao_factory;

  if (kv.exists_class("wfn_world:periodic_ao_factory")) {
    pao_factory = kv.class_ptr<PeriodicAOFactory<Tile, Policy>>(
        "wfn_world:periodic_ao_factory");
  } else {
    pao_factory = std::make_shared<PeriodicAOFactory<Tile, Policy>>(kv);
    std::shared_ptr<DescribedClass> pao_factory_base = pao_factory;
    KeyVal &kv_nonconst = const_cast<KeyVal &>(kv);
    kv_nonconst.keyval("wfn_world")
        .assign("periodic_ao_factory", pao_factory_base);
  }
  return pao_factory;
}

namespace detail {

/*!
 * \brief This shifts the origin of a Basis object
 *
 * \note All functions in the basis will be shifted
 * \param basis the original Basis object
 * \param shift the 3D vector of the shift
 * \return the shared pointer of shifted Basis object
 */
std::shared_ptr<Basis> shift_basis_origin(Basis &basis, Vector3d shift);

/*!
 * \brief This shifts the origin of a Basis object by multiple vectors,
 * and returns a compound Basis that combines all shifted bases
 *
 * \param basis the original Basis object
 * \param shift_base the base position where all shifting vectors start
 * \param nshift the range of included lattices
 * \param dcell the direct unit cell params
 * \return the shared pointer of the compound Basis object
 */
std::shared_ptr<Basis> shift_basis_origin(Basis &basis, Vector3d shift_base,
                                          Vector3i nshift, Vector3d dcell);

}  // namespace detail

template <typename Tile, typename Policy>
class PeriodicAOFactory : public AOFactory<TA::TensorD, Policy> {
 public:
  using TArray = TA::DistArray<Tile, Policy>;
  using Op = std::function<Tile(TA::TensorZ &&)>;

  PeriodicAOFactory() = default;
  PeriodicAOFactory(PeriodicAOFactory &&) = default;
  PeriodicAOFactory &operator=(PeriodicAOFactory &&) = default;

  /*!
   * \brief KeyVal constructor for PeriodicAOFactory
   * \param kv the KeyVal object
   */
  PeriodicAOFactory(const KeyVal &kv) : AOFactory<TA::TensorD, Policy>(kv) {
    std::string prefix = "";
    if (kv.exists("wfn_world") || kv.exists_class("wfn_world"))
      prefix = "wfn_world:";

    // Molecule was already created at this path, bypass registry and construct
    // UnitCell
    unitcell_ = kv.class_ptr<UnitCell>(prefix + "molecule", true);
    dcell_ = unitcell_->dcell();

    R_max_ =
        decltype(R_max_)(kv.value<std::vector<int>>(prefix + "rmax").data());
    RD_max_ =
        decltype(RD_max_)(kv.value<std::vector<int>>(prefix + "rdmax").data());
    RJ_max_ =
        decltype(RJ_max_)(kv.value<std::vector<int>>(prefix + "rjmax").data());

    using ::mpqc::lcao::detail::direct_ord_idx;
    R_size_ = 1 + direct_ord_idx(R_max_(0), R_max_(1), R_max_(2), R_max_);
    RJ_size_ = 1 + direct_ord_idx(RJ_max_(0), RJ_max_(1), RJ_max_(2), RJ_max_);
    RD_size_ = 1 + direct_ord_idx(RD_max_(0), RD_max_(1), RD_max_(2), RD_max_);

    set_oper(Tile());

    print_detail_ = kv.value<bool>("print_detail", false);
  }

  /// set oper based on Tile type
  template <typename T = Tile>
  void set_oper(typename std::enable_if<std::is_same<T, TA::TensorZ>::value,
                                        T>::type &&t) {
    op_ = TA::Noop<TA::TensorZ, true>();
  }

  ~PeriodicAOFactory() noexcept = default;

  /// wrapper to compute function
  TArray compute(const std::wstring &);

  /*!
   * \brief This computes integral by Formula
   *
   * This function will look into registry first
   * if Formula computed, it will return it from registry
   * if not, it will compute it
   * \param formula the desired Formula type
   * \return the TA::DistArray object
   */
  TArray compute(const Formula &formula);

  /// @return a shared ptr to the UnitCell object
  std::shared_ptr<UnitCell> unitcell() const { return unitcell_; }

  /// @return the range of expansion of Bloch Gaussians in AO Gaussians
  Vector3i R_max() { return R_max_; }

  /// @return the range of Coulomb operation
  Vector3i RJ_max() { return RJ_max_; }

  /// @return the range of density representation
  Vector3i RD_max() { return RD_max_; }

  /// @return the cardinal number of lattices included in Bloch Gaussian
  /// expansion
  int64_t R_size() { return R_size_; }

  /// @return the cardinal number of lattices included in Coulomb operation
  int64_t RJ_size() { return RJ_size_; }

  /// @return the cardinal number of lattices included in density representation
  int64_t RD_size() { return RD_size_; }

  /// @return UnitCell object
  UnitCell &unitcell() { return *unitcell_; }

  /*!
   * \brief This sets the density for coulomb and exchange computations
   * in PeriodicAOFactory
   *
   * \param D is the density feeded to PeriodicAOFactory
   */
  void set_density(TArray D) { D_ = D; }

  /// @return density matrix
  TArray get_density() { return D_; }

  /*!
   * \brief This computes sparse complex array
   *
   * \param world MADNESS world
   * \param engine a utility::TSPool object
   * that is initialized with Operator and bases
   * \param bases std::array of Basis
   * \param p_screen Screener
   * \return the integral sparse array
   */
  template <typename U = Policy>
  TA::DistArray<
      Tile, typename std::enable_if<std::is_same<U, TA::SparsePolicy>::value,
                                    TA::SparsePolicy>::type>
  compute_integrals(madness::World &world, ShrPool<libint2::Engine> &engine,
                    BasisVector const &bases,
                    std::shared_ptr<Screener> p_screen =
                        std::make_shared<Screener>(Screener{})) {
    detail::integral_engine_precision = 0.0;
    auto result = sparse_complex_integrals(world, engine, bases, p_screen, op_);
    return result;
  }

 private:
  /// parse one body formula and set engine_pool and basis array for periodic
  /// system
  void parse_one_body_periodic(
      const Formula &formula,
      std::shared_ptr<utility::TSPool<libint2::Engine>> &engine_pool,
      BasisVector &bases, const Molecule &shifted_mol);

  /// parse two body formula and set engine_pool and basis array for periodic
  /// system
  void parse_two_body_periodic(
      const Formula &formula,
      std::shared_ptr<utility::TSPool<libint2::Engine>> &engine_pool,
      BasisVector &bases, Vector3d shift_coul, bool if_coulomb,
      std::shared_ptr<Screener> &p_screener);

  /*!
   * \brief Construct sparse complex integral tensors in parallel.
   *
   * \param shr_pool should be a std::shared_ptr to an IntegralTSPool
   * \param bases should be a std::array of Basis, which will be copied.
   * \param op needs to be a function or functor that takes a TA::TensorZ && and
   * returns any valid tile type. Op is copied so it can be moved.
   * ```
   * auto t = [](TA::TensorZ &&ten){return std::move(ten);};
   * ```
   *
   * \param screen should be a std::shared_ptr to a Screener.
   */
  template <typename E>
  TA::DistArray<Tile, TA::SparsePolicy> sparse_complex_integrals(
      madness::World &world, ShrPool<E> shr_pool, BasisVector const &bases,
      std::shared_ptr<Screener> screen = std::make_shared<Screener>(Screener{}),
      Op op = TA::Noop<TA::TensorZ, true>());

  std::shared_ptr<UnitCell> unitcell_;  ///> UnitCell private member

  Op op_;

  TArray D_;  ///> Density

  Vector3i R_max_ = {
      0, 0, 0};  ///> range of expansion of Bloch Gaussians in AO Gaussians
  Vector3i RJ_max_ = {0, 0, 0};       ///> range of Coulomb operation
  Vector3i RD_max_ = {0, 0, 0};       ///> range of density representation
  Vector3d dcell_ = {0.0, 0.0, 0.0};  ///> direct unit cell params (in a.u.)

  int64_t
      R_size_;  ///> cardinal # of lattices included in Bloch Gaussian expansion
  int64_t RJ_size_;  ///> cardinal # of lattices included in Coulomb operation
  int64_t
      RD_size_;  ///> cardinal # of lattices included in density representation

  bool print_detail_;  ///> if true, print a lot more details
  FormulaRegistry<TArray> ao_formula_registry_;
};

template <typename Tile, typename Policy>
typename PeriodicAOFactory<Tile, Policy>::TArray
PeriodicAOFactory<Tile, Policy>::compute(const std::wstring &formula_string) {
  auto formula = Formula(formula_string);
  return compute(formula);
}

template <typename Tile, typename Policy>
typename PeriodicAOFactory<Tile, Policy>::TArray
PeriodicAOFactory<Tile, Policy>::compute(const Formula &formula) {
  TArray result;
  BasisVector bs_array;
  std::shared_ptr<utility::TSPool<libint2::Engine>> engine_pool;
  double size = 0.0;

  // retrieve the integral if it is in registry
  auto iter = ao_formula_registry_.find(formula);

  if (iter != ao_formula_registry_.end()) {
    result = iter->second;
    utility::print_par(this->world(), "Retrieved Periodic AO Integral: ",
                       utility::to_string(formula.string()));
    double size = mpqc::detail::array_size(result);
    utility::print_par(this->world(), " Size: ", size, " GB\n");
    return result;
  }

  if (formula.oper().is_fock()) {
    utility::print_par(this->world(),
                       "\nComputing Fock Integral for Periodic System: ",
                       utility::to_string(formula.string()), "\n");

    auto v_formula = formula;
    v_formula.set_operator_type(Operator::Type::Nuclear);

    auto t_formula = formula;
    t_formula.set_operator_type(Operator::Type::Kinetic);

    auto j_formula = formula;
    j_formula.set_operator_type(Operator::Type::J);

    auto k_formula = formula;
    k_formula.set_operator_type(Operator::Type::K);

    auto v = this->compute(v_formula);
    auto t = this->compute(t_formula);
    auto j = this->compute(j_formula);
    auto k = this->compute(k_formula);

    auto time0 = mpqc::now(this->world(), false);
    result("rho, sigma") = v("rho, sigma") + t("rho, sigma");
    result("rho, sigma") += 2.0 * j("rho, sigma") - k("rho, sigma");
    auto time1 = mpqc::now(this->world(), false);
    auto time = mpqc::duration_in_s(time0, time1);

    size = mpqc::detail::array_size(result);
    utility::print_par(this->world(), " Size: ", size, " GB");
    utility::print_par(this->world(), " Time: ", time, " s\n");

    ao_formula_registry_.insert(formula, result);

  } else if (formula.rank() == 2) {
    utility::print_par(this->world(),
                       "\nComputing One Body Integral for Periodic System: ",
                       utility::to_string(formula.string()), "\n");

    auto time0 = mpqc::now(this->world(), false);
    if (formula.oper().type() == Operator::Type::Kinetic ||
        formula.oper().type() == Operator::Type::Overlap) {
      parse_one_body_periodic(formula, engine_pool, bs_array, *unitcell_);
      result = compute_integrals(this->world(), engine_pool, bs_array);
    } else if (formula.oper().type() == Operator::Type::Nuclear) {
      for (auto RJ = 0; RJ < RJ_size_; ++RJ) {
        using ::mpqc::lcao::detail::direct_vector;
        using ::mpqc::lcao::detail::shift_mol_origin;
        auto shift_mol = direct_vector(RJ, RJ_max_, dcell_);
        auto shifted_mol = shift_mol_origin(*unitcell_, shift_mol);
        parse_one_body_periodic(formula, engine_pool, bs_array, *shifted_mol);
        if (RJ == 0)
          result = compute_integrals(this->world(), engine_pool, bs_array);
        else
          result("mu, nu") +=
              compute_integrals(this->world(), engine_pool, bs_array)("mu, nu");
      }
    } else
      throw std::runtime_error("Rank-2 operator type not supported");
    auto time1 = mpqc::now(this->world(), false);
    auto time = mpqc::duration_in_s(time0, time1);

    size = mpqc::detail::array_size(result);
    utility::print_par(this->world(), " Size: ", size, " GB");
    utility::print_par(this->world(), " Time: ", time, " s\n");

    ao_formula_registry_.insert(formula, result);

  } else if (formula.rank() == 4) {
    auto time_4idx = 0.0;
    auto time_contr = 0.0;
    auto time = 0.0;

    std::shared_ptr<Screener> p_screener =
        std::make_shared<Screener>(Screener{});

    if (print_detail_) {
      utility::print_par(this->world(),
                         "\nComputing Two Body Integral for Periodic System: ",
                         utility::to_string(formula.string()), "\n");
    }

    if (formula.oper().type() == Operator::Type::J) {
      auto time_j0 = mpqc::now(this->world(), false);

      auto j_formula = formula;
      j_formula.set_operator_type(Operator::Type::Coulomb);

      for (auto RJ = 0; RJ < RJ_size_; ++RJ) {
        using ::mpqc::lcao::detail::direct_vector;
        auto vec_RJ = direct_vector(RJ, RJ_max_, dcell_);
        parse_two_body_periodic(j_formula, engine_pool, bs_array, vec_RJ, true,
                                p_screener);

        auto time_g0 = mpqc::now(this->world(), false);
        auto J =
            compute_integrals(this->world(), engine_pool, bs_array, p_screener);
        auto time_g1 = mpqc::now(this->world(), false);

        time_4idx += mpqc::duration_in_s(time_g0, time_g1);

        auto time_contr0 = mpqc::now(this->world(), false);
        if (RJ == 0)
          result("mu, nu") = J("mu, nu, lambda, rho") * D_("lambda, rho");
        else
          result("mu, nu") += J("mu, nu, lambda, rho") * D_("lambda, rho");
        auto time_contr1 = mpqc::now(this->world(), false);
        time_contr += mpqc::duration_in_s(time_contr0, time_contr1);
      }
      auto time_j1 = mpqc::now(this->world(), false);
      time = mpqc::duration_in_s(time_j0, time_j1);

    } else if (formula.oper().type() == Operator::Type::K) {
      auto time_k0 = mpqc::now(this->world(), false);

      auto k_formula = formula;
      k_formula.set_operator_type(Operator::Type::Coulomb);
      for (auto RJ = 0; RJ < RJ_size_; ++RJ) {
        using ::mpqc::lcao::detail::direct_vector;
        auto vec_RJ = direct_vector(RJ, RJ_max_, dcell_);
        parse_two_body_periodic(k_formula, engine_pool, bs_array, vec_RJ, false,
                                p_screener);

        auto time_g0 = mpqc::now(this->world(), false);
        auto K =
            compute_integrals(this->world(), engine_pool, bs_array, p_screener);
        auto time_g1 = mpqc::now(this->world(), false);
        time_4idx += mpqc::duration_in_s(time_g0, time_g1);

        auto time_contr0 = mpqc::now(this->world(), false);
        if (RJ == 0)
          result("mu, nu") = K("mu, lambda, nu, rho") * D_("lambda, rho");
        else
          result("mu, nu") += K("mu, lambda, nu, rho") * D_("lambda, rho");
        auto time_contr1 = mpqc::now(this->world(), false);
        time_contr += mpqc::duration_in_s(time_contr0, time_contr1);
      }
      auto time_k1 = mpqc::now(this->world(), false);
      time = mpqc::duration_in_s(time_k0, time_k1);

    } else
      throw std::runtime_error("Rank-4 operator type not supported");

    if (print_detail_) {
      size = mpqc::detail::array_size(result);
      utility::print_par(this->world(), " Size: ", size, " GB\n");
      utility::print_par(this->world(), " \t4-index g tensor time: ", time_4idx,
                         " s\n");
      utility::print_par(this->world(), " \tg*D contraction time: ", time_contr,
                         " s\n");
      utility::print_par(this->world(), " \ttotal time: ", time, " s\n");
    }

  } else
    throw std::runtime_error("Operator rank not supported");

  return result;
}

template <typename Tile, typename Policy>
void PeriodicAOFactory<Tile, Policy>::parse_one_body_periodic(
    const Formula &formula,
    std::shared_ptr<utility::TSPool<libint2::Engine>> &engine_pool,
    BasisVector &bases, const Molecule &shifted_mol) {
  auto bra_indices = formula.bra_indices();
  auto ket_indices = formula.ket_indices();

  TA_ASSERT(bra_indices.size() == 1);
  TA_ASSERT(ket_indices.size() == 1);

  auto bra_index = bra_indices[0];
  auto ket_index = ket_indices[0];

  TA_ASSERT(bra_index.is_ao());
  TA_ASSERT(ket_index.is_ao());

  auto bra_basis = detail::index_to_basis(this->orbital_basis_registry(),bra_index);
  auto ket_basis = detail::index_to_basis(this->orbital_basis_registry(),ket_index);

  TA_ASSERT(bra_basis != nullptr);
  TA_ASSERT(ket_basis != nullptr);

  // Form a compound ket basis by shifting origins from -Rmax to Rmax
  Vector3d zero_shift_base(0.0, 0.0, 0.0);
  ket_basis =
      detail::shift_basis_origin(*ket_basis, zero_shift_base, R_max_, dcell_);

  bases = BasisVector{{*bra_basis, *ket_basis}};

  auto oper_type = formula.oper().type();
  engine_pool = make_engine_pool(
      detail::to_libint2_operator(oper_type),
      utility::make_array_of_refs(*bra_basis, *ket_basis), libint2::BraKet::x_x,
      detail::to_libint2_operator_params(oper_type, shifted_mol));
}

template <typename Tile, typename Policy>
void PeriodicAOFactory<Tile, Policy>::parse_two_body_periodic(
    const Formula &formula,
    std::shared_ptr<utility::TSPool<libint2::Engine>> &engine_pool,
    BasisVector &bases, Vector3d shift_coul, bool if_coulomb,
    std::shared_ptr<Screener> &p_screener) {
  auto bra_indices = formula.bra_indices();
  auto ket_indices = formula.ket_indices();

  TA_ASSERT(bra_indices.size() == 2);
  TA_ASSERT(ket_indices.size() == 2);

  auto bra_index0 = bra_indices[0];
  auto bra_index1 = bra_indices[0];
  auto ket_index0 = ket_indices[0];
  auto ket_index1 = ket_indices[0];

  TA_ASSERT(bra_index0.is_ao());
  TA_ASSERT(bra_index1.is_ao());
  TA_ASSERT(ket_index0.is_ao());
  TA_ASSERT(ket_index1.is_ao());

  auto bra_basis0 = detail::index_to_basis(this->orbital_basis_registry(), bra_index0);
  auto bra_basis1 = detail::index_to_basis(this->orbital_basis_registry(), bra_index1);
  auto ket_basis0 = detail::index_to_basis(this->orbital_basis_registry(), ket_index0);
  auto ket_basis1 = detail::index_to_basis(this->orbital_basis_registry(), ket_index1);

  TA_ASSERT(bra_basis0 != nullptr);
  TA_ASSERT(bra_basis1 != nullptr);
  TA_ASSERT(ket_basis0 != nullptr);
  TA_ASSERT(ket_basis1 != nullptr);

  // Form a compound index basis
  Vector3d zero_shift_base(0.0, 0.0, 0.0);
  if (if_coulomb) {
    bra_basis1 = detail::shift_basis_origin(*bra_basis1, zero_shift_base,
                                            R_max_, dcell_);
    ket_basis0 = detail::shift_basis_origin(*ket_basis0, shift_coul);
  } else {
    bra_basis1 = detail::shift_basis_origin(*bra_basis1, shift_coul);
    ket_basis0 = detail::shift_basis_origin(*ket_basis0, zero_shift_base,
                                            R_max_, dcell_);
  }
  ket_basis1 =
      detail::shift_basis_origin(*ket_basis1, shift_coul, RD_max_, dcell_);

  if (formula.notation() == Formula::Notation::Chemical)
    bases = BasisVector{{*bra_basis0, *bra_basis1, *ket_basis0, *ket_basis1}};
  else
    throw "Physical notation not supported!";
  //    bases = BasisVector{{*bra_basis0, *ket_basis0, *bra_basis1,
  //    *ket_basis1}};

  auto oper_type = formula.oper().type();
  engine_pool = make_engine_pool(
      detail::to_libint2_operator(oper_type),
      utility::make_array_of_refs(bases[0], bases[1], bases[2], bases[3]),
      libint2::BraKet::xx_xx,
      detail::to_libint2_operator_params(oper_type, *unitcell_));

  p_screener = detail::make_screener(this->world(),engine_pool, bases,this->screen_, this->screen_threshold_);
}

template <typename Tile, typename Policy>
template <typename E>
TA::DistArray<Tile, TA::SparsePolicy>
PeriodicAOFactory<Tile, Policy>::sparse_complex_integrals(
    madness::World &world, ShrPool<E> shr_pool, BasisVector const &bases,
    std::shared_ptr<Screener> screen, Op op) {
  auto time0 = mpqc::now(this->world(), true);
  // Build the Trange and Shape Tensor
  auto trange = detail::create_trange(bases);
  const auto tvolume = trange.tiles_range().volume();
  std::vector<std::pair<unsigned long, Tile>> tiles(tvolume);
  TA::TensorF tile_norms(trange.tiles_range(), 0.0);

  // Copy the Bases for the Integral Builder
  auto shr_bases = std::make_shared<BasisVector>(bases);

  // Make a pointer to an Integral builder.  Doing this because we want to use
  // it in Tasks.
  auto builder_ptr = std::make_shared<IntegralBuilder<Tile, E>>(
      std::move(shr_pool), std::move(shr_bases), std::move(screen),
      std::move(op));

  auto task_f = [=](int64_t ord, detail::IdxVec idx, TA::Range rng,
                    TA::TensorF *tile_norms_ptr, Tile *out_tile) {

    // This is why builder was made into a shared_ptr.
    auto &builder = *builder_ptr;
    auto ta_tile = builder.integrals(idx, std::move(rng));

    const auto tile_volume = ta_tile.range().volume();
    const auto tile_norm = ta_tile.norm();

    // Keep tile if it was significant.
    bool save_norm =
        tile_norm >= tile_volume * TA::SparseShape<float>::threshold();
    if (save_norm) {
      *out_tile = builder.op(std::move(ta_tile));

      auto &norms = *tile_norms_ptr;
      norms[ord] = tile_norm;
    }
  };

  auto pmap = TA::SparsePolicy::default_pmap(world, tvolume);

  auto time_f0 = mpqc::now(this->world(), true);
  for (auto const ord : *pmap) {
    tiles[ord].first = ord;
    detail::IdxVec idx = trange.tiles_range().idx(ord);
    world.taskq.add(task_f, ord, idx, trange.make_tile_range(ord), &tile_norms,
                    &tiles[ord].second);
  }
  world.gop.fence();
  auto time_f1 = mpqc::now(this->world(), true);
  auto time_f = mpqc::duration_in_s(time_f0, time_f1);

  TA::SparseShape<float> shape(world, tile_norms, trange);
  TA::DistArray<Tile, TA::SparsePolicy> out(world, trange, shape, pmap);

  detail::set_array(tiles, out);
  out.truncate();
  auto time1 = mpqc::now(this->world(), true);
  auto time = mpqc::duration_in_s(time0, time1);

  if (print_detail_) {
    utility::print_par(this->world(), " \tsum of task_f time: ", time_f, " s\n",
                       " \ttotal compute time: ", time, " s\n\n");
  }

  return out;
}

/// Make PeriodicAOFactory printable
template <typename Tile, typename Policy>
std::ostream &operator<<(std::ostream &os,
                         PeriodicAOFactory<Tile, Policy> &pao) {
  os << "\nPeriodicAOFactory computational parameters:" << std::endl;
  os << "\tR_max (range of expansion of Bloch Gaussians in AO Gaussians): ["
     << pao.R_max().transpose() << "]" << std::endl;
  os << "\tRj_max (range of Coulomb operation): [" << pao.RJ_max().transpose()
     << "]" << std::endl;
  os << "\tRd_max (Range of density representation): ["
     << pao.RD_max().transpose() << "]" << std::endl;
  return os;
}

extern template class PeriodicAOFactory<TA::TensorZ, TA::SparsePolicy>;

}  // namespace gaussian
}  // namespace lcao
}  // namespace mpqc
#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_LCAO_FACTORY_PERIODIC_AO_FACTORY_H_
