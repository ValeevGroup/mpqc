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
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vectord;

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

template <typename Tile, typename Policy>
using PeriodicAOFactoryBase =
    Factory<TA::DistArray<Tile, Policy>, TA::DistArray<Tile, Policy>>;
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
class PeriodicAOFactory : public PeriodicAOFactoryBase<Tile, Policy> {
 public:
  using TArray = TA::DistArray<Tile, Policy>;
  using DirectTArray = DirectArray<Tile, Policy>;
  using Op = std::function<Tile(TA::TensorD &&)>;

 public:
  PeriodicAOFactory() = default;
  PeriodicAOFactory(PeriodicAOFactory &&) = default;
  PeriodicAOFactory &operator=(PeriodicAOFactory &&) = default;

  /*!
   * \brief KeyVal constructor for PeriodicAOFactory
   * \param kv the KeyVal object
   */
  PeriodicAOFactory(const KeyVal &kv)
      : PeriodicAOFactoryBase<Tile, Policy>(kv) {
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

    auto default_precision = std::numeric_limits<double>::epsilon();
    precision_ = kv.value<double>(prefix + "precision", default_precision);
    detail::integral_engine_precision = precision_;
    ExEnv::out0() << indent << "Precision = " << precision_ << "\n";

    screen_ = kv.value<std::string>(prefix + "screen", "schwarz");
    screen_threshold_ = kv.value<double>(prefix + "threshold", 1.0e-10);

    //    auto convert_op = [](TA::TensorD &&arg) -> TA::TensorZ {
    //      return TA::TensorZ(arg.range(), arg.data());
    //    };
    //    op_ = convert_op;
    detail::set_oper(op_);

    print_detail_ = kv.value<bool>("print_detail", false);

    auto orbital_space_registry =
        std::make_shared<OrbitalSpaceRegistry<TArray>>();

    this->set_orbital_registry(orbital_space_registry);

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
  TArray compute(const Formula &formula) override;

  /// wrapper to compute_direct
  TArray compute_direct(const std::wstring &);

  /// This computes integral direct by Formula
  TArray compute_direct(const Formula &formula) override;

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
  TArray compute_integrals(madness::World &world,
                           ShrPool<libint2::Engine> &engine,
                           BasisVector const &bases,
                           std::shared_ptr<Screener> p_screen =
                               std::make_shared<Screener>(Screener{})) {
    auto result = sparse_complex_integrals(world, engine, bases, p_screen, op_);
    return result;
  }

  /// This computes sparse complex array using integral direct
  DirectTArray compute_direct_integrals(
      madness::World &world, ShrPool<libint2::Engine> &engine,
      BasisVector const &bases, std::shared_ptr<Screener> p_screen =
                                    std::make_shared<Screener>(Screener{})) {
    auto result =
        direct_sparse_complex_integrals(world, engine, bases, p_screen, op_);
    return result;
  }

 public:
  /// @return a shared ptr to the UnitCell object
  std::shared_ptr<UnitCell> unitcell() const { return unitcell_; }

  /// @return screen method
  const std::string &screen() const { return screen_; }

  /// @return screen threshold
  double screen_threshold() const { return screen_threshold_; }

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

  /// This sets the density for coulomb and exchange computations
  void set_density(TArray D) { D_ = D; }

  /// @return density matrix
  TArray get_density() { return D_; }

 private:
  /// compute integrals that has one dimension
  TArray compute1(const Formula &formula);

  /// compute integrals that has two dimension for periodic systems
  TArray compute2(const Formula &formula);

  /// compute integrals that has two dimension for periodic systems
  TArray compute3(const Formula &formula);

  /// compute integrals that has four dimension for periodic systems
  TArray compute4(const Formula &formula);

  /// compute integrals that has two dimension for periodic systems
  //  TArray compute_direct2(const Formula &formula);

  /// compute integrals that has two dimension for periodic systems
  TArray compute_direct4(const Formula &formula);

  /// parse 1-body 1-center formula written as <X|U> and set engine_pool and
  /// basis array
  void parse_one_body_one_center(
      const Formula &formula,
      std::shared_ptr<utility::TSPool<libint2::Engine>> &engine_pool,
      BasisVector &bases);

  /// parse 1-body 2-center formula and set engine_pool and basis array for
  /// periodic
  /// system
  void parse_one_body_two_center_periodic(
      const Formula &formula,
      std::shared_ptr<utility::TSPool<libint2::Engine>> &engine_pool,
      BasisVector &bases, const Molecule &shifted_mol);

  /// parse 2-body 2-center formula and set engine_pool and basis array for
  /// periodic
  /// system. The ket basis is shifted by \c shift
  void parse_two_body_two_center_periodic(
      const Formula &formula,
      std::shared_ptr<utility::TSPool<libint2::Engine>> &engine_pool,
      BasisVector &bases, const Vector3d &shift);

  /// parse 2-body 3-center formula and set engine_pool and basis array for
  /// periodic
  /// system. The ket basis is shifted by \c shift
  void parse_two_body_three_center_periodic(
      const Formula &formula,
      std::shared_ptr<utility::TSPool<libint2::Engine>> &engine_pool,
      BasisVector &bases, const Vector3d &shift,
      std::shared_ptr<Screener> &p_screener);

  /// parse 2-body 4-center formula and set engine_pool and basis array for
  /// periodic
  /// system
  void parse_two_body_four_center_periodic(
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
  TArray sparse_complex_integrals(madness::World &world, ShrPool<E> shr_pool,
                                  BasisVector const &bases,
                                  std::shared_ptr<Screener> screen, Op op);

  /// This constructs direct sparse complex integral tensors in parallel
  template <typename E>
  DirectTArray direct_sparse_complex_integrals(madness::World &world,
                                               ShrPool<E> shr_pool,
                                               BasisVector const &bases,
                                               std::shared_ptr<Screener> screen,
                                               Op op);

 private:
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
  std::string screen_;
  double precision_;
  double screen_threshold_;
  std::vector<DirectTArray> gj_;
  std::vector<DirectTArray> gk_;
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

  // retrieve the integral if it is in registry
  auto iter = this->registry_.find(formula);

  if (iter != this->registry_.end()) {
    result = iter->second;
    utility::print_par(this->world(), "Retrieved Periodic AO Integral: ",
                       utility::to_string(formula.string()));
    double size = mpqc::detail::array_size(result);
    utility::print_par(this->world(), " Size: ", size, " GB\n");
    return result;
  }

  if (formula.rank() == 1) {
    result = compute1(formula);
    this->registry_.insert(formula, result);
  } else if (formula.rank() == 2) {
    result = compute2(formula);
    this->registry_.insert(formula, result);
  } else if (formula.rank() == 3) {
    result = compute3(formula);
    this->registry_.insert(formula, result);
  } else if (formula.rank() == 4) {
    result = compute4(formula);
  } else
    throw std::runtime_error("Operator rank not supported");

  return result;
}

template <typename Tile, typename Policy>
typename PeriodicAOFactory<Tile, Policy>::TArray
PeriodicAOFactory<Tile, Policy>::compute1(const Formula &formula) {
  BasisVector bs_array;
  TArray result;

  std::shared_ptr<utility::TSPool<libint2::Engine>> engine_pool;
  auto &world = this->world();

  double size = 0.0;

  ExEnv::out0() << "\nComputing One Center Integral for Periodic System: "
                << utility::to_string(formula.string()) << std::endl;
  exit(19);

  auto basis =
      detail::index_to_basis(*this->basis_registry(), formula.bra_indices()[0]);
  auto emptybasis = Basis();
  auto oper_type = formula.oper().type();

  engine_pool = make_engine_pool(
      detail::to_libint2_operator(oper_type),
      utility::make_array_of_refs(*basis, emptybasis), libint2::BraKet::x_x,
      detail::to_libint2_operator_params(oper_type, *unitcell_));
  bs_array = BasisVector{{*basis, emptybasis}};

  result = compute_integrals(world, engine_pool, bs_array);

  return result;
}

template <typename Tile, typename Policy>
typename PeriodicAOFactory<Tile, Policy>::TArray
PeriodicAOFactory<Tile, Policy>::compute2(const Formula &formula) {
  BasisVector bs_array;
  TArray result;
  std::shared_ptr<utility::TSPool<libint2::Engine>> engine_pool;
  auto &world = this->world();
  double size = 0.0;

  ExEnv::out0() << "\nComputing Two Center Integral for Periodic System: "
                << utility::to_string(formula.string()) << std::endl;

  auto time0 = mpqc::now(world, this->accurate_time());

  if (formula.oper().type() == Operator::Type::Identity) {
    // Identity matrix

    auto bra_index = formula.bra_indices()[0];
    auto ket_index = formula.ket_indices()[0];
    auto bra_basis = this->basis_registry()->retrieve(bra_index);
    auto ket_basis = this->basis_registry()->retrieve(ket_index);
    auto bra_tr = bra_basis->create_trange1();
    auto ket_tr = ket_basis->create_trange1();
    // create diagonal array
    result = array_ops::create_diagonal_array_from_eigen<Tile, Policy>(
        world, bra_tr, ket_tr, 1.0);
    result.truncate();

  } else if (formula.bra_indices()[0] == OrbitalIndex(L"U") ||
             formula.ket_indices()[0] == OrbitalIndex(L"U")) {
    // compute one center integrals written as <X|U> or <U|X>
    parse_one_body_one_center(formula, engine_pool, bs_array);
    auto result_2D = compute_integrals(world, engine_pool, bs_array);
    result = result_2D;

  } else if (formula.oper().type() == Operator::Type::Kinetic ||
             formula.oper().type() == Operator::Type::Overlap) {
    parse_one_body_two_center_periodic(formula, engine_pool, bs_array,
                                       *unitcell_);
    result = compute_integrals(world, engine_pool, bs_array);
  } else if (formula.oper().type() == Operator::Type::Nuclear) {
    for (auto RJ = 0; RJ < RJ_size_; ++RJ) {
      using ::mpqc::lcao::detail::direct_vector;
      using ::mpqc::lcao::detail::shift_mol_origin;
      auto shift_mol = direct_vector(RJ, RJ_max_, dcell_);
      auto shifted_mol = shift_mol_origin(*unitcell_, shift_mol);
      parse_one_body_two_center_periodic(formula, engine_pool, bs_array,
                                         *shifted_mol);
      if (RJ == 0)
        result = compute_integrals(world, engine_pool, bs_array);
      else
        result("mu, nu") +=
            compute_integrals(world, engine_pool, bs_array)("mu, nu");
    }
  } else if (formula.oper().type() == Operator::Type::Coulomb) {
    for (auto RJ = 0; RJ < RJ_size_; ++RJ) {
      using ::mpqc::lcao::detail::direct_vector;
      auto vec_RJ = direct_vector(RJ, RJ_max_, dcell_);
      parse_two_body_two_center_periodic(formula, engine_pool, bs_array,
                                         vec_RJ);
      if (RJ == 0)
        result = compute_integrals(world, engine_pool, bs_array);
      else
        result("mu, nu") +=
            compute_integrals(world, engine_pool, bs_array)("mu, nu");
    }
  } else
    throw std::runtime_error("Rank-2 operator type not supported");
  auto time1 = mpqc::now(world, this->accurate_time());
  auto time = mpqc::duration_in_s(time0, time1);

  size = mpqc::detail::array_size(result);
  utility::print_par(world, " Size: ", size, " GB");
  utility::print_par(world, " Time: ", time, " s\n");

  return result;
}

template <typename Tile, typename Policy>
typename PeriodicAOFactory<Tile, Policy>::TArray
PeriodicAOFactory<Tile, Policy>::compute3(const Formula &formula) {
  TArray result;
  BasisVector bs_array;
  std::shared_ptr<utility::TSPool<libint2::Engine>> engine_pool;
  auto &world = this->world();
  std::shared_ptr<Screener> p_screener = std::make_shared<Screener>(Screener{});

  double size = 0.0;

  ExEnv::out0() << "\nComputing Two Center Integral for Periodic System: "
                << utility::to_string(formula.string()) << std::endl;

  auto time0 = mpqc::now(world, this->accurate_time());
  if (formula.oper().type() == Operator::Type::Coulomb) {
    for (auto RJ = 0; RJ < RJ_size_; ++RJ) {
      using ::mpqc::lcao::detail::direct_vector;
      auto vec_RJ = direct_vector(RJ, RJ_max_, dcell_);
      parse_two_body_three_center_periodic(formula, engine_pool, bs_array,
                                           vec_RJ, p_screener);
      auto g = compute_integrals(world, engine_pool, bs_array, p_screener);
      if (RJ == 0)
        result("mu, nu, K") = g("mu, nu, K");
      else
        result("mu, nu, K") += g("mu, nu, K");
    }
  } else
    throw std::runtime_error("Rank-3 operator type not supported");
  auto time1 = mpqc::now(world, this->accurate_time());
  auto time = mpqc::duration_in_s(time0, time1);

  size = mpqc::detail::array_size(result);
  utility::print_par(world, " Size: ", size, " GB");
  utility::print_par(world, " Time: ", time, " s\n");

  return result;
}

template <typename Tile, typename Policy>
typename PeriodicAOFactory<Tile, Policy>::TArray
PeriodicAOFactory<Tile, Policy>::compute4(const Formula &formula) {
  BasisVector bs_array;
  TArray result;
  std::shared_ptr<utility::TSPool<libint2::Engine>> engine_pool;
  auto &world = this->world();
  double size = 0.0;

  auto time_4idx = 0.0;
  auto time_contr = 0.0;
  auto time = 0.0;

  std::shared_ptr<Screener> p_screener = std::make_shared<Screener>(Screener{});

  if (print_detail_) {
    utility::print_par(world,
                       "\nComputing Two Body Integral for Periodic System: ",
                       utility::to_string(formula.string()), "\n");
  }

  if (formula.oper().type() == Operator::Type::J) {
    auto time_j0 = mpqc::now(world, this->accurate_time());

    auto j_formula = formula;
    j_formula.set_operator_type(Operator::Type::Coulomb);

    for (auto RJ = 0; RJ < RJ_size_; ++RJ) {
      using ::mpqc::lcao::detail::direct_vector;
      auto vec_RJ = direct_vector(RJ, RJ_max_, dcell_);
      parse_two_body_four_center_periodic(j_formula, engine_pool, bs_array,
                                          vec_RJ, true, p_screener);

      auto time_g0 = mpqc::now(world, this->accurate_time());
      auto g = compute_integrals(world, engine_pool, bs_array, p_screener);
      auto time_g1 = mpqc::now(world, this->accurate_time());

      if (print_detail_) {
        double size = mpqc::detail::array_size(g);
        ExEnv::out0() << " Size of 4-index g tensor:" << size << " GB"
                      << std::endl;
      }
      time_4idx += mpqc::duration_in_s(time_g0, time_g1);

      auto time_contr0 = mpqc::now(world, this->accurate_time());
      if (RJ == 0)
        result("mu, nu") = g("mu, nu, lambda, rho") * D_("lambda, rho");
      else
        result("mu, nu") += g("mu, nu, lambda, rho") * D_("lambda, rho");
      auto time_contr1 = mpqc::now(world, this->accurate_time());
      time_contr += mpqc::duration_in_s(time_contr0, time_contr1);
    }
    auto time_j1 = mpqc::now(world, this->accurate_time());
    time = mpqc::duration_in_s(time_j0, time_j1);

  } else if (formula.oper().type() == Operator::Type::K) {
    auto time_k0 = mpqc::now(world, this->accurate_time());

    auto k_formula = formula;
    k_formula.set_operator_type(Operator::Type::Coulomb);
    for (auto RJ = 0; RJ < RJ_size_; ++RJ) {
      using ::mpqc::lcao::detail::direct_vector;
      auto vec_RJ = direct_vector(RJ, RJ_max_, dcell_);
      parse_two_body_four_center_periodic(k_formula, engine_pool, bs_array,
                                          vec_RJ, false, p_screener);

      auto time_g0 = mpqc::now(world, this->accurate_time());
      auto g = compute_integrals(world, engine_pool, bs_array, p_screener);
      auto time_g1 = mpqc::now(world, this->accurate_time());

      if (print_detail_) {
        double size = mpqc::detail::array_size(g);
        ExEnv::out0() << " Size of 4-index g tensor:" << size << " GB"
                      << std::endl;
      }
      time_4idx += mpqc::duration_in_s(time_g0, time_g1);

      auto time_contr0 = mpqc::now(world, this->accurate_time());
      if (RJ == 0)
        result("mu, nu") = g("mu, lambda, nu, rho") * D_("lambda, rho");
      else
        result("mu, nu") += g("mu, lambda, nu, rho") * D_("lambda, rho");
      auto time_contr1 = mpqc::now(world, this->accurate_time());
      time_contr += mpqc::duration_in_s(time_contr0, time_contr1);
    }
    auto time_k1 = mpqc::now(world, this->accurate_time());
    time = mpqc::duration_in_s(time_k0, time_k1);

  } else
    throw std::runtime_error("Rank-4 operator type not supported");

  if (print_detail_) {
    size = mpqc::detail::array_size(result);
    utility::print_par(world, " Size: ", size, " GB\n");
    utility::print_par(world, " \t4-index g tensor time: ", time_4idx, " s\n");
    utility::print_par(world, " \tg*D contraction time: ", time_contr, " s\n");
    utility::print_par(world, " \ttotal time: ", time, " s\n");
  }

  return result;
}

template <typename Tile, typename Policy>
typename PeriodicAOFactory<Tile, Policy>::TArray
PeriodicAOFactory<Tile, Policy>::compute_direct(
    const std::wstring &formula_string) {
  auto formula = Formula(formula_string);
  return compute_direct(formula);
}

template <typename Tile, typename Policy>
typename PeriodicAOFactory<Tile, Policy>::TArray
PeriodicAOFactory<Tile, Policy>::compute_direct(const Formula &formula) {
  TA_USER_ASSERT(lcao::detail::if_all_ao(formula),
                 "AOFactory only accept AO index!\n");
  TArray result;
  // retrieve the integral if it is in registry
  auto iter = this->registry_.find(formula);

  if (iter != this->registry_.end()) {
    result = iter->second;
    utility::print_par(this->world(), "Retrieved Periodic AO Integral: ",
                       utility::to_string(formula.string()));
    double size = mpqc::detail::array_size(result);
    utility::print_par(this->world(), " Size: ", size, " GB\n");
    return result;
  } else {
    if (formula.rank() == 2) {
      result = compute2(formula);
      this->registry_.insert(formula, result);
    } else if (formula.rank() == 4) {
      result = compute_direct4(formula);
    } else
      throw std::runtime_error("Operator rank not supported");
  }

  return result;
}

template <typename Tile, typename Policy>
typename PeriodicAOFactory<Tile, Policy>::TArray
PeriodicAOFactory<Tile, Policy>::compute_direct4(const Formula &formula) {
  BasisVector bs_array;
  TArray result;
  std::shared_ptr<utility::TSPool<libint2::Engine>> engine_pool;
  auto &world = this->world();
  double size = 0.0;

  auto time_4idx = 0.0;
  auto time_contr = 0.0;
  auto time = 0.0;

  std::shared_ptr<Screener> p_screener = std::make_shared<Screener>(Screener{});

  if (print_detail_) {
    utility::print_par(world,
                       "\nComputing Two Body Integral for Periodic System: ",
                       utility::to_string(formula.string()), "\n");
  }

  if (formula.oper().type() == Operator::Type::J) {
    auto time_j0 = mpqc::now(world, this->accurate_time());

    auto j_formula = formula;
    j_formula.set_operator_type(Operator::Type::Coulomb);

    if (gj_.empty()) {
      gj_ = std::vector<DirectTArray>(RJ_size_, DirectTArray());
    }

    for (auto RJ = 0; RJ < RJ_size_; ++RJ) {
      using ::mpqc::lcao::detail::direct_vector;

      DirectTArray &g = gj_[RJ];
      if (!g.array().is_initialized()) {
        auto vec_RJ = direct_vector(RJ, RJ_max_, dcell_);
        parse_two_body_four_center_periodic(j_formula, engine_pool, bs_array,
                                            vec_RJ, true, p_screener);

        // std::cout << "In compute 4 J\n" <<
        // p_screener->norm_estimate(bs_array) << std::endl;

        auto time_g0 = mpqc::now(world, this->accurate_time());
        g = compute_direct_integrals(world, engine_pool, bs_array, p_screener);
        auto time_g1 = mpqc::now(world, this->accurate_time());

        if (print_detail_) {
          double size = mpqc::detail::array_size(g.array());
          ExEnv::out0() << " Size of 4-index g tensor:" << size << " GB"
                        << std::endl;
        }
        time_4idx += mpqc::duration_in_s(time_g0, time_g1);
      }

      auto time_contr0 = mpqc::now(world, this->accurate_time());
      if (RJ == 0)
        result("mu, nu") = g("mu, nu, lambda, rho") * D_("lambda, rho");
      else
        result("mu, nu") += g("mu, nu, lambda, rho") * D_("lambda, rho");
      auto time_contr1 = mpqc::now(world, this->accurate_time());
      time_contr += mpqc::duration_in_s(time_contr0, time_contr1);
    }
    auto time_j1 = mpqc::now(world, this->accurate_time());
    time = mpqc::duration_in_s(time_j0, time_j1);

  } else if (formula.oper().type() == Operator::Type::K) {
    auto time_k0 = mpqc::now(world, this->accurate_time());

    auto k_formula = formula;
    k_formula.set_operator_type(Operator::Type::Coulomb);

    if (gk_.empty()) {
      gk_ = std::vector<DirectTArray>(RJ_size_, DirectTArray());
    }

    for (auto RJ = 0; RJ < RJ_size_; ++RJ) {
      using ::mpqc::lcao::detail::direct_vector;
      DirectTArray &g = gk_[RJ];
      if (!g.array().is_initialized()) {
        auto vec_RJ = direct_vector(RJ, RJ_max_, dcell_);
        parse_two_body_four_center_periodic(k_formula, engine_pool, bs_array,
                                            vec_RJ, false, p_screener);

        // std::cout << "In compute 4 K\n" <<
        // p_screener->norm_estimate(bs_array) << std::endl;

        auto time_g0 = mpqc::now(world, this->accurate_time());
        g = compute_direct_integrals(world, engine_pool, bs_array, p_screener);
        auto time_g1 = mpqc::now(world, this->accurate_time());

        if (print_detail_) {
          double size = mpqc::detail::array_size(g.array());
          ExEnv::out0() << " Size of 4-index g tensor:" << size << " GB"
                        << std::endl;
        }
        time_4idx += mpqc::duration_in_s(time_g0, time_g1);
      }

      auto time_contr0 = mpqc::now(world, this->accurate_time());
      if (RJ == 0)
        result("mu, nu") = g("mu, lambda, nu, rho") * D_("lambda, rho");
      else
        result("mu, nu") += g("mu, lambda, nu, rho") * D_("lambda, rho");
      auto time_contr1 = mpqc::now(world, this->accurate_time());
      time_contr += mpqc::duration_in_s(time_contr0, time_contr1);
    }
    auto time_k1 = mpqc::now(world, this->accurate_time());
    time = mpqc::duration_in_s(time_k0, time_k1);

  } else
    throw ProgrammingError("Rank-4 operator type not supported", __FILE__,
                           __LINE__);

  if (print_detail_) {
    size = mpqc::detail::array_size(result);
    utility::print_par(world, " Size: ", size, " GB\n");
    utility::print_par(world, " \t4-index g tensor time: ", time_4idx, " s\n");
    utility::print_par(world, " \tg*D contraction time: ", time_contr, " s\n");
    utility::print_par(world, " \ttotal time: ", time, " s\n");
  }

  return result;
}

template <typename Tile, typename Policy>
void PeriodicAOFactory<Tile, Policy>::parse_one_body_one_center(
    const Formula &formula,
    std::shared_ptr<utility::TSPool<libint2::Engine>> &engine_pool,
    BasisVector &bases) {
  auto bra_indices = formula.bra_indices();
  auto ket_indices = formula.ket_indices();

  TA_ASSERT(bra_indices.size() == 1);
  TA_ASSERT(ket_indices.size() == 1);

  auto bra_index = bra_indices[0];
  auto ket_index = ket_indices[0];

  TA_ASSERT(bra_index == OrbitalIndex(L"U") || ket_index == OrbitalIndex(L"U"));

  if (bra_index == OrbitalIndex(L"U"))
    TA_ASSERT(ket_index.is_ao());
  else
    TA_ASSERT(bra_index.is_ao());

  const auto &basis_registry = *this->basis_registry();

  auto bra_basis = detail::index_to_basis(basis_registry, bra_index);
  auto ket_basis = detail::index_to_basis(basis_registry, ket_index);

  TA_ASSERT(bra_basis != nullptr);
  TA_ASSERT(ket_basis != nullptr);

  bases = BasisVector{{*bra_basis, *ket_basis}};

  auto oper_type = formula.oper().type();
  engine_pool = make_engine_pool(
      detail::to_libint2_operator(oper_type),
      utility::make_array_of_refs(*bra_basis, *ket_basis), libint2::BraKet::x_x,
      detail::to_libint2_operator_params(oper_type, *unitcell_));
}

template <typename Tile, typename Policy>
void PeriodicAOFactory<Tile, Policy>::parse_one_body_two_center_periodic(
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

  const auto &basis_registry = *this->basis_registry();

  auto bra_basis = detail::index_to_basis(basis_registry, bra_index);
  auto ket_basis = detail::index_to_basis(basis_registry, ket_index);

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
void PeriodicAOFactory<Tile, Policy>::parse_two_body_two_center_periodic(
    const Formula &formula,
    std::shared_ptr<utility::TSPool<libint2::Engine>> &engine_pool,
    BasisVector &bases, const Vector3d &shift) {
  auto bra_indices = formula.bra_indices();
  auto ket_indices = formula.ket_indices();

  TA_ASSERT(bra_indices.size() == 1);
  TA_ASSERT(ket_indices.size() == 1);

  auto bra_index = bra_indices[0];
  auto ket_index = ket_indices[0];

  TA_ASSERT(bra_index.is_ao());
  TA_ASSERT(ket_index.is_ao());

  const auto &basis_registry = *this->basis_registry();

  auto bra_basis = detail::index_to_basis(basis_registry, bra_index);
  auto ket_basis = detail::index_to_basis(basis_registry, ket_index);

  TA_ASSERT(bra_basis != nullptr);
  TA_ASSERT(ket_basis != nullptr);

  // Shift ket basis
  ket_basis = detail::shift_basis_origin(*ket_basis, shift);

  bases = BasisVector{{*bra_basis, *ket_basis}};

  auto oper_type = formula.oper().type();
  engine_pool = make_engine_pool(
      detail::to_libint2_operator(oper_type),
      utility::make_array_of_refs(*bra_basis, *ket_basis),
      libint2::BraKet::xs_xs,
      detail::to_libint2_operator_params(oper_type, *unitcell_));
}

template <typename Tile, typename Policy>
void PeriodicAOFactory<Tile, Policy>::parse_two_body_three_center_periodic(
    const Formula &formula,
    std::shared_ptr<utility::TSPool<libint2::Engine>> &engine_pool,
    BasisVector &bases, const Vector3d &shift,
    std::shared_ptr<Screener> &p_screener) {
  auto bra_indices = formula.bra_indices();
  auto ket_indices = formula.ket_indices();

  TA_ASSERT(bra_indices.size() == 2);
  TA_ASSERT(ket_indices.size() == 1);

  auto bra_index0 = bra_indices[0];
  auto bra_index1 = bra_indices[1];
  auto ket_index0 = ket_indices[0];

  TA_ASSERT(bra_index0.is_ao());
  TA_ASSERT(bra_index1.is_ao());
  TA_ASSERT(ket_index0.is_ao());

  const auto &basis_registry = *this->basis_registry();

  auto bra_basis0 = detail::index_to_basis(basis_registry, bra_index0);
  auto bra_basis1 = detail::index_to_basis(basis_registry, bra_index1);
  auto ket_basis0 = detail::index_to_basis(basis_registry, ket_index0);

  TA_ASSERT(bra_basis0 != nullptr);
  TA_ASSERT(bra_basis1 != nullptr);
  TA_ASSERT(ket_basis0 != nullptr);

  // Form a compound index bra1 basis
  Vector3d zero_shift_base(0.0, 0.0, 0.0);
  bra_basis1 =
      detail::shift_basis_origin(*bra_basis1, zero_shift_base, R_max_, dcell_);
  // Shift ket basis
  ket_basis0 = detail::shift_basis_origin(*ket_basis0, shift);

  if (formula.notation() == Formula::Notation::Chemical)
    bases = BasisVector{{*bra_basis0, *bra_basis1, *ket_basis0}};
  else
    throw "Physical notation not supported!";

  auto oper_type = formula.oper().type();
  engine_pool = make_engine_pool(
      detail::to_libint2_operator(oper_type),
      utility::make_array_of_refs(*bra_basis0, *bra_basis1, *ket_basis0),
      libint2::BraKet::xx_xs,
      detail::to_libint2_operator_params(oper_type, *unitcell_));

  if (!screen_.empty()) {
    auto screen_engine_pool = make_engine_pool(
        detail::to_libint2_operator(oper_type),
        utility::make_array_of_refs(*bra_basis0, *bra_basis1, *ket_basis0),
        libint2::BraKet::xx_xx,
        detail::to_libint2_operator_params(oper_type, *unitcell_));

    p_screener =
        detail::make_screener(this->world(), screen_engine_pool, bases,
                              this->screen(), this->screen_threshold());
  }
}

template <typename Tile, typename Policy>
void PeriodicAOFactory<Tile, Policy>::parse_two_body_four_center_periodic(
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

  const auto &basis_registry = *this->basis_registry();

  auto bra_basis0 = detail::index_to_basis(basis_registry, bra_index0);
  auto bra_basis1 = detail::index_to_basis(basis_registry, bra_index1);
  auto ket_basis0 = detail::index_to_basis(basis_registry, ket_index0);
  auto ket_basis1 = detail::index_to_basis(basis_registry, ket_index1);

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

  p_screener = detail::make_screener(this->world(), engine_pool, bases,
                                     this->screen(), this->screen_threshold());
}

template <typename Tile, typename Policy>
template <typename E>
typename PeriodicAOFactory<Tile, Policy>::TArray
PeriodicAOFactory<Tile, Policy>::sparse_complex_integrals(
    madness::World &world, ShrPool<E> shr_pool, BasisVector const &bases,
    std::shared_ptr<Screener> screen, Op op) {
  auto time0 = mpqc::now(world, true);
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

  auto time_f0 = mpqc::now(world, true);
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

template <typename Tile, typename Policy>
template <typename E>
typename PeriodicAOFactory<Tile, Policy>::DirectTArray
PeriodicAOFactory<Tile, Policy>::direct_sparse_complex_integrals(
    madness::World &world, ShrPool<E> shr_pool, BasisVector const &bases,
    std::shared_ptr<Screener> screen, Op op) {
  // Build the Trange and Shape Tensor
  auto trange = detail::create_trange(bases);

  auto pmap = Policy::default_pmap(world, trange.tiles_range().volume());

  TA::TensorF tile_norms = screen->norm_estimate(world, bases, *pmap);

  TA::SparseShape<float> shape(world, tile_norms, trange);

  // Copy the Bases for the Integral Builder
  auto shr_bases = std::make_shared<BasisVector>(bases);

  // Make a pointer to an Integral builder.  Doing this because we want to use
  // it in Tasks.
  auto builder = make_direct_integral_builder(world, std::move(shr_pool),
                                              std::move(shr_bases),
                                              std::move(screen), std::move(op));

  auto direct_array = DirectArray<Tile, Policy, E>(std::move(builder));
  auto builder_ptr = direct_array.builder();
  using DirectTileType = DirectTile<Tile, E>;

  auto task_f = [=](int64_t ord, detail::IdxVec idx, TA::Range rng) {

    return DirectTileType(idx, std::move(rng), std::move(builder_ptr));

  };

  TA::DistArray<DirectTileType, TA::SparsePolicy> out(world, trange,
                                                      std::move(shape), pmap);

  for (auto const &ord : *pmap) {
    if (!out.is_zero(ord)) {
      detail::IdxVec idx = trange.tiles_range().idx(ord);
      auto range = trange.make_tile_range(ord);
      auto tile_fut = world.taskq.add(task_f, ord, idx, range);
      out.set(ord, tile_fut);
    }
  }
  world.gop.fence();

  direct_array.set_array(std::move(out));
  return direct_array;
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

#if TA_DEFAULT_POLICY == 1
extern template class PeriodicAOFactory<TA::TensorD, TA::SparsePolicy>;
#endif

}  // namespace gaussian
}  // namespace lcao
}  // namespace mpqc
#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_LCAO_FACTORY_PERIODIC_AO_FACTORY_H_
