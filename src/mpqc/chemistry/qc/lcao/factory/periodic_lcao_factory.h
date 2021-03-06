#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_LCAO_FACTORY_PERIODIC_LCAO_FACTORY_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_LCAO_FACTORY_PERIODIC_LCAO_FACTORY_H_

#include "mpqc/chemistry/qc/lcao/factory/periodic_ao_factory.h"

namespace mpqc {
namespace lcao {

template <typename Tile, typename Policy>
class PeriodicLCAOFactory;

template <typename Tile, typename Policy>
using PeriodicLCAOFactoryBase =
    Factory<TA::DistArray<Tile, Policy>, TA::DistArray<Tile, Policy>>;

/*!
 * \brief This constructs a PeriodicLCAOFactory object
 */
template <typename Tile, typename Policy>
std::shared_ptr<PeriodicLCAOFactory<Tile, Policy>>
construct_periodic_lcao_factory(const KeyVal &kv) {
  std::shared_ptr<PeriodicLCAOFactory<Tile, Policy>> periodic_lcao_factory;
  if (kv.exists_class("wfn_world:periodic_lcao_factory")) {
    periodic_lcao_factory = kv.class_ptr<PeriodicLCAOFactory<Tile, Policy>>(
        "wfn_world:periodic_lcao_factory");
  } else {
    periodic_lcao_factory =
        std::make_shared<PeriodicLCAOFactory<Tile, Policy>>(kv);
    std::shared_ptr<DescribedClass> ao_factory_base = periodic_lcao_factory;
    KeyVal &kv_nonconst = const_cast<KeyVal &>(kv);
    kv_nonconst.keyval("wfn_world")
        .assign("periodic_lcao_factory", ao_factory_base);
  }
  return periodic_lcao_factory;
}

template <typename Tile, typename Policy>
class PeriodicLCAOFactory : public PeriodicLCAOFactoryBase<Tile, Policy> {
 public:
  using TArray = TA::DistArray<Tile, Policy>;
  using AOFactoryType = gaussian::PeriodicAOFactory<Tile, Policy>;

  /*!
   * \brief KeyVal constructor for PeriodicLCAOFactory
   * \param kv the KeyVal object
   */
  PeriodicLCAOFactory(const KeyVal &kv)
      : PeriodicLCAOFactoryBase<Tile, Policy>(kv),
        pao_factory_(
            gaussian::construct_periodic_ao_factory<TA::TensorD, Policy>(kv)) {
    std::string prefix = "";
    if (kv.exists("wfn_world") || kv.exists_class("wfn_world"))
      prefix = "wfn_world:";

    // Molecule was already created at this path, bypass registry and construct
    // UnitCell
    unitcell_ = kv.class_ptr<UnitCell>(prefix + "atoms", true);

    dcell_ = unitcell_->dcell();

    R_max_ = pao_factory_->R_max();
    RJ_max_ = pao_factory_->RJ_max();
    RD_max_ = pao_factory_->RD_max();
    R_size_ = pao_factory_->R_size();
    RJ_size_ = pao_factory_->RJ_size();
    RD_size_ = pao_factory_->RD_size();

    // seperate k_points from zRHF::k_points
    using ::mpqc::detail::k_ord_idx;
    nk_ = decltype(nk_)(kv.value<std::vector<int>>("ref:k_points").data());
    k_size_ = 1 + k_ord_idx(nk_(0) - 1, nk_(1) - 1, nk_(2) - 1, nk_);

    auto orbital_space_registry =
        std::make_shared<OrbitalSpaceRegistry<TArray>>();

    this->set_orbital_registry(orbital_space_registry);
  }

  TArray compute_direct(const Formula &formula) override {
    throw ProgrammingError("Not Implemented!", __FILE__, __LINE__);
  }

  /*!
   * \brief This computes integrals by Formula
   * \param formula the desired Formula type
   * \return the TA::DistArray object
   */
  TArray compute(const Formula &formula) override;

  using PeriodicLCAOFactoryBase<Tile, Policy>::compute;
  using PeriodicLCAOFactoryBase<Tile, Policy>::compute_direct;

  /// return reference to PeriodicAOFactory object
  AOFactoryType &pao_factory() const { return *pao_factory_; }

 private:
  /// compute integrals that has two dimensions
  TArray compute2(const Formula &formula_string);

  /// compute fock integrals <p|F|q>
  TArray compute_fock_component(Formula &ao_formula,
                                std::shared_ptr<gaussian::Basis> bra_basis,
                                std::shared_ptr<gaussian::Basis> ket_basis);

  /// compute integrals that has four dimensions
  TArray compute4(const Formula &formula_string);

  /// find the correct range of summation for σ in <μ0 νR|ρR_j σ(R_j+R)>
  std::vector<int64_t> restricted_lattice_range(int64_t R, int64_t RJ);

  std::shared_ptr<AOFactoryType> pao_factory_;
  std::shared_ptr<UnitCell> unitcell_;

  /// unitcell and pbc information
  Vector3i R_max_;   ///> range of expansion of Bloch Gaussians in AO Gaussians
  Vector3i RJ_max_;  ///> range of Coulomb operation
  Vector3i RD_max_;  ///> range of density representation
  Vector3i nk_ = {1, 1, 1};  ///> # of k points in each direction
  Vector3d dcell_;           ///> direct unit cell params (in a.u.)
  int64_t
      R_size_;  ///> cardinal # of lattices included in Bloch Gaussian expansion
  int64_t RJ_size_;  ///> cardinal # of lattices included in Coulomb operation
  int64_t
      RD_size_;  ///> cardinal # of lattices included in density representation
  int64_t k_size_;  ///> cardinal # of k points
};

template <typename Tile, typename Policy>
typename PeriodicLCAOFactory<Tile, Policy>::TArray
PeriodicLCAOFactory<Tile, Policy>::compute(const Formula &formula) {
  auto iter = this->registry_.find(formula);

  auto &world = this->world();

  TArray result;

  if (iter != this->registry_.end()) {
    result = iter->second;
    utility::print_par(world, "Retrieved Periodic LCAO Integral: ",
                       utility::to_string(formula.string()));
    double size = mpqc::detail::array_size(result);
    utility::print_par(world, " Size: ", size, " GB\n");
    return result;
  }

  if (formula.rank() == 2) {
    result = compute2(formula);
    this->registry_.insert(formula, result);
  } else if (formula.rank() == 4) {
    result = compute4(formula);
    this->registry_.insert(formula, result);
  }

  madness::print_meminfo(
      world.rank(),
      "Periodic LCAOFactory: " + utility::to_string(formula.string()));
  return result;
}

template <typename Tile, typename Policy>
typename PeriodicLCAOFactory<Tile, Policy>::TArray
PeriodicLCAOFactory<Tile, Policy>::compute2(const Formula &formula) {
  auto &world = this->world();
  auto time0 = mpqc::now(world, this->accurate_time_);

  TArray result;

  // get AO formula
  auto ao_formula = detail::lcao_to_ao(formula, this->orbital_registry());

  // get AO bases
  auto bra_index = ao_formula.bra_indices()[0];
  auto ket_index = ao_formula.ket_indices()[0];

  const auto &basis_registry = *pao_factory_->basis_registry();

  auto bra_basis = gaussian::detail::index_to_basis(basis_registry, bra_index);
  auto ket_basis = gaussian::detail::index_to_basis(basis_registry, ket_index);

  TA_ASSERT(bra_basis != nullptr);
  TA_ASSERT(ket_basis != nullptr);

  // compute periodic AO integrals
  TArray pao_ints;
  if (ao_formula.oper().is_fock()) {
    auto v_formula = ao_formula;
    v_formula.set_operator_type(Operator::Type::Nuclear);
    auto t_formula = ao_formula;
    t_formula.set_operator_type(Operator::Type::Kinetic);
    auto j_formula = ao_formula;
    j_formula.set_operator_type(Operator::Type::J);
    auto k_formula = ao_formula;
    k_formula.set_operator_type(Operator::Type::K);

    auto v = compute_fock_component(v_formula, bra_basis, ket_basis);
    auto t = compute_fock_component(t_formula, bra_basis, ket_basis);
    auto j = compute_fock_component(j_formula, bra_basis, ket_basis);
    auto k = compute_fock_component(k_formula, bra_basis, ket_basis);

    pao_ints("p, q") = v("p, q") + t("p, q") + 2.0 * j("p, q") - k("p, q");
    // Symmetrize AO fock
    pao_ints("p, q") = 0.5 * (pao_ints("p, q") + pao_ints("q, p"));
  } else {
    using ::mpqc::detail::direct_vector;
    using ::mpqc::lcao::gaussian::detail::shift_basis_origin;
    using ::mpqc::lcao::gaussian::detail::to_libint2_operator;
    using ::mpqc::lcao::gaussian::detail::to_libint2_operator_params;
    using ::mpqc::lcao::gaussian::make_engine_pool;
    using ::mpqc::utility::make_array_of_refs;

    // shift AO bases and compute AO integrals
    for (auto R = 0; R < R_size_; ++R) {
      auto vec_R = direct_vector(R, R_max_, dcell_);

      auto bra = bra_basis;
      auto ket = shift_basis_origin(*ket_basis, vec_R);

      auto bases = mpqc::lcao::gaussian::BasisVector{{*bra, *ket}};
      auto engine_pool = mpqc::lcao::gaussian::make_engine_pool(
          to_libint2_operator(ao_formula.oper().type()),
          make_array_of_refs(bases[0], bases[1]), libint2::BraKet::x_x,
          to_libint2_operator_params(ao_formula.oper().type(), *unitcell_));

      auto ao_int = pao_factory_->compute_integrals(world, engine_pool, bases);

      if (R == 0)
        pao_ints("p, q") = ao_int("p, q");
      else
        pao_ints("p, q") += ao_int("p, q");
    }
  }
  auto aobuild_time1 = mpqc::now(world, this->accurate_time_);
  auto aobuild_duration = mpqc::duration_in_s(time0, aobuild_time1);

  // get MO coefficients
  auto left_index = formula.bra_indices()[0];
  if (left_index.is_lcao()) {
    auto &left = this->orbital_registry().retrieve(left_index);
    result("i, q") = pao_ints("p, q") * left("p, i");
  }

  auto right_index = formula.ket_indices()[0];
  if (right_index.is_lcao()) {
    auto &right = this->orbital_registry().retrieve(right_index);
    result("p, i") = result("p, q") * right("q, i");
  }

  auto time1 = mpqc::now(world, this->accurate_time_);
  auto trans_duration = mpqc::duration_in_s(aobuild_time1, time1);
  auto duration = mpqc::duration_in_s(time0, time1);

  result.truncate();

  ExEnv::out0() << " Transformed Gamma-Point Periodic LCAO Integral: "
                << utility::to_string(formula.string()) << std::endl;
  double size = mpqc::detail::array_size(result);
  ExEnv::out0() << " Size: " << size << " GB" << std::endl;
  ExEnv::out0() << " Time: " << duration << " s" << std::endl;
  ExEnv::out0() << "    PAO build time: " << aobuild_duration << " s"
                << std::endl;
  ExEnv::out0() << "    PAO->CO transform time: " << trans_duration << " s"
                << std::endl;

  return result;
}

template <typename Tile, typename Policy>
typename PeriodicLCAOFactory<Tile, Policy>::TArray
PeriodicLCAOFactory<Tile, Policy>::compute_fock_component(
    Formula &ao_formula, std::shared_ptr<gaussian::Basis> bra_basis,
    std::shared_ptr<gaussian::Basis> ket_basis) {
  ExEnv::out0() << indent << "Computing Periodic AO integral "
                << utility::to_string(ao_formula.string())
                << " for PeriodicLCAOFactory." << std::endl;

  auto &world = this->world();
  auto time0 = mpqc::now(world, this->accurate_time_);
  TArray result_ta;
  TArray ao_int;

  using ::mpqc::detail::direct_vector;
  using ::mpqc::detail::shift_mol_origin;
  using ::mpqc::lcao::gaussian::detail::shift_basis_origin;
  using ::mpqc::lcao::gaussian::detail::to_libint2_operator;
  using ::mpqc::lcao::gaussian::detail::to_libint2_operator_params;
  using ::mpqc::lcao::gaussian::make_engine_pool;
  using ::mpqc::utility::make_array_of_refs;

  for (auto R = 0; R < R_size_; ++R) {
    auto vec_R = direct_vector(R, R_max_, dcell_);
    auto bra = bra_basis;
    auto ket = shift_basis_origin(*ket_basis, vec_R);

    if (ao_formula.oper().type() == Operator::Type::Kinetic) {
      auto bases = mpqc::lcao::gaussian::BasisVector{{*bra, *ket}};

      auto engine_pool = mpqc::lcao::gaussian::make_engine_pool(
          to_libint2_operator(ao_formula.oper().type()),
          make_array_of_refs(bases[0], bases[1]), libint2::BraKet::x_x,
          to_libint2_operator_params(ao_formula.oper().type(), *unitcell_));

      ao_int = pao_factory_->compute_integrals(world, engine_pool, bases);

    } else if (ao_formula.oper().type() == Operator::Type::Nuclear) {
      auto bases = mpqc::lcao::gaussian::BasisVector{{*bra, *ket}};

      for (auto RJ = 0; RJ < RJ_size_; ++RJ) {
        auto shift_mol = direct_vector(RJ, RJ_max_, dcell_);
        auto shifted_mol = shift_mol_origin(*unitcell_, shift_mol);
        auto engine_pool = mpqc::lcao::gaussian::make_engine_pool(
            to_libint2_operator(ao_formula.oper().type()),
            make_array_of_refs(bases[0], bases[1]), libint2::BraKet::x_x,
            to_libint2_operator_params(ao_formula.oper().type(), *shifted_mol));
        if (RJ == 0)
          ao_int = pao_factory_->compute_integrals(world, engine_pool, bases);
        else
          ao_int("p, q") +=
              pao_factory_->compute_integrals(world, engine_pool, bases)("p, q");
      }
    } else if (ao_formula.oper().type() == Operator::Type::J) {
      // change operator type to Coulomb for computing integrals
      ao_formula.set_operator_type(Operator::Type::Coulomb);
      auto D = pao_factory_->get_density();

      for (auto RJ = 0; RJ < RJ_size_; ++RJ) {
        auto vec_RJ = direct_vector(RJ, RJ_max_, dcell_);

        auto bra0 = bra_basis;
        auto bra1 = shift_basis_origin(*ket_basis, vec_R);
        auto ket0 = shift_basis_origin(*bra_basis, vec_RJ);
        auto ket1 = shift_basis_origin(*ket_basis, vec_RJ, RD_max_, dcell_);

        // construct bases, engine_pool, and screener for AO integral computing
        auto bases =
            mpqc::lcao::gaussian::BasisVector{{*bra0, *bra1, *ket0, *ket1}};
        auto engine_pool = mpqc::lcao::gaussian::make_engine_pool(
            to_libint2_operator(ao_formula.oper().type()),
            make_array_of_refs(bases[0], bases[1], bases[2], bases[3]),
            libint2::BraKet::xx_xx,
            to_libint2_operator_params(ao_formula.oper().type(), *unitcell_));
        auto p_screener = gaussian::detail::make_screener(
            world, engine_pool, bases, pao_factory_->screen(),
            pao_factory_->screen_threshold());

        // compute AO-based integrals
        auto J = pao_factory_->compute_integrals(world, engine_pool, bases,
                                                p_screener);
        // sum over RJ
        if (RJ == 0)
          ao_int("p, q") = J("p, q, r, s") * D("r, s");
        else
          ao_int("p, q") += J("p, q, r, s") * D("r, s");
      }
      // change operator type back to J for future iterations
      ao_formula.set_operator_type(Operator::Type::J);

    } else if (ao_formula.oper().type() == Operator::Type::K) {
      // change operator type to Coulomb for computing integrals
      ao_formula.set_operator_type(Operator::Type::Coulomb);
      auto D = pao_factory_->get_density();

      for (auto RJ = 0; RJ < RJ_size_; ++RJ) {
        auto vec_RJ = direct_vector(RJ, RJ_max_, dcell_);

        auto bra0 = bra_basis;
        auto bra1 = shift_basis_origin(*ket_basis, vec_RJ);
        auto ket0 = shift_basis_origin(*ket_basis, vec_R);
        auto ket1 = shift_basis_origin(*ket_basis, vec_RJ, RD_max_, dcell_);

        // construct bases, engine_pool, and screener for AO integral computing
        auto bases =
            mpqc::lcao::gaussian::BasisVector{{*bra0, *bra1, *ket0, *ket1}};
        auto engine_pool = mpqc::lcao::gaussian::make_engine_pool(
            to_libint2_operator(ao_formula.oper().type()),
            make_array_of_refs(bases[0], bases[1], bases[2], bases[3]),
            libint2::BraKet::xx_xx,
            to_libint2_operator_params(ao_formula.oper().type(), *unitcell_));
        auto p_screener = gaussian::detail::make_screener(
            world, engine_pool, bases, pao_factory_->screen(),
            pao_factory_->screen_threshold());

        // compute AO-based integrals
        auto K = pao_factory_->compute_integrals(world, engine_pool, bases,
                                                p_screener);
        // sum over RJ
        if (RJ == 0)
          ao_int("p, q") = K("p, r, q, s") * D("r, s");
        else
          ao_int("p, q") += K("p, r, q, s") * D("r, s");
      }
      // change operator type back to K for future iterations
      ao_formula.set_operator_type(Operator::Type::K);
    }

    if (R == 0)
      result_ta("p, q") = ao_int("p, q");
    else
      result_ta("p, q") += ao_int("p, q");
  }

  auto time1 = mpqc::now(world, this->accurate_time_);
  auto duration = mpqc::duration_in_s(time0, time1);

  ExEnv::out0() << "    Time: " << duration << " s" << std::endl;

  return result_ta;
}

template <typename Tile, typename Policy>
typename PeriodicLCAOFactory<Tile, Policy>::TArray
PeriodicLCAOFactory<Tile, Policy>::compute4(const Formula &formula) {
  auto &world = this->world();
  auto time0 = mpqc::now(world, this->accurate_time_);

  TArray result;

  // get AO formula
  auto ao_formula = detail::lcao_to_ao(formula, this->orbital_registry());

  // get AO bases
  auto bra_index0 = ao_formula.bra_indices()[0];
  auto bra_index1 = ao_formula.bra_indices()[1];
  auto ket_index0 = ao_formula.ket_indices()[0];
  auto ket_index1 = ao_formula.ket_indices()[1];

  const auto &basis_registry = *pao_factory_->basis_registry();

  auto bra_basis0 =
      gaussian::detail::index_to_basis(basis_registry, bra_index0);
  auto bra_basis1 =
      gaussian::detail::index_to_basis(basis_registry, bra_index1);
  auto ket_basis0 =
      gaussian::detail::index_to_basis(basis_registry, ket_index0);
  auto ket_basis1 =
      gaussian::detail::index_to_basis(basis_registry, ket_index1);

  TA_ASSERT(bra_basis0 != nullptr);
  TA_ASSERT(bra_basis1 != nullptr);
  TA_ASSERT(ket_basis0 != nullptr);
  TA_ASSERT(ket_basis1 != nullptr);

  using ::mpqc::detail::direct_vector;
  using ::mpqc::lcao::gaussian::detail::shift_basis_origin;
  using ::mpqc::lcao::gaussian::detail::to_libint2_operator;
  using ::mpqc::lcao::gaussian::detail::to_libint2_operator_params;
  using ::mpqc::lcao::gaussian::make_engine_pool;
  using ::mpqc::utility::make_array_of_refs;

  double pao_build_time = 0.0;
  double sum_r_time = 0.0;
  // shift AO bases and compute AO integrals
  TArray pao_ints;
  auto sum_count = 0;
  for (auto R = 0; R < R_size_; ++R) {
    auto vec_R = direct_vector(R, R_max_, dcell_);
    for (auto RJ = 0; RJ < RJ_size_; ++RJ) {
      auto vec_RJ = direct_vector(RJ, RJ_max_, dcell_);

      auto selected_RD_list = restricted_lattice_range(R, RJ);
      for (auto RD_idx = 0ul; RD_idx < selected_RD_list.size(); ++RD_idx) {
        auto RD = selected_RD_list[RD_idx];
        auto vec_RD = direct_vector(RD, RD_max_, dcell_);

        auto t_pao0 = mpqc::now(world, this->accurate_time_);
        auto bra0 = bra_basis0;
        auto bra1 = shift_basis_origin(*bra_basis1, vec_RJ);
        auto ket0 = shift_basis_origin(*ket_basis0, vec_R);
        auto ket1 = shift_basis_origin(*ket_basis1, vec_RJ + vec_RD);

        // compute as in chemists' notation
        // construct bases, engine_pool, and screener for AO integral computing
        auto bases =
            mpqc::lcao::gaussian::BasisVector{{*bra0, *ket0, *bra1, *ket1}};
        auto engine_pool = mpqc::lcao::gaussian::make_engine_pool(
            to_libint2_operator(ao_formula.oper().type()),
            make_array_of_refs(bases[0], bases[1], bases[2], bases[3]),
            libint2::BraKet::xx_xx,
            to_libint2_operator_params(ao_formula.oper().type(), *unitcell_));
        auto p_screener = gaussian::detail::make_screener(
            world, engine_pool, bases, pao_factory_->screen(),
            pao_factory_->screen_threshold());

        // compute AO-based integrals
        auto ao_int = pao_factory_->compute_integrals(world, engine_pool, bases,
                                                     p_screener);
        auto t_pao1 = mpqc::now(world, this->accurate_time_);
        pao_build_time += mpqc::duration_in_s(t_pao0, t_pao1);

        // return as in physicists' notation
        // sum over R, RJ, and RD
        if (sum_count == 0)
          pao_ints("p, q, r, s") = ao_int("p, r, q, s");
        else
          pao_ints("p, q, r, s") += ao_int("p, r, q, s");
        auto t_sum_r1 = mpqc::now(world, this->accurate_time_);
        sum_r_time += mpqc::duration_in_s(t_pao1, t_sum_r1);

        sum_count++;
      }
    }
  }

  auto aobuild_time1 = mpqc::now(world, this->accurate_time_);
  auto aobuild_duration = mpqc::duration_in_s(time0, aobuild_time1);

  // get MO coefficients
  auto left_index1 = formula.bra_indices()[0];
  if (left_index1.is_lcao()) {
    auto &left1 = this->orbital_registry().retrieve(left_index1);
    result("i, q, r, s") = pao_ints("p, q, r, s") * left1("p, i");
  }

  auto left_index2 = formula.bra_indices()[1];
  if (left_index2.is_lcao()) {
    auto &left2 = this->orbital_registry().retrieve(left_index2);
    result("p, i, r, s") = result("p, q, r, s") * left2("q, i");
  }

  auto right_index1 = formula.ket_indices()[0];
  if (right_index1.is_lcao()) {
    auto &right1 = this->orbital_registry().retrieve(right_index1);
    result("p, q, i, s") = result("p, q, r, s") * right1("r, i");
  }

  auto right_index2 = formula.ket_indices()[1];
  if (right_index2.is_lcao()) {
    auto &right2 = this->orbital_registry().retrieve(right_index2);
    result("p, q, r, i") = result("p, q, r, s") * right2("s, i");
  }

  auto time1 = mpqc::now(world, this->accurate_time_);
  auto trans_duration = mpqc::duration_in_s(aobuild_time1, time1);
  auto duration = mpqc::duration_in_s(time0, time1);

  result.truncate();

  ExEnv::out0() << " Transformed Gamma-Point Periodic LCAO Integral: "
                << utility::to_string(formula.string()) << std::endl;
  double size = mpqc::detail::array_size(result);
  ExEnv::out0() << " Size: " << size << " GB" << std::endl;
  ExEnv::out0() << " Time: " << duration << " s" << std::endl;
  ExEnv::out0() << "    PAO build total time:  " << aobuild_duration << " s"
                << std::endl;
  ExEnv::out0() << "        part 1: " << pao_build_time << " s" << std::endl;
  ExEnv::out0() << "        part 2: " << sum_r_time << " s" << std::endl;
  ExEnv::out0() << "    PAO->CO transform time: " << trans_duration << " s"
                << std::endl;

  return result;
}

template <typename Tile, typename Policy>
std::vector<int64_t>
PeriodicLCAOFactory<Tile, Policy>::restricted_lattice_range(int64_t R,
                                                            int64_t RJ) {
  std::vector<int64_t> result;
  using ::mpqc::detail::direct_vector;

  auto vec_max = direct_vector(RJ_size_ - 1, RJ_max_, dcell_);
  auto distance_max = vec_max.norm();
  auto vec_R = direct_vector(R, R_max_, dcell_);
  auto vec_RJ = direct_vector(RJ, RJ_max_, dcell_);

  for (int64_t R_select = 0; R_select < RD_size_; ++R_select) {
    auto vec_R_select = direct_vector(R_select, RD_max_, dcell_);
    auto distance = (vec_RJ + vec_R_select - vec_R).norm();
    if (distance <= distance_max) result.push_back(R_select);
  }

  return result;
}

#if TA_DEFAULT_POLICY == 1
extern template class PeriodicLCAOFactory<TA::TensorZ, TA::SparsePolicy>;
#endif

}  // namespace lcao
}  // namespace mpqc
#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_LCAO_FACTORY_PERIODIC_LCAO_FACTORY_H_
