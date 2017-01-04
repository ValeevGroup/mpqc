#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_PERIODIC_LCAO_FACTORY_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_PERIODIC_LCAO_FACTORY_H_

#include "mpqc/chemistry/qc/integrals/lcao_factory.h"
#include "mpqc/chemistry/qc/integrals/periodic_ao_factory.h"

namespace mpqc {
namespace lcao {

template <typename Tile, typename Policy>
class PeriodicLCAOFactory;

namespace detail {
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

}  // namespace detail

template <typename Tile, typename Policy>
class PeriodicLCAOFactory : public LCAOFactory<TA::TensorD, Policy> {
 public:
  using TArray = TA::DistArray<Tile, Policy>;
  using AOFactoryType = gaussian::PeriodicAOFactory<Tile, Policy>;

  /*!
   * \brief KeyVal constructor for PeriodicLCAOFactory
   * \param kv the KeyVal object
   */
  PeriodicLCAOFactory(const KeyVal &kv)
      : LCAOFactory<TA::TensorD, Policy>(kv),
        pao_factory_(
            *gaussian::construct_periodic_ao_factory<Tile, Policy>(kv)),
        unitcell_(
            *kv.keyval("wfn_world:molecule").template class_ptr<UnitCell>()),
        orbital_space_registry_(
            std::make_shared<OrbitalSpaceRegistry<TArray>>()) {
    dcell_ = unitcell_.dcell();

    R_max_ = pao_factory_.R_max();
    RJ_max_ = pao_factory_.RJ_max();
    RD_max_ = pao_factory_.RD_max();
    R_size_ = pao_factory_.R_size();
    RJ_size_ = pao_factory_.RJ_size();
    RD_size_ = pao_factory_.RD_size();

    // seperate k_points from zRHF::k_points
    nk_ = decltype(nk_)(kv.value<std::vector<int>>("k_points").data());
    k_size_ = 1 + detail::k_ord_idx(nk_(0) - 1, nk_(1) - 1, nk_(2) - 1, nk_);
  }

  /// wrapper to compute function
  TArray compute(const std::wstring &formula_string);

  /*!
   * \brief This computes integrals by Formula
   * \param formula the desired Formula type
   * \return the TA::DistArray object
   */
  TArray compute(const Formula &formula);

  /// return reference to PeriodicAOFactory object
  AOFactoryType &pao_factory() const { return pao_factory_; }

 private:
  /// compute integrals that has two dimensions
  TArray compute2(const Formula &formula_string);

  /// compute integrals that has four dimensions
  TArray compute4(const Formula &formula_string);

  AOFactoryType &pao_factory_;
  UnitCell &unitcell_;
  std::shared_ptr<OrbitalSpaceRegistry<TArray>> orbital_space_registry_;

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
PeriodicLCAOFactory<Tile, Policy>::compute(const std::wstring &formula_string) {
  Formula formula(formula_string);
  return compute(formula);
}

template <typename Tile, typename Policy>
typename PeriodicLCAOFactory<Tile, Policy>::TArray
PeriodicLCAOFactory<Tile, Policy>::compute(const Formula &formula) {
  TArray result;

  ExEnv::out0() << "Computing Periodic LCAO integrals " << utility::to_string(formula.string()) << " ..." << std::endl;

  if (formula.rank() == 2) {
    result = compute2(formula);
  } else if (formula.rank() == 4) {
    result = compute4(formula);
  }

  return result;
}

template <typename Tile, typename Policy>
typename PeriodicLCAOFactory<Tile, Policy>::TArray
PeriodicLCAOFactory<Tile, Policy>::compute2(const Formula &formula) {
  TArray result;

  return result;
}

template <typename Tile, typename Policy>
typename PeriodicLCAOFactory<Tile, Policy>::TArray
PeriodicLCAOFactory<Tile, Policy>::compute4(const Formula &formula) {
  mpqc::time_point time0 = mpqc::now(this->world_, this->accurate_time_);

  TArray result;

  // get AO formula
  auto ao_formula = this->mo_to_ao(formula);

  // get AO bases
  auto bra_index0 = ao_formula.bra_indices()[0];
  auto bra_index1 = ao_formula.bra_indices()[1];
  auto ket_index0 = ao_formula.ket_indices()[0];
  auto ket_index1 = ao_formula.ket_indices()[1];

  auto bra_basis0 = pao_factory_.index_to_basis(bra_index0);
  auto bra_basis1 = pao_factory_.index_to_basis(bra_index1);
  auto ket_basis0 = pao_factory_.index_to_basis(ket_index0);
  auto ket_basis1 = pao_factory_.index_to_basis(ket_index1);

  TA_ASSERT(bra_basis0 != nullptr);
  TA_ASSERT(bra_basis1 != nullptr);
  TA_ASSERT(ket_basis0 != nullptr);
  TA_ASSERT(ket_basis1 != nullptr);

  using ::mpqc::lcao::detail::direct_vector;
  using ::mpqc::lcao::gaussian::detail::shift_basis_origin;
  using ::mpqc::lcao::gaussian::detail::to_libint2_operator;
  using ::mpqc::lcao::gaussian::detail::to_libint2_operator_params;
  using ::mpqc::lcao::gaussian::make_engine_pool;
  using ::mpqc::utility::make_array_of_refs;

  // shift AO bases and compute AO integrals
  TArray pao_ints;
  auto sum_count = 0;
  for (auto R = 0; R < R_size_; ++R) {
    auto vec_R = direct_vector(R, R_max_, dcell_);
    for (auto RJ = 0; RJ < RJ_size_; ++RJ) {
      auto vec_RJ = direct_vector(RJ, RJ_max_, dcell_);
      for (auto RD = 0; RD < RD_size_; ++RD) {
        auto vec_RD = direct_vector(RD, RD_max_, dcell_);

        auto bra0 = bra_basis0;
        auto bra1 = shift_basis_origin(*bra_basis1, vec_R);
        auto ket0 = shift_basis_origin(*ket_basis0, vec_RJ);
        auto ket1 = shift_basis_origin(*ket_basis1, vec_RJ + vec_RD);

        auto bases =
            mpqc::lcao::gaussian::BasisVector{{*bra0, *bra1, *ket0, *ket1}};
        auto engine_pool = mpqc::lcao::gaussian::make_engine_pool(
            to_libint2_operator(ao_formula.oper().type()),
            make_array_of_refs(bases[0], bases[1], bases[2], bases[3]),
            libint2::BraKet::xx_xx,
            to_libint2_operator_params(ao_formula.oper().type(), pao_factory_,
                                       unitcell_));
        auto ao_int =
            pao_factory_.compute_integrals(this->world_, engine_pool, bases);

        if (sum_count == 0)
          pao_ints("p, q, r, s") = ao_int("p, q, r, s");
        else
          pao_ints("p, q, r, s") += ao_int("p, q, r, s");

        sum_count++;
      }
    }
  }

  // get MO coefficients
  auto left_index1 = ao_formula.bra_indices()[0];
  if (left_index1.is_mo()) {
    auto &left1 = orbital_space_registry_->retrieve(left_index1);
    result("i, q, r, s") = pao_ints("p, q, r, s") * left1("p, i");
  }

  auto left_index2 = ao_formula.bra_indices()[1];
  if (left_index2.is_mo()) {
    auto &left2 = orbital_space_registry_->retrieve(left_index2);
    result("p, i, r, s") = result("p, q, r, s") * left2("q, i");
  }

  auto right_index1 = ao_formula.ket_indices()[0];
  if (right_index1.is_mo()) {
    auto &right1 = orbital_space_registry_->retrieve(right_index1);
    result("p, q, i, s") = result("p, q, r, s") * right1("r, i");
  }

  auto right_index2 = ao_formula.ket_indices()[1];
  if (right_index2.is_mo()) {
    auto &right2 = orbital_space_registry_->retrieve(right_index2);
    result("p, q, r, i") = result("p, q, r, s") * right2("s, i");
  }

  result.truncate();

  mpqc::time_point time1 = mpqc::now(this->world_, this->accurate_time_);
  auto duration = mpqc::duration_in_s(time0, time1);
  double size = mpqc::detail::array_size(result);

  ExEnv::out0() << " Transformed Gamma-Point Periodic LCAO Integral: "
                << utility::to_string(formula.string()) << std::endl;
  ExEnv::out0() << " Size: " << size << " GB" << std::endl;
  ExEnv::out0() << " Time: " << duration << " s" << std::endl;

  return result;
}

}  // namespace lcao
}  // namespace mpqc
#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_PERIODIC_LCAO_FACTORY_H_
