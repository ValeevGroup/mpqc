//
// Created by Chong Peng on 2/13/17.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_LCAO_FACTORY_FACTORY_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_LCAO_FACTORY_FACTORY_H_

#include "mpqc/chemistry/qc/lcao/expression/formula_registry.h"
#include "mpqc/chemistry/qc/lcao/wfn/wfn_world.h"
#include "mpqc/util/keyval/keyval.h"

namespace mpqc {
namespace lcao {

template <typename Array, typename DirectArray=Array>
class Factory : virtual public DescribedClass {
 public:
  Factory() = default;

  Factory(std::shared_ptr<WavefunctionWorld> wfn_world)
      : registry_(),
        direct_registry_(),
        orbital_space_registry_(nullptr),
        accurate_time_(false),
        wfn_world_(wfn_world) {}

  // clang-format off
  /**
   * @brief The KeyVal constructor.
   *
   * @param kv the KeyVal object, it will be queried for the following keywords:
   *
   *  | Keyword | Type | Default| Description |
   *  |---------|------|--------|-------------|
   *  |wfn_world| WavefunctionWorld | none | WavefunctionWorld object |
   *  |accurate_time|bool|false|if true, do fence at timing|
   */
  // clang-format on

  Factory(const KeyVal& kv)
      : Factory(kv.class_ptr<WavefunctionWorld>("wfn_world")) {
    accurate_time_ = kv.value<bool>("accurate_time", false);
  }

  /// @return MADNESS world
  madness::World& world() { return wfn_world_->world(); }

  /// @return a pointer to the molecule in the world
  const std::shared_ptr<Molecule>& atoms() const { return wfn_world_->atoms(); }

  /// @return a pointer to Basis Registry
  std::shared_ptr<gaussian::OrbitalBasisRegistry> const basis_registry() {
    return wfn_world_->basis_registry();
  }

  /// return const orbital registry
  const OrbitalSpaceRegistry<Array>& orbital_registry() const {
    TA_USER_ASSERT(orbital_space_registry_ != nullptr,
                   "OrbitalSpaceRegistry not initialized! \n");
    return *orbital_space_registry_;
  }

  /// return orbital registry
  OrbitalSpaceRegistry<Array>& orbital_registry() {
    TA_USER_ASSERT(orbital_space_registry_ != nullptr,
                   "OrbitalSpaceRegistry not initialized! \n");
    return *orbital_space_registry_;
  }

  void set_orbital_registry(
      const std::shared_ptr<OrbitalSpaceRegistry<Array>>& obs_registry) {
    orbital_space_registry_ = obs_registry;
  }

  /// return const registry
  const FormulaRegistry<Array>& registry() const { return registry_; }

  /// return registry
  FormulaRegistry<Array>& registry() { return registry_; }

  /// return const direct registry
  template <typename T = Array, typename U = DirectArray>
  const typename std::enable_if<!std::is_same<T, U>::value,
                                FormulaRegistry<DirectArray>>::type&
  direct_registry() const {
    return direct_registry_;
  }

  /// return direct_registry
  template <typename T = Array, typename U = DirectArray>
  const typename std::enable_if<!std::is_same<T, U>::value,
                                FormulaRegistry<DirectArray>>::type&
  direct_registry() {
    return direct_registry_;
  }

  /// wrapper to compute function
  Array compute(const std::wstring& string) {
    auto formula = Formula(string);
    return compute(formula);
  }

  /// wrapper to compute direct function
  template <typename T = Array, typename U = DirectArray>
  typename std::enable_if<!std::is_same<T, U>::value, DirectArray>::type
  compute_direct(const std::wstring& str) {
    auto formula = Formula(str);
    return compute_direct(formula);
  }

  /// compute with str and return expression
  TA::expressions::TsrExpr<Array, true> operator()(const std::wstring& str) {
    auto formula = Formula(str);
    auto array = compute(formula);
    auto& result = registry_.retrieve(formula);
    return result(formula.to_ta_expression());
  };

  /// obsolete Factory
  virtual void obsolete() {
    registry_.purge();
    direct_registry_.purge();
    if (orbital_space_registry_ != nullptr) {
      orbital_space_registry_->clear();
    }
  }

  bool accurate_time() const { return accurate_time_; }

  /// abstract functions
  /// compute array
  virtual Array compute(const Formula& formula) = 0;

  /// compute direct array
  //  virtual typename std::enable_if<!std::is_same<Array, DirectArray>::value,
  //                                  DirectArray>::type
  virtual DirectArray compute_direct(const Formula& formula) = 0;

  /// purge formula that contain Operator described by string \c str
  /// from registry and direct_registry
  virtual void purge_operator(const std::wstring& str) {
    Operator oper(str);
    Operator::Type oper_type = oper.type();

    registry_.purge_operator(oper_type);
    direct_registry_.purge_operator(oper_type);
  }

  /// purge formula that contain index described by string \c idx_str
  /// from registry and direct_registry
  virtual void purge_index(const std::wstring& idx_str) {
    OrbitalIndex index(idx_str);
    registry_.purge_index(index);
    direct_registry_.purge_index(index);
  }

  /// purge formula described by string \c str
  /// from mo_registry
  virtual void purge_formula(const std::wstring& str) {
    registry_.purge_formula(str);
    direct_registry_.purge_formula(str);
  }

 protected:
  /// registry for Array
  FormulaRegistry<Array> registry_;
  /// registry for DirectArray
  FormulaRegistry<DirectArray> direct_registry_;
  /// registry for Orbital Space
  std::shared_ptr<OrbitalSpaceRegistry<Array>> orbital_space_registry_;
  /// if do fence when time
  bool accurate_time_;

 private:
  std::shared_ptr<WavefunctionWorld> wfn_world_;
};

}  // namespace lcao
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_LCAO_FACTORY_FACTORY_H_
