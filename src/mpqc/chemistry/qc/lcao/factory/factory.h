//
// Created by Chong Peng on 2/13/17.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_LCAO_FACTORY_FACTORY_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_LCAO_FACTORY_FACTORY_H_

#include "mpqc/chemistry/qc/lcao/wfn/wfn_world.h"
#include "mpqc/util/keyval/keyval.h"
#include "mpqc/chemistry/qc/lcao/expression/formula_registry.h"

namespace mpqc {
namespace lcao {

template <typename Array, typename DirectArray>
class Factory : virtual public DescribedClass {
 public:

  Factory() = default;

  Factory(std::shared_ptr<WavefunctionWorld> wfn_world)
      : wfn_world_(wfn_world),
        registry_(),
        direct_registry_() {}

  /// @return MADNESS world
  madness::World& world() { return wfn_world_->world(); }

  /// @return a pointer to the molecule in the world
  const std::shared_ptr<Molecule>& atoms() const { return wfn_world_->atoms(); }

  /// @return a pointer to Basis Registry
  std::shared_ptr<gaussian::OrbitalBasisRegistry> const basis_registry() {
    return wfn_world_->basis_registry();
  }

  /// return const registry
  const FormulaRegistry<Array>& registry() const {
    return registry_;
  }

  /// return registry
  FormulaRegistry<Array>& registry() {
    return registry_;
  }

  /// return const registry
  const FormulaRegistry<DirectArray>& direct_registry() const {
    return direct_registry_;
  }

  /// return registry
  FormulaRegistry<DirectArray>& direct_registry() {
    return direct_registry_;
  }

  /// wrapper to compute function
  Array compute(const std::wstring& string){
    auto formula = Formula(string);
    return compute(formula);
  }

  /// wrapper to compute direct function
  DirectArray compute_direct(const std::wstring& str) {
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

  /// abstract functions
  /// compute array
  virtual Array compute(const Formula& formula) = 0;

  /// compute direct array
  virtual DirectArray compute_direct(const Formula& formula) = 0;


 protected:
  FormulaRegistry<Array> registry_;
  FormulaRegistry<DirectArray> direct_registry_;

 private:
  std::shared_ptr<WavefunctionWorld> wfn_world_;
};

}  // namespace lcao
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_LCAO_FACTORY_FACTORY_H_
