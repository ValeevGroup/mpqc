#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_PERIODIC_LCAO_FACTORY_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_PERIODIC_LCAO_FACTORY_H_

#include "mpqc/chemistry/qc/integrals/lcao_factory.h"
#include "mpqc/chemistry/qc/integrals/periodic_ao_factory.h"

namespace mpqc {
namespace lcao {

namespace detail {}  // namespace detail

template <typename Tile, typename Policy>
class PeriodicLCAOFactory;

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
class PeriodicLCAOFactory : public LCAOFactory<TA::TensorD, Policy> {
 public:
  using TArray = TA::DistArray<Tile, Policy>;
  using AOFactoryType = gaussian::PeriodicAOFactory<Tile, Policy>;

  /*!
   * \brief KeyVal constructor for PeriodicLCAOFactory
   * \param kv the KeyVal object
   */
  PeriodicLCAOFactory(const KeyVal &kv) : LCAOFactory<TA::TensorD, Policy>(kv) {
    std::string prefix = "";
    if (kv.exists("wfn_world") || kv.exists_class("wfn_world")) {
      prefix = "wfn_world:";
    }
  }

  /// wrapper to compute function
  TArray compute(const std::wstring &formula_string);

  /*!
   * \brief This computes integrals by Formula
   * \param formula the desired Formula type
   * \return the TA::DistArray object
   */
  TArray compute(const Formula &formula);

 private:
  /// compute integrals that has two dimensions
  TArray compute2(const Formula &formula_string);

  /// compute integrals that has four dimensions
  TArray compute4(const Formula &formula_string);
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

    ExEnv::out0() << "Computing Periodic LCAO integrals ..." << std::endl;

    if (formula.rank() == 2) {
        result = compute2(formula);
        mo_formula_registry_.insert(formula, result);
    } else if (formula.rank() == 4) {
        result = compute4(formula);
        mo_formula_registry_.insert(formula, result);
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
    TArray result;

    return result;
}

}  // namespace lcao
}  // namespace mpqc
#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_INTEGRALS_PERIODIC_LCAO_FACTORY_H_
