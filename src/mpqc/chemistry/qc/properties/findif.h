#ifndef SRC_MPQC_CHEMISTRY_QC_PROPERTIES_FINDIF_H_
#define SRC_MPQC_CHEMISTRY_QC_PROPERTIES_FINDIF_H_

#include "mpqc/chemistry/qc/properties/property.h"
#include "mpqc/math/function/findif.h"

namespace mpqc {

/// Evaluates FiniteDifferenceDerivative with respect to molecular coordinates
template <size_t Order, typename Value>
class MolecularFiniteDifferenceDerivative
    : public math::FiniteDifferenceDerivative<Order, Value,
                                              MolecularCoordinates>,
      public Property {
 public:
  // clang-format off
  /**
   * @brief The KeyVal constructor
   * @param kv the KeyVal object to be queried
   *
   * The KeyVal object will be queried for all keywords of the KeyVal ctor of the math::FiniteDifferenceDerivative class,
   * as well as the following keywords:
   * | KeyWord | Type | Default| Description |
   * |---------|------|--------|-------------|
   * | coords | MolecularCoordinates | CartMolecularCoordinates (must specify the \c atoms keyword) | specifies which molecular coordinates will be used to compute the derivative. The only valid choice is CartMolecularCoordinates |
   * | atoms | Molecule | none | if \c coords is not specified, must specify this keyword to be able to construct CartMolecularCoordinates |
   */
  // clang-format on
  MolecularFiniteDifferenceDerivative(const KeyVal& kv)
      : math::FiniteDifferenceDerivative<Order, Value, MolecularCoordinates>(
            kv, make_coords(kv), default_target_precision_) {}

  constexpr static double default_target_precision_ = 1e-6;

  void write(KeyVal& kv) const override {
    auto kv_val = kv.keyval("value");
    this->get_value().value()->write(kv_val);
  }

 private:
  /// overrides Property::evaluate()
  void evaluate() override { this->value(); }

  /// constructs molecular coordinates
  static std::shared_ptr<MolecularCoordinates> make_coords(const KeyVal& kv) {
    // looks for
    // 1. coords keyword
    // 2. make CartesianMolecularCoordinates if atoms or molecule keyword given
    auto coords = kv.class_ptr<MolecularCoordinates>("coords");
    if (!coords) {
      std::shared_ptr<Molecule> atoms = kv.class_ptr<Molecule>("atoms");
      if (!atoms) atoms = kv.class_ptr<Molecule>("molecule");
      if (atoms)
        coords = std::make_shared<CartMolecularCoordinates>(atoms);
      else {
        throw InputError(
            "MolecularFiniteDifferenceDerivative expects keywords \"coords\" "
            "or \"atoms\", and neither was given",
            __FILE__, __LINE__);
      }
    }
    return coords;
  }
};

}  // namespace mpqc

#endif /* SRC_MPQC_CHEMISTRY_QC_PROPERTIES_FINDIF_H_ */
