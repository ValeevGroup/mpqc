#ifndef SRC_MPQC_CHEMISTRY_QC_PROPERTIES_FINDIF_H_
#define SRC_MPQC_CHEMISTRY_QC_PROPERTIES_FINDIF_H_

#include "mpqc/chemistry/qc/properties/property.h"
#include "mpqc/math/function/findif.h"

namespace mpqc {

/// Evaluates FiniteDifferenceDerivative with respect to molecular coordinates
template <size_t Order, typename Value>
class MolecularFiniteDifferenceDerivative
    : public math::FiniteDifferenceDerivative<Order, Value, MolecularCoordinates>,
      public Property {
 public:
  MolecularFiniteDifferenceDerivative(const KeyVal& kv)
      : math::FiniteDifferenceDerivative<Order, Value, MolecularCoordinates>(
            kv, (kv.class_ptr<MolecularCoordinates>("coords")
                     ? kv.class_ptr<MolecularCoordinates>("coords")
                     : std::make_shared<CartMolecularCoordinates>(
                           kv.class_ptr<Molecule>("molecule"))),
                           default_target_precision_) {}

  constexpr static double default_target_precision_ = 1e-6;

  void write(KeyVal& kv) const override {
    auto kv_val = kv.keyval("value");
    this->get_value().value()->write(kv_val);
  }

 private:
  /// overrides Property::evaluate()
  void evaluate() override { this->value(); }
};

}

#endif /* SRC_MPQC_CHEMISTRY_QC_PROPERTIES_FINDIF_H_ */
