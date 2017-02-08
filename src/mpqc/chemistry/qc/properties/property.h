#ifndef SRC_MPQC_CHEMISTRY_QC_PROPERTIES_PROPERTY_H_
#define SRC_MPQC_CHEMISTRY_QC_PROPERTIES_PROPERTY_H_

#include <atomic>
#include <cstdint>
#include <memory>
#include <vector>

#include "mpqc/chemistry/molecule/coords.h"
#include "mpqc/chemistry/molecule/molecule.h"
#include "mpqc/chemistry/qc/wfn/wfn.h"
#include "mpqc/math/function/function.h"
#include "mpqc/math/function/taylor.h"
#include "mpqc/util/misc/task.h"
#include "mpqc/util/misc/provider.h"

/// top-level MPQC namespace
namespace mpqc {

// using TiledArray::detail::scalar_type;

/// this is the base for all properties that MPQC can compute via the input.
/// MPQC main will read KeyVal and search for a Property object, compute it
/// using the given wave function
class Property : public Task {
 public:
  /// evaluates this object
  virtual void evaluate() = 0;

  /// prints this object to \c kv
  /// @param kv the output KeyVal object
  virtual void write(KeyVal& kv) const = 0;

 private:
  /// implements Task::execute()
  void execute() override final { evaluate(); }
};

/**
 * \brief WavefunctionProperty computes a Taylor expansion of a molecular
 * property using a visiting Wavefunction .
 *
 * This is to be used as a base class for ALL properties of Wavefunction
 * classes.
 */

template <typename Value>
class WavefunctionProperty
    : public math::TaylorExpansionFunction<Value, MolecularCoordinates>,
      public Property {
 public:
  using base_type = math::TaylorExpansionFunction<Value, MolecularCoordinates>;
  using typename base_type::value_type;  // i.e.
                                         // TaylorExpansionCoefficients<Value>
  using typename base_type::function_base_type;

  // clang-format off
  /**
   * @brief The KeyVal constructor
   * @param kv the KeyVal object to be queried
   *
   * The KeyVal object will be queried for all keywords of the auxiliary KeyVal ctor of the MolecularTaylorExpansion class,
   * as well as the following keywords:
   * | KeyWord | Type | Default| Description |
   * |---------|------|--------|-------------|
   * | wfn | Wavefunction | none | the Wavefunction that will compute this  |
   * | coords | MolecularCoordinates | CartMolecularCoordinates from the Molecule object of the Wavefunction object | the molecular coordinates  |
   */
  // clang-format on
  WavefunctionProperty(const KeyVal& kv)
      : WavefunctionProperty(kv, base_type::default_precision_) {}

 protected:
  const std::shared_ptr<Wavefunction>& wfn() const { return wfn_; }

  virtual void do_evaluate() = 0;

  // clang-format off
  /**
   * @brief The auxiliary KeyVal constructor
   *
   * Like standard KeyVal ctor, except accepts extra defaults from derived classes. Queries same keywords.
   *
   * @param kv the KeyVal object to be queried
   * @param default_precision the default precision to be used by MolecularTaylorExpansion
   */
  // clang-format on
  WavefunctionProperty(const KeyVal& kv, double default_precision)
      : base_type(kv, (kv.class_ptr<MolecularCoordinates>("coord")
                           ? kv.class_ptr<MolecularCoordinates>("coord")
                           : std::make_shared<CartMolecularCoordinates>(
                                 kv.class_ptr<Wavefunction>("wfn")->atoms())),
                  default_precision) {
    wfn_ = kv.class_ptr<Wavefunction>("wfn");
    if (wfn_ == nullptr)
      throw InputError(
          "WavefunctionProperty did not receive a Wavefunction object",
          __FILE__, __LINE__, "wfn");
  }

  void write(KeyVal& kv) const override {
    auto kv_val = kv.keyval("value");
    this->get_value().value()->write(kv_val);
  }

 private:
  std::shared_ptr<Wavefunction> wfn_;

  void evaluate() override { this->compute(); }

  void compute() override {
    auto original_value_timestamp = this->get_value().timestamp();
    this->do_evaluate();
    auto updated_value_timestamp = this->get_value().timestamp();
    if (original_value_timestamp == updated_value_timestamp)
      throw ProgrammingError(
          "WavefunctionProperty::compute(): Wavefunction forgot to call "
          "set_value?",
          __FILE__, __LINE__);
  }
};

}  // namespace mpqc

#endif  // SRC_MPQC_CHEMISTRY_QC_PROPERTIES_PROPERTY_H_
