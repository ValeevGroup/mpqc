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

/// top-level MPQC namespace
namespace mpqc {

// using TiledArray::detail::scalar_type;

/// computing a runtime-typed property for a runtime-typed wave functions calls
/// for a visitor pattern
/// we want:
/// - an extensible list of properties (this means not treating energy special
/// as MPQC did
///   - this also means not too much recompilation when adding a property
/// - an extensible list of wave functions:
///   - not too much work to add a new Wavefunction
///   - ideally no work is function want a new Wavefunction that computes energy
/// - avoid N^2 methods (e.g. CCSD::compute_energy(),
/// CCSD_T::compute_electric_dipole_moment() )
/// - allow computation of derivatives of properties
/// - allow precision tracking

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

////////////////////////////////////////////////////////////////////////

/// \brief Base for classes that provide \c Properties .

/// This provides to the class that inherits this an ability to visit
/// each property \c P in \c Properties by overloading
/// the corresponding \c P::Provider::can_evaluate and \c
/// P::Provider::evaluate methods.
/// @tparam Properties the property type list
template <typename... Properties>
class Provides : public Properties::Provider... {};

/// @return true if Provider can provide Property
template <typename Property, typename Provider>
bool provides(const std::shared_ptr<Provider>& provider) {
  return std::dynamic_pointer_cast<Property::Provider>(provider) == nullptr;
}

namespace detail {
/// has_provider<T>::value is true if T::Provider is a valid type
template <typename T>
struct has_provider {
  template <typename U>
  static std::true_type test(typename U::Provider*);
  template <typename U>
  static std::false_type test(...);
  static constexpr const bool value = decltype(test<T>(nullptr))::value;
};

/// to_pointer<T> converts obj to a pointer:
/// - if T is a pointer type, will return obj
/// - if T is a value or a reference, will take the address
/// - if T is a shared_ptr, will return the corresponding raw pointer
template <typename T>
T* to_pointer(T* obj) {
  return obj;
}
template <typename T>
auto to_pointer(const T& obj)
    -> std::enable_if_t<utility::meta::is_shared_ptr<T>::value,
                        decltype(obj.get())> {
  return obj.get();
}
template <typename T>
typename std::decay<T>::type* to_pointer(T& obj) {
  return &obj;
}

/// obtains a description of the object pointer to by \c obj_ptr
template <typename T>
std::string description(T* obj_ptr) {
  auto dc_obj_ptr = dynamic_cast<DescribedClass*>(obj_ptr);
  if (dc_obj_ptr != nullptr) {
    return std::string("class ") + dc_obj_ptr->class_key();
  } else {
    std::ostringstream oss;
    oss << "object @ " << obj_ptr << " of type " << typeid(*obj_ptr).name();
    return oss.str();
  }
}

}  // namespace detail

/// Evaluates \c property using \c provider
template <typename Property, typename Provider>
std::enable_if_t<detail::has_provider<Property>::value, Property&> operator<<(
    Property& property, Provider& provider) {
  auto provider_ptr = detail::to_pointer(provider);
  auto* evaluator = dynamic_cast<typename Property::Provider*>(provider_ptr);
  if (evaluator == nullptr) {
    std::ostringstream oss;
    oss << detail::description(provider_ptr) << " does not compute "
        << detail::description(&property)
        << ", needs to derive from an appropriate Provides<>";
    throw ProgrammingError(oss.str().c_str(), __FILE__, __LINE__);
  }
  evaluator->evaluate(&property);
  return property;
}

}  // namespace mpqc

#endif  // SRC_MPQC_CHEMISTRY_QC_PROPERTIES_PROPERTY_H_
