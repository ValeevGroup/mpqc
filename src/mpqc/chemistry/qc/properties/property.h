#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_PROPERTIES_FUNCTION_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_PROPERTIES_FUNCTION_H_

#include <atomic>
#include <cstdint>
#include <memory>
#include <vector>

#include <TiledArray/type_traits.h>

#include "mpqc/chemistry/molecule/molecule.h"
#include "mpqc/chemistry/qc/wfn/wfn.h"
#include "mpqc/util/misc/exception.h"
#include "mpqc/util/misc/task.h"

namespace mpqc {

using TiledArray::detail::scalar_type;

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

/// Need to recompute is detected by comparing timestamps of property parameters
/// and value
/// timestamps are produced by this factory method
class TimestampFactory {
 public:
  using timestamp_type = uint64_t;
  static timestamp_type make() {
    static std::atomic<timestamp_type> current_timestamp{0};
    return current_timestamp++;
  }
};

/// Timestampable<T> is a proxy to T that keeps the timestamp of the last
/// modification.
/// Stores T on heap to support the default_initialized state.
template <typename T>
class Timestampable {
 public:
  using timestamp_type = TimestampFactory::timestamp_type;

  Timestampable() : timestamp_(get_timestamp()) {}
  Timestampable(const T& val)
      : value_(std::make_shared<T>(val)), timestamp_(get_timestamp()) {}

  /// copy ctor keeps the timestamp, deep copies value
  Timestampable(const Timestampable& other)
      : timestamp_(other.timestamp_), value_(std::make_shared<T>(T(other))) {}

  /// retrieve the value_ as a conversion to T
  operator T() const {
    assert(value_ != nullptr);
    return *value_;
  }

  /// @return the current timestamp
  const timestamp_type& timestamp() const { return timestamp_; }

 private:
  std::shared_ptr<T> value_;
  timestamp_type timestamp_;

  static timestamp_type get_timestamp() { return TimestampFactory::make(); }
};

template <typename T>
std::ostream& operator<<(std::ostream& os, const Timestampable<T>& x) {
  os << T(x);
  return os;
}

template <typename Function>
class FunctionVisitorBase;

/// property = function of a sequence of Parameters yielding a Value
template <typename Value, typename Parameters>
class Function {
 public:
  typedef Value value_type;
  typedef Parameters parameters_type;

  Function() = default;

  Function(const std::shared_ptr<Parameters>& params)
      : value_(), params_(params), obsolete_(true) {}

  /// Timestampable<Value> will be convert to temporary Value,
  /// has to return by value here
  const Value value() {
    if (must_compute()) compute();
    return value_;
  }

  // if value has been computed less recently than the last update of
  // parameters
  // or obsolete
  bool must_compute() const {
    if (value_.timestamp() < params_.timestamp() || obsolete_)
      return true;
    else
      return false;
  }

  std::shared_ptr<Parameters> params() const { return params_; }

 protected:
  /// Direct access to the value of this function, bypasses timestamp check
  /// @return the current value, i.e. \c value_
  const Timestampable<Value>& get_value() const { return value_; }
  /// Sets the value of this function, used by the compute() method of derived
  /// classes
  /// or by Function visitors via FunctionVisitorBase::set_value() .
  /// @param v the value to be returned by Function::value()
  void set_value(Value v) { value_ = Timestampable<Value>(v); }

  /// evaluates \c value , implemented by the derived class
  virtual void compute() = 0;

  // allow classes derived from FunctionVisitorBase<Function> to set the value
  friend class FunctionVisitorBase<Function<Value, Parameters>>;

 private:
  Timestampable<Value> value_;
  Timestampable<std::shared_ptr<Parameters>> params_;
  bool obsolete_;
};

/// FunctionVisitorBase makes possible for Visitors of a child of Function
/// to call Function::set_value
template <typename Function>
class FunctionVisitorBase {
 protected:
  static void set_value(Function* f,
                        const typename Function::value_type& value) {
    f->set_value(value);
  }
};

/// molecular coordinates + taylor expansion params
/// @todo use MolecularCoordinates in MolecularTaylorExpansionParams
class MolecularTaylorExpansionParams {
 public:
  MolecularTaylorExpansionParams(
      const std::shared_ptr<Molecule>& molecule,
      std::vector<double> taylor_expansion_precision)
      : molecule_(molecule),
        precision_(taylor_expansion_precision.begin(),
                   taylor_expansion_precision.end()) {}
  // # of molecular coords + 1 precision * number of derivatives
  size_t nparams() const {
    return molecule_->atoms().size() * 3 + precision_.size();
  }

  /// @return the order of Taylor expansion
  size_t order() const { return precision_.size() - 1; }
  /// @param ord the derivatiove order
  /// @return the target precision for derivatives order \c ord
  double target_precision(size_t ord) const { return precision_.at(ord); }

 private:
  std::shared_ptr<Molecule> molecule_;
  std::vector<double>
      precision_;  // target precision for each order of expansion
};

/// N-th order Taylor expansion of a function of \c K variables

/// @tparam Value can be complex valued or a vector (e.g. expansion of a dipole
/// moment)
template <typename Value>
class TaylorExpansion {
 public:
  TaylorExpansion() : nvars_(0), derivs_() {}

  TaylorExpansion(Value real)
      : nvars_(1), derivs_(1, std::vector<Value>(1, real)) {}

  TaylorExpansion(size_t nvars, size_t order)
      : nvars_(nvars), derivs_(order + 1) {
    derivs_[0].resize(1, Value(0));
    for (int d = 1; d <= order; ++d) {
      derivs_[d].resize((derivs_[d - 1].size() * (nvars_ + d - 1)) / d,
                        Value(0));
    }
  }

  /// access a packed set of unique derivatives of order \c N. Derivative \f$
  /// {i_1, i_2, i_3, ... i_N} \f$
  /// where \f$ i_1 \leq i_2 \leq i_3 ... \leq i_N \f$
  /// is stored in position TODO derive the general formula
  const std::vector<Value>& derivs(size_t N) const { return derivs_.at(N); }
  std::vector<Value>& derivs(size_t N) { return derivs_.at(N); }

  /// @return the order of the expansion, i.e. the maximum derivative order
  /// @throw ProgrammingError if this object is not initialized
  size_t order() const {
    if (derivs_.empty())
      throw ProgrammingError(
          "TaylorExpansion::order called, but the object is not initialized",
          __FILE__, __LINE__);
    return derivs_.size() - 1;
  }

  /// Print to an output stream
  /// @param os the output stream
  /// @throw ProgrammingError if this object is not initialized
  void print(std::ostream& os) const {
    if (derivs_.empty())
      throw ProgrammingError(
          "TaylorExpansion::print called, but the object is not initialized",
          __FILE__, __LINE__);
    os << indent << "value: " << printf("%20.15lf", derivs_.at(0).at(0))
       << std::endl;
    assert(order() == 0);
  }

 private:
  /// Number of variables
  size_t nvars_;
  /// values of unique derivatives of each order
  /// # of unique derivs of order \f$ O = n \times (n+1) \times (n+O-1) / (1
  /// \times 2 \times O)
  std::vector<std::vector<Value>> derivs_;
};

template <typename Value>
std::ostream& operator<<(std::ostream& os, const TaylorExpansion<Value>& x) {
  x.print(os);
  return os;
}

/// molecular taylor expansion = local expansion of a function of molecular
/// coordinates
/// NB does not depend on fields, need MolecularPropertyInEMField for
/// something like that
template <typename Value>
class MolecularTaylorExpansion
    : public Function<TaylorExpansion<Value>, MolecularTaylorExpansionParams> {
 public:
  typedef Function<TaylorExpansion<Value>, MolecularTaylorExpansionParams>
      function_base_type;

  // make_molecular_params needs to create a MolecularCoordinates object and
  // set callbacks so that Molecule can update its timestamp when
  // Molecule::set_coordinates is called

  MolecularTaylorExpansion(
      const std::shared_ptr<Molecule>& molecule_,
      std::vector<typename scalar_type<Value>::type>
          taylor_expansion_precision)
      : Function<TaylorExpansion<Value>, MolecularTaylorExpansionParams>(
            std::make_shared<MolecularTaylorExpansionParams>(
                molecule_, taylor_expansion_precision)) {}

  /// @return the order of Taylor expansion
  size_t order() const {
    return std::static_pointer_cast<MolecularTaylorExpansionParams>(this->params())
        ->order();
  }

  /// @param ord the derivatiove order
  /// @return the target precision for derivatives order \c ord
  double target_precision(size_t ord) const {
    return std::static_pointer_cast<MolecularTaylorExpansionParams>(this->params())
        ->target_precision(ord);
  }
};

/// this is the base for all properties that MPQC can compute via the input.
/// MPQC main will read KeyVal and search for a Property object, compute it
/// using the given wave function
class Property : public Task {
 public:
  // in a visitor pattern this is the "accept" method
  // the argument, Wavefunction*, does not appear here, it will be a member
  // of the derived class
  virtual void evaluate() = 0;

 private:
  /// implements Task::execute()
  void execute() override final { evaluate(); }
};

/**
 *
 * wave function property = molecular property computable from a wave function
 * base class for ALL properties of wave functions PropertyBase
 */

template <typename Value>
class WavefunctionProperty : public MolecularTaylorExpansion<Value>,
                             public Property {
 public:
  using typename MolecularTaylorExpansion<Value>::value_type;
  using typename MolecularTaylorExpansion<Value>::function_base_type;
  WavefunctionProperty(Wavefunction* wfn_ptr,
                       std::vector<typename scalar_type<Value>::type>
                           taylor_expansion_precision)
      : MolecularTaylorExpansion<Value>(wfn_ptr->atoms(),
                                        taylor_expansion_precision),
        wfn_(wfn_ptr) {
    if (wfn_ptr == nullptr)
      throw ProgrammingError("WavefunctionProperty ctor received null wfn ptr",
                             __FILE__, __LINE__);
  }

 protected:
  Wavefunction* wfn() { return wfn_; }

  virtual void do_evaluate() = 0;

 private:
  Wavefunction* wfn_;

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

    // report the result
    ExEnv::out0() << indent << "Property \"" << this->class_key()
                  << "\" computed with Wavefunction \"" << wfn_->class_key()
                  << "\":" << std::endl
                  << incindent << this->get_value() << std::endl
                  << decindent;
  }
};

////////////////////////////////////////////////////////////////////////

/// \brief Base for classes that can evaluate \c Properties .

/// This provides to the class that inherits this an ability to visit
/// each property \c P in \c Properties by overloading
/// the corresponding \c P::Evaluator::can_evaluate and \c P::Evaluator::evaluate
/// methods.
/// @tparam Properties the property type list
template <typename... Properties>
class CanEvaluate;

template <typename Property0, typename... RestOfProperties>
class CanEvaluate<Property0, RestOfProperties...>
    : public Property0::Evaluator, public CanEvaluate<RestOfProperties...> {};

template <typename Property>
class CanEvaluate<Property> : public Property::Evaluator {};

}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_PROPERTIES_ENERGY_H_
