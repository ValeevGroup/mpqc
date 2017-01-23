#ifndef SRC_MPQC_CHEMISTRY_QC_PROPERTIES_PROPERTY_H_
#define SRC_MPQC_CHEMISTRY_QC_PROPERTIES_PROPERTY_H_

#include <atomic>
#include <cstdint>
#include <memory>
#include <vector>

#include <TiledArray/type_traits.h>

#include "mpqc/chemistry/molecule/coords.h"
#include "mpqc/chemistry/molecule/molecule.h"
#include "mpqc/chemistry/qc/wfn/wfn.h"
#include "mpqc/util/misc/exception.h"
#include "mpqc/util/misc/task.h"

/// top-level MPQC namespace
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

/// \brief Produces timestamps

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

/// Function maps Parameters to a Value

/// Function keeps parameters and value as part of its state.
/// It recomputes the value as necessary by keeping track of timestamps
/// on parameters and value. The value can also be made obsolete explicitly.
/// Lastly, Function provides support for the Visitor pattern.
///
/// @tparam Value the type of values computed by Function
/// @tparam Parameters the type of parameters used by Function
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

  static const typename Function::value_type get_value(Function* f) {
    return f->get_value();
  }
};

/// N-th order Taylor expansion of a function of \c K variables

/// @tparam Value can be complex valued or a vector (e.g. expansion of a dipole
/// moment)
template <typename Value>
class TaylorExpansionCoefficients {
 public:
  TaylorExpansionCoefficients() : nvars_(0), derivs_() {}

  TaylorExpansionCoefficients(Value real)
      : nvars_(1), derivs_(1, std::vector<Value>(1, real)) {}

  TaylorExpansionCoefficients(size_t nvars, size_t order)
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
          "TaylorExpansionCoefficients::order called, but the object is not "
          "initialized",
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
std::ostream& operator<<(std::ostream& os,
                         const TaylorExpansionCoefficients<Value>& x) {
  x.print(os);
  return os;
}

/// TaylorExpansionFunction<Value> is a Function that returns
/// TaylorExpansionCoefficients<Value>
template <typename Value, typename Parameters>
class TaylorExpansionFunction
    : public Function<TaylorExpansionCoefficients<Value>, Parameters> {
 public:
  typedef Function<TaylorExpansionCoefficients<Value>, Parameters>
      function_base_type;

  /// @return the order of Taylor expansion
  size_t order() const { return precision_.size() - 1; }

  /// @param ord the derivative order
  /// @return the target precision for derivatives order \c ord
  double target_precision(size_t ord) const { return precision_.at(ord); }

 protected:
  /// the default target precision to use if not provided to the constructor
  static constexpr double default_precision_ = 1e-8;

  // clang-format off
  /**
   * @brief auxiliary KeyVal constructor
   *
   * The KeyVal object will be query the following keywords:
   * | KeyWord | Type | Default| Description |
   * |---------|------|--------|-------------|
   * | precision | {array<real> \| real} | {none \| 1e-8} | target precision for {each \| all} derivative orders|
   * | deriv_order | int | 0 | the highest derivative order; only queried if precision is not given or is given as an array |
   *
   * @note that if keyword \c precision is not given then will use the \c default_precision parameter.
   *
   * @param kv the KeyVal object to be queried
   * @param params the Parameters object that this function depends on
   * @param default_precision the default precision to use
   */
  // clang-format on
  TaylorExpansionFunction(const KeyVal& kv, std::shared_ptr<Parameters> params,
                          double default_precision = default_precision_)
      : TaylorExpansionFunction(params, init_precision(kv, default_precision)) {
  }

 private:
  std::vector<double>
      precision_;  // target precision for each order of expansion

  TaylorExpansionFunction(std::shared_ptr<Parameters> params,
                          std::vector<double> taylor_expansion_precision)
      : function_base_type(params), precision_(taylor_expansion_precision) {
    if (taylor_expansion_precision.empty())
      throw ProgrammingError("empty precision list", __FILE__, __LINE__);
  }

  std::vector<double> init_precision(const KeyVal& kv,
                                     double default_precision) {
    if (kv.exists("precision") &&
        kv.count("precision") > 0) {  // given an array of precisions
      return kv.value<std::vector<double>>("precision");
    } else {
      auto deriv_order = kv.value<size_t>("deriv_order", 0);
      auto precision = kv.value<double>("precision", default_precision);
      return std::vector<double>(deriv_order + 1, precision);
    }
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
 * \brief WavefunctionProperty computes a Taylor expansion of a molecular
 * property
 *        using a visiting Wavefunction .
 *
 * This is to be used as a base class for ALL properties of Wavefunction
 * classes.
 */

template <typename Value>
class WavefunctionProperty
    : public TaylorExpansionFunction<Value, MolecularCoordinates>,
      public Property {
 public:
  using base_type = TaylorExpansionFunction<Value, MolecularCoordinates>;
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
  std::shared_ptr<Wavefunction> wfn() const { return wfn_; }

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
      : base_type(kv,
                  (kv.class_ptr<MolecularCoordinates>("coord")
                       ? kv.class_ptr<MolecularCoordinates>("coord")
                       : std::make_shared<CartMolecularCoordinates>(
                             kv.class_ptr<lcao::Wavefunction>("wfn")->atoms())),
                  default_precision) {
    wfn_ = kv.class_ptr<lcao::Wavefunction>("wfn");
    if (wfn_ == nullptr)
      throw InputError(
          "WavefunctionProperty did not receive a Wavefunction object",
          __FILE__, __LINE__, "wfn");
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
/// the corresponding \c P::Evaluator::can_evaluate and \c
/// P::Evaluator::evaluate
/// methods.
/// @tparam Properties the property type list
template <typename... Properties>
class CanEvaluate : public Properties::Evaluator... {};

}  // namespace mpqc

#endif  // SRC_MPQC_CHEMISTRY_QC_PROPERTIES_PROPERTY_H_
