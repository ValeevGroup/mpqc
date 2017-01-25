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
  explicit Timestampable(T&& val)
      : value_(std::make_shared<T>(val)), timestamp_(get_timestamp()) {}
  explicit Timestampable(std::shared_ptr<T> val)
      : value_(val), timestamp_(get_timestamp()) {}

  /// copy ctor keeps the timestamp, deep copies value
  Timestampable(const Timestampable& other)
      : timestamp_(other.timestamp_), value_(std::make_shared<T>(T(other))) {}

  Timestampable& operator=(const T& val) {
    value_ = std::make_shared<T>(val);
    timestamp_ = get_timestamp();
    return *this;
  }
  Timestampable& operator=(std::shared_ptr<T> val) {
    value_ = val;
    timestamp_ = get_timestamp();
    return *this;
  }

  /// retrieve a (non-const) reference value_ updates the timestamp
  operator T&() {
    assert(value_ != nullptr);
    timestamp_ = get_timestamp();
    return *value_;
  }

  /// retrieve a const reference value_
  operator const T&() const {
    assert(value_ != nullptr);
    return *value_;
  }

  /// retrieve a shared pointer, updates the timestamp
  operator std::shared_ptr<T>() {
    assert(value_ != nullptr);
    timestamp_ = get_timestamp();
    return value_;
  }

  /// retrieve a const shared pointer
  operator std::shared_ptr<const T>() const {
    assert(value_ != nullptr);
    return value_;
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
  os << static_cast<const T&>(x);
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
/// @note both Value and Parameters are managed via shared pointers
///
/// @tparam Value the type of value computed by Function; can be an abstract
/// class
/// @tparam Parameters the type of parameters used by Function; can be an
/// abstract class
template <typename Value, typename Parameters>
class Function {
 public:
  typedef Value value_type;
  typedef Parameters parameters_type;

  Function() = default;
  virtual ~Function() = default;

  explicit Function(std::shared_ptr<Parameters> params)
      : value_(), params_(params), obsolete_(true) {}

  /// Timestampable<Value> will be convert to temporary Value,
  /// has to return by value here
  std::shared_ptr<const Value> value() {
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

  std::shared_ptr<const Parameters> params() const { return params_; }

  virtual void set_params(std::shared_ptr<Parameters> params) {
    params_ = params;
  }

 protected:
  /// Direct access to the value of this function, bypasses timestamp check
  /// @return the current value, i.e. \c value_
  const Timestampable<Value>& get_value() const { return value_; }
  /// Sets the value of this function, used by the compute() method of derived
  /// classes
  /// or by Function visitors via FunctionVisitorBase::set_value() .
  /// @param v the value to be returned by Function::value()
  void set_value(Value v) { value_ = Timestampable<Value>(std::move(v)); }

  /// evaluates \c value , implemented by the derived class
  virtual void compute() = 0;

  // allow classes derived from FunctionVisitorBase<Function> to set the value
  friend class FunctionVisitorBase<Function<Value, Parameters>>;

 private:
  Timestampable<Value> value_;
  Timestampable<Parameters> params_;
  bool obsolete_;
};

/// FunctionVisitorBase makes possible for Visitors of a child of Function
/// to call Function::set_value
template <typename Function>
class FunctionVisitorBase {
 protected:
  static void set_value(Function* f,
                        typename Function::value_type value) {
    f->set_value(std::move(value));
  }

  static const typename Function::value_type& get_value(Function* f) {
    return f->get_value();
  }
};

namespace detail {
namespace function {

template <typename T,
          typename = typename std::enable_if<!std::is_abstract<T>::value>::type>
std::shared_ptr<typename std::decay<T>::type> clone(T* other) {
  return std::make_shared<typename std::decay<T>::type>(*other);
}

template <typename T,
          typename = typename std::enable_if<!std::is_abstract<T>::value>::type>
std::shared_ptr<typename std::decay<T>::type> clone(std::shared_ptr<T> other) {
  return std::make_shared<typename std::decay<T>::type>(*other);
}

template <typename T>
std::shared_ptr<typename std::decay<T>::type> clone(T* other) {
  return other->clone();
}

template <typename T>
std::shared_ptr<typename std::decay<T>::type> clone(std::shared_ptr<T> other) {
  return other->clone();
}

}  // namespace function
}  // namespace detail

/// N-th order Taylor expansion of a function of \c K variables

/// @tparam Value can be complex valued or a vector (e.g. expansion of a dipole
/// moment)
template <typename Value>
class TaylorExpansionCoefficients {
 public:
  TaylorExpansionCoefficients() : nvars_(0), derivs_() {}

  TaylorExpansionCoefficients(Value val)
      : nvars_(0), derivs_(1, std::vector<Value>(1, val)) {}

  TaylorExpansionCoefficients(std::vector<std::vector<Value>> coefs)
      : nvars_(coefs.size() > 0 ? coefs[1].size() : 0), derivs_(std::move(coefs)) {}

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

  /// evaluates this object
  virtual void evaluate() = 0;

  /// prints this object to \c os
  /// @param os the output stream
  virtual void print(std::ostream& os = ExEnv::out0()) const = 0;

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

  void print(std::ostream& os = ExEnv::out0()) const override {
    os << indent << "Property \"" << this->class_key() << "\":\n" << incindent;
    os << indent << "wfn:\n" << incindent;
    wfn_->print(os);
    os << decindent;
    os << "value = " << this->get_value() << std::endl;
    os << decindent;
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

/// \brief Base for classes that can evaluate \c Properties .

/// This provides to the class that inherits this an ability to visit
/// each property \c P in \c Properties by overloading
/// the corresponding \c P::Evaluator::can_evaluate and \c
/// P::Evaluator::evaluate methods.
/// @tparam Properties the property type list
template <typename... Properties>
class CanEvaluate : public Properties::Evaluator... {};

////////////////////////////////////////////////////////////////////////

/// Optimizer<Real, Params> seeks stationary points of
/// TaylorExpansionFunction<Real, Params>

/// It is itself a Function, namely Function<Params, std::tuple<>>.
template <typename Real, typename Params>
class Optimizer : public Function<Params, std::tuple<>>,
                  virtual public DescribedClass {
 public:
  typedef TaylorExpansionFunction<Real, Params> function_type;

  Optimizer(const KeyVal& kv) {
    // obtain target precision
    precision_ = kv.value<double>("precision", 1e-4);

    // obtain target precision
    function_ = kv.class_ptr<function_type, std::true_type>("function");
  }

  std::shared_ptr<function_type> function() const { return function_; }
  double precision() const { return precision_; }

 private:
  std::shared_ptr<function_type> function_;
  double precision_;
};

/// It is itself a Function, namely Function<Params, OptimizerParams>.
template <typename Real, typename Params>
class QuasiNewtonOptimizer : public Optimizer<Real, Params> {
 public:
  using base_type = Optimizer<Real, Params>;
  using typename base_type::function_type;
  using base_type::function;

  QuasiNewtonOptimizer(const KeyVal& kv) : Optimizer<Real, Params>(kv) {
    // if function cannot compute gradients, look for numerical gradient
    if (this->function()->order() < 1) {
      gradient_ = kv.class_ptr<function_type>("gradient");
      if (!gradient_)
        throw InputError(
            "QuasiNewtonOptimizer: function cannot compute gradients, and "
            "gradient object not provided",
            __FILE__, __LINE__, "gradient");
      if (gradient_->order() < 1)
        throw InputError(
            "QuasiNewtonOptimizer: the gradient object cannot compute "
            "gradients",
            __FILE__, __LINE__, "gradient");
    }

    if (this->function()->order() < 2) {
      guess_hessian_ = kv.class_ptr<function_type>("guess_hessian");
      if (!guess_hessian_)
        throw InputError(
            "QuasiNewtonOptimizer: function cannot compute hessians, and "
            "guess_hessian object not provided",
            __FILE__, __LINE__, "guess_hessian");
    }

    // hessian_update_ = kv.class_ptr<HessianUpdate>("hessian_update");

    max_iter_ = kv.value<size_t>("max_iter", 40);
  }

 private:
  size_t max_iter_;
  std::shared_ptr<function_type> gradient_;
  std::shared_ptr<function_type> guess_hessian_;
  // std::shared_ptr<HessianUpdate> hessian_update_;

  using Vector = Eigen::Matrix<Real, Eigen::Dynamic, 1>;
  using Matrix = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>;

  void compute() {
    bool converged = false;

    using TA::detail::size;
    const auto nparams = size(*(this->function()->params()));

    guess_hessian_->set_params(this->function()->params());
    auto hess_vec = guess_hessian_->value().derivs(2);
    // convert hessian to Eigen Matrix
    assert(hess_vec.size() == nparams * (nparams + 1) / 2);
    Matrix hess(nparams, nparams);
    for (auto r = 0, rc = 0; r != nparams; ++r) {
      for (auto c = 0; c != r; ++c, ++rc) {
        hess(r, c) = hess(c, r) = hess_vec[rc];
      }
    }

    do {
      auto current_params = this->function()->params();
      if (gradient_) gradient_->set_params(current_params);

      // compute value and gradient
      const auto& val_and_grad =
          gradient_.nonnull() ? gradient_->value() : this->function()->value();

      // extract gradient
      auto grad_vec = val_and_grad.derivs(1);
      Vector grad = Eigen::Map<Vector>(&grad_vec[0], grad_vec.size());

      // check convergence
      converged = grad.norm() < this->precision();

      if (!converged) {
        // compute step
        Vector step = -hess.inverse() * grad;

        // update params
        // this->function()->set_params(this->function()->params() + step);
        assert(false && "not yet implemented");

        // update hessian
        // hessian_update_->update(hess, grad);
      }
    } while (!converged);
  }
};

/// Computes finite-difference approximation to function derivatives
template <size_t Order, typename Value, typename Parameters>
class FiniteDifferenceDerivative
    : public TaylorExpansionFunction<Value, Parameters>,
      virtual public DescribedClass {
 public:
  static_assert(Order == 1,
                "Only 1st-order FiniteDifferenceDerivative implemented");

  FiniteDifferenceDerivative(const KeyVal& kv)
      :   FiniteDifferenceDerivative(kv,
                                     nullptr,
                                     1e-6)
  {}

 protected:
  // tuple of coordinate indices, in decreasing order, specifies each derivative
  using DerivIdx = std::array<size_t, Order>;
  // each displacement = set of displacement values for each of N coordinates in
  // the derivative
  using Displacement = std::array<double, Order>;
  // sequence of {displacement,coefficient} pairs, with coefficients including
  // deltas
  using LinearCombinationOfDisplacements =
      std::vector<std::pair<Displacement, double>>;

  FiniteDifferenceDerivative(const KeyVal& kv,
                             std::shared_ptr<Parameters> params,
                             double default_target_precision)
      : TaylorExpansionFunction<Value, Parameters>(kv, params, default_target_precision),
        delta_(kv.value<double>("delta", 1e-2)),
        error_order_(kv.value<size_t>("error_order", 0)) {
    function_ = kv.class_ptr<function_type, std::true_type>("function");
    if (function_ == nullptr)
      throw InputError(
          "FiniteDifferenceDerivative was not given a Function to "
          "differentiate",
          __FILE__, __LINE__, "function");
  }

  /// derived classes may need to customize how displacements are computed
  /// @param idx the index of the derivative
  /// @return the stensil formula as a LinearCombinationOfDisplacements
  virtual LinearCombinationOfDisplacements generate_displacements(
      DerivIdx idx) const {
    /// manual implementations of
    /// FiniteDifferenceDerivative::generate_displacements for different orders
    LinearCombinationOfDisplacements result;
    switch (Order) {
      case 1: {
        const auto one_over_delta_ = 1 / delta_;
        result.emplace_back(
            std::make_pair(Displacement{{delta_}}, one_over_delta_));
        result.emplace_back(
            std::make_pair(Displacement{{-delta_}}, -one_over_delta_));
      } break;
      default:
        assert(false && "unreachable");
    }
    return result;
  }

  /// the function to differentiate
  using function_type = TaylorExpansionFunction<Value, Parameters>;
  std::shared_ptr<const function_type> function() const { return function_; }
  double delta() const { return delta_; }
  size_t error_order() const { return error_order_; }

 private:
  std::shared_ptr<function_type> function_;
  /// the list of displacements for each coordinate
  std::map<DerivIdx, LinearCombinationOfDisplacements> disps_;
  /// the displacement size
  double delta_;
  /// the accuracy of the highest-order derivative in terms of powers of \c
  /// delta_
  /// \c error_order=0 means derivative is accurate to \c delta_^2 ,
  /// \c error_order=1 -- \c \delta^4 , etc.
  size_t error_order_;

  const std::shared_ptr<function_type>& function() { return function_; }

  void compute() override {
    // grab ref parameters to reset later
    using detail::function::clone;
    auto function_ref_params = clone(function()->params());

    auto ref_params = this->params();
    using TA::detail::size;
    const auto nparams = size(*ref_params);

    /////////////////////////
    // create displacements
    /////////////////////////

    // loop over unique derivatives
    for (size_t p = 0; p != nparams; ++p) {
      // for each derivative generate displacements
      if (error_order_ > 0) {
        throw FeatureNotImplemented(
            "FiniteDifferenceDerivative: only leading order stensil "
            "implemented yet",
            __FILE__, __LINE__);
      }

      auto idx = DerivIdx{{p}};
      disps_.insert(std::make_pair(idx, generate_displacements(idx)));
    }

    /////////////////////////
    // compute
    /////////////////////////
    size_t coord = 0;
    std::vector<Value> grad_vec(nparams);
    for (const auto& disp : disps_) {
      Value result = 0;
      const auto& idx = disp.first;
      ExEnv::out0() << indent << "displacement " << disp.first[0] << std::endl;
      using detail::function::clone;
      auto disp_params = clone(ref_params);
      for (const auto& term : disp.second) {
        // apply the displacement
        increment(disp_params.get(), idx, term.first);
        function()->set_params(disp_params);
        // compute
        auto val = function()->value();
        result += term.second * val->derivs(0)[0];  // have energies only
        // revert the displacement
        decrement(disp_params.get(), idx, term.first);
      }

      grad_vec[coord] = result;  // computing gradient only
      ++coord;
    }

    // reset the parameters
    function()->set_params(std::const_pointer_cast<Parameters>(ref_params));

    // finalize: set_value
    // N.B. return 0 for energy
    std::vector<std::vector<Value>> energy_and_grad(2);
    energy_and_grad[0].push_back(static_cast<Value>(0));
    energy_and_grad[1] = std::move(grad_vec);
    TaylorExpansionCoefficients<Value> value(energy_and_grad);
    this->set_value(value);
  }

};

/// Evaluates FiniteDifferenceDerivative with respect to molecular coordinates
template <size_t Order, typename Value>
class MolecularFiniteDifferenceDerivative
    : public FiniteDifferenceDerivative<Order, Value, MolecularCoordinates>,
      public Property {
 public:
  MolecularFiniteDifferenceDerivative(const KeyVal& kv)
      : FiniteDifferenceDerivative<Order, Value, MolecularCoordinates>(
            kv, (kv.class_ptr<MolecularCoordinates>("coords")
                     ? kv.class_ptr<MolecularCoordinates>("coords")
                     : std::make_shared<CartMolecularCoordinates>(
                           kv.class_ptr<Molecule>("molecule"))),
                           default_target_precision_) {}

  constexpr static double default_target_precision_ = 1e-6;

  void print(std::ostream& os = ExEnv::out0()) const override {
    os << indent << "Property " << this->class_key() << ":" << std::endl << incindent;
    auto func = this->function();
    // if function is a Property, use its print method
    auto func_prop = std::dynamic_pointer_cast<const Property>(func);
    if (func_prop) {
      os << indent << "function:\n" << incindent;
      func_prop->print(os);
      os << decindent;
    }
    else {
      os << indent << "function = unknown type\n";
    }
    os << indent << "delta = " << this->delta() << std::endl;
    os << indent << "error_order = " << this->error_order() << std::endl;
    os << indent << "value = " << this->get_value() << std::endl;
    os << decindent;
  }

 private:
  /// overrides Property::evaluate()
  void evaluate() override { this->value(); }
};

}  // namespace mpqc

#endif  // SRC_MPQC_CHEMISTRY_QC_PROPERTIES_PROPERTY_H_
