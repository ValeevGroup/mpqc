#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_PROPERTIES_FUNCTION_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_PROPERTIES_FUNCTION_H_

#include <atomic>
#include <cstdint>
#include <memory>
#include <vector>

#include <TiledArray/type_traits.h>

#include <mpqc/chemistry/molecule/molecule.h>
#include <mpqc/util/misc/exception.h>

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
    std::atomic<timestamp_type> current_timestamp{timestamp_type{0}};
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
  Timestampable(const Timestampable& other) : timestamp_(other.timestamp_),
      value_(std::make_shared<T>(T(other))) {}

  /// retrieve the value_ as a conversion to T
  operator T() {
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

/// property = function of a sequence of Parameters yielding a Value
template <typename Value, typename Parameters>
class Function {
 public:
  typedef Value value_type;
  typedef Parameters parameters_type;

  Function(const std::shared_ptr<Parameters>& params)
      : params_(params), obsolete_(true) {}

  const Value& value() {
    if (must_compute()) compute();
    return value_;
  }

  // if value has been computed less recently than the last update of
  // parameters
  // or obsolete
  bool must_compute() const {
    if (value_->timestamp() < params_->timestamp() || obsolete_)
      return true;
    else
      return false;
  }

  const std::shared_ptr<Parameters>& params() const { return params_; }

 protected:
  /// implement in derived class
  virtual void compute() = 0;

  const Timestampable<Value>& get_value() const { return value_; }
  void set_value(Value v) { value_ = v; }

 private:
  Timestampable<Value> value_;
  std::shared_ptr<Parameters> params_;
  bool obsolete_;
};

/// molecular coordinates + taylor expansion params
class MolecularTaylorExpansionParams {
 public:
  MolecularTaylorExpansionParams(
      const std::shared_ptr<Molecule>& molecule,
      std::initializer_list<double> taylor_expansion_precision)
      : molecule_(molecule),
        precision_(taylor_expansion_precision.begin(),
                   taylor_expansion_precision.end()) {}
  // 3 coords * # of atoms + 1 precision * number of derivatives
  size_t nparams() const { return molecule_->atoms().size() * 3 + precision_.size(); }

  int order() const { return precision_.size() - 1; }
  double desired_precision(int order) const { return precision_.at(order); }

 private:
  std::shared_ptr<Molecule> molecule_;
  std::vector<double>
      precision_;  // precision desired for each order of expansion
};

/// N-th order Taylor expansion of a function of \c K variables

/// @tparam Real the real type
template <typename Real>
class TaylorExpansion {
 public:
  TaylorExpansion() : nvars_(0), order_(0), derivs_(order_+1, std::vector<Real>(1, Real(0))) {}
  TaylorExpansion(size_t nvars, int order) : nvars_(nvars), order_(order), derivs_(order_+1) {
    derivs_[0].resize(1, Real(0));
    for(int d=1; d <= order_; ++d) {
      derivs_[d].resize( (derivs_[d-1].size() * (nvars_ + d - 1)) / d , Real(0));
    }
  }

  /// access a packed set of unique derivatives of order \c N. Derivative \f$ {i_1, i_2, i_3, ... i_N} \f$
  /// where \f$ i_1 \leq i_2 \leq i_3 ... \leq i_N \f$
  /// is stored in position TODO derive the general formula
  const std::vector<Real> & derivs(int N) const { return derivs_.at(N); }
  std::vector<Real>& derivs(int N) { return derivs_.at(N); }

 private:
  /// Number of variables
  size_t nvars_;
  /// expansion order
  int order_;
  /// values of unique derivatives of each order
  /// # of unique derivs of order \f$ O = n \times (n+1) \times (n+O-1) / (1 \times 2 \times O)
  std::vector<std::vector<Real>> derivs_;
};

/// molecular taylor expansion = local expansion of a function of molecular coordinates
/// NB does not depend on fields, need MolecularPropertyInEMField for
/// something like that
template <typename Value>
class MolecularTaylorExpansion
    : public Function<TaylorExpansion<Value>, MolecularTaylorExpansionParams> {
 public:
  // make_molecular_params needs to create a MolecularCoordinates object and
  // set
  // callbacks
  // so that Molecule can update its timestamp when Molecule::set_coordinates
  // is
  // called
  MolecularTaylorExpansion(
      const std::shared_ptr<Molecule>& molecule_,
      std::initializer_list<typename scalar_type<Value>::type> taylor_expansion_precision)
      : Function<TaylorExpansion<Value>, MolecularTaylorExpansionParams>(
            std::make_shared<MolecularTaylorExpansionParams>(
                molecule_, taylor_expansion_precision)) {}
};

/// this is the base for all properties that MPQC can compute via the input.
/// MPQC main will read KeyVal and search for a PropertyBase object, compute it
/// using the given wave function
class PropertyBase : public DescribedClass {
 public:
  // in a visitor pattern this is the "accept" method
  // the argument, Wavefunction*, does not appear here, it will be a member
  // of the is derived class
  virtual void evaluate() = 0;
};


/// wave function property = molecular property computable from a wave function
/// base class for ALL properties of wave functions
/// PropertyBase
template <typename Value>
class WavefunctionProperty : public MolecularTaylorExpansion<Value>, public PropertyBase {
 public:
  using typename MolecularTaylorExpansion<Value>::value_type;
  WavefunctionProperty(
      Wavefunction* wfn_ptr,
      std::initializer_list<typename scalar_type<Value>::type> taylor_expansion_precision)
      : MolecularTaylorExpansion<Value>(wfn_ptr->atoms(),
                                 taylor_expansion_precision),
        wfn_(wfn_ptr) {
    if (wfn_ptr == nullptr)
      throw ProgrammingError("WavefunctionProperty ctor received null wfn ptr",
                             __FILE__, __LINE__);
  }

 protected:
  Wavefunction* wfn() { return wfn_; }

 private:
  Wavefunction* wfn_;

  void compute() {
    auto original_value_timestamp = this->get_value().timestamp();
    this->evaluate();
    auto updated_value_timestamp = this->get_value().timestamp();
    if (original_value_timestamp == updated_value_timestamp)
      throw ProgrammingError(
          "WavefunctionProperty::compute() forgot to call set_value?", __FILE__,
          __LINE__);
  }
};

////////////////////////////////////////////////////////////////////////
// to add a wavefunction property class P:
// - derive class P from WavefunctionProperty<T> and override
// P::evaluate()
// - define class P::EvaluatorBase to be used as a public base classes that
// can compute it

class Energy : public WavefunctionProperty<double> {
 public:
  /// every class that can evaluate Energy (e.g. Wavefunction) will publicly
  /// inherit from Energy::EvaluatorBase
  class EvaluatorBase {
   public:
    /// EvaluatorBase::can_evaluate returns true if \c energy can be computed.
    /// For example, if \c energy demands taylor expansion to 1st order
    /// but this wave function does not have analytic nuclear gradients,
    /// will return false.
    virtual bool can_evaluate(Energy* energy) = 0;
    /// EvaluatorBase::evaluate computes the taylor expansion of the energy
    /// and uses set_value to assign the values to \c energy
    virtual void evaluate(Energy* energy) = 0;
  };

  friend class EvaluatorBase;  // hence EnergyEvaluator can call
                               // Function::set_value()

  Energy(Wavefunction* wfn_ptr,
         std::initializer_list<double> taylor_expansion_precision)
      : WavefunctionProperty<double>(wfn_ptr, taylor_expansion_precision) {}

  void evaluate() override {
    auto evaluator = dynamic_cast<EvaluatorBase*>(wfn());
    if (evaluator == nullptr) {
      std::ostringstream oss;
      // TODO Must implement DescribedClass::key() instead of using RTTI's
      // unpredicable output
      oss << "Wavefunction " << typeid(*wfn()).name()
          << " cannot compute Energy" << std::endl;
      throw InputError(oss.str().c_str(), __FILE__, __LINE__);
    }
    evaluator->evaluate(this);
  }
};

#if 0
  namespace lcao {

  // there is a way to automate making lists of properties
  // that a given Wavefunction supports
  class CCSD : public Wavefunction, public Energy::EvaluatorBase {
   public:
    bool can_evaluate(Energy* energy) override {
    }
    void evaluate(Energy* energy) override {
    }
  };

  // there is a way to automate making lists of properties
  // that a given Wavefunction supports
  class GF2F12 : public Wavefunction, public GFRealPole::EvaluatorBase {
   public:
    bool can_evaluate(GFRealPole* pole) override {
    }
    void evaluate(GFRealPole* pole) override {
    }
  };
  }  // namespace lcao

  namespace mra {

  class SCF : public Wavefunction, public Energy::EvaluatorBase {
   public:
    bool can_evaluate(Energy* energy) override {
    }
    void evaluate(Energy* energy) override {
    }
  };
  }  // namespace mra
#endif

}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_PROPERTIES_ENERGY_H_
