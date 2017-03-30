#ifndef SRC_MPQC_MATH_FUNCTION_FUNCTION_H_
#define SRC_MPQC_MATH_FUNCTION_FUNCTION_H_

#include "mpqc/util/misc/timestamp.h"

namespace mpqc {
namespace math {

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
  template <typename X> using Timestampable = ::mpqc::utility::Timestampable<X>;

  Function() = default;
  virtual ~Function() { }

  explicit Function(std::shared_ptr<Parameters> params)
      : value_(), params_(params), obsolete_(true) {}

  /// Timestampable<Value> will be converted to temporary Value,
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

}
}  // namespace mpqc

#endif /* SRC_MPQC_MATH_FUNCTION_FUNCTION_H_ */
