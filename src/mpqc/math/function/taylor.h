/*
 * taylor.h
 *
 *  Created on: Jan 26, 2017
 *      Author: evaleev
 */

#ifndef SRC_MPQC_MATH_FUNCTION_TAYLOR_H_
#define SRC_MPQC_MATH_FUNCTION_TAYLOR_H_

#include "mpqc/math/function/function.h"
#include "mpqc/util/misc/exception.h"

namespace mpqc {
namespace math {

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

  /// useful shorthand for derivs(0)[0]
  /// @ return the function value
  const Value& value() const { return derivs(0)[0]; }
  /// useful shorthand for derivs(1)
  /// @ return the function gradient
  const std::vector<Value>& gradient() const { return derivs(1); }
  /// useful shorthand for derivs(2)
  /// @ return the function gradient
  const std::vector<Value>& hessian() const { return derivs(2); }

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
    os << indent << "value = " << printf("%20.15lf", derivs_.at(0).at(0))
       << std::endl;
    if (order() > 0) {
      typedef Eigen::Matrix<Value, Eigen::Dynamic, 1> Vector;
      Eigen::Map<const Vector> grad_map(&derivs_.at(1)[0], derivs_.at(1).size());
      Vector grad = grad_map;
      os << indent << "gradient = " << grad << std::endl;
    }
    assert(order() < 2);
  }

  void write(KeyVal& kv) const {
    kv.assign("value", this->derivs(0).at(0));
    if (order() > 0) {
        kv.assign("gradient", this->derivs(1));
    }
    assert(order() < 2);
  }

 private:
  /// Number of variables
  size_t nvars_;
  /// values of unique derivatives of each order
  /// # of unique derivs of order \f$ O = n \times (n+1) \times (n+O-1) / (1
  /// \times 2 \times O) \f$
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

}
}
#endif /* SRC_MPQC_MATH_FUNCTION_TAYLOR_H_ */
