/*
 * findif.h
 *
 *  Created on: Jan 26, 2017
 *      Author: evaleev
 */

#ifndef SRC_MPQC_MATH_FUNCTION_FINDIF_H_
#define SRC_MPQC_MATH_FUNCTION_FINDIF_H_

#include "mpqc/math/function/taylor.h"
#include "mpqc/util/keyval/keyval.h"

namespace mpqc {
namespace math {

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

}
}



#endif /* SRC_MPQC_MATH_FUNCTION_FINDIF_H_ */
