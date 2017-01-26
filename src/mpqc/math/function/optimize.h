/*
 * optimize.h
 *
 *  Created on: Jan 26, 2017
 *      Author: evaleev
 */

#ifndef SRC_MPQC_MATH_FUNCTION_OPTIMIZE_H_
#define SRC_MPQC_MATH_FUNCTION_OPTIMIZE_H_

#include "mpqc/math/function/function.h"
#include "mpqc/util/keyval/keyval.h"

namespace mpqc {
namespace math {


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

}
}


#endif /* SRC_MPQC_MATH_FUNCTION_OPTIMIZE_H_ */
