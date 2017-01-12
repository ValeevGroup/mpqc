//
// Created by Chong Peng on 1/11/17.
//

#ifndef MPQC_ENERGY_H
#define MPQC_ENERGY_H

#include "property.h"

namespace mpqc{

/**
 *
 * to add a wavefunction property class P:
 * - derive class P from WavefunctionProperty<T> and override
 * P::evaluate()
 * - define class P::EvaluatorBase to be used as a public base classes that
 * can compute it
 */
class Energy : public WavefunctionProperty<double> {
public:
  using typename WavefunctionProperty<double>::function_base_type;

  /**
   *  every class that can evaluate Energy (e.g. Wavefunction) will publicly
   *  inherit from Energy::EvaluatorBase
   */
  class Evaluator : public FunctionVisitorBase<function_base_type> {
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

  Energy(Wavefunction* wfn_ptr,
         std::initializer_list<double> taylor_expansion_precision)
      : WavefunctionProperty<double>(wfn_ptr, taylor_expansion_precision) {}

  void evaluate() override;
};

}

#endif //MPQC_ENERGY_H
