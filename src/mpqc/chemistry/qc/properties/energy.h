//
// Created by Chong Peng on 1/11/17.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_PROPERTIES_ENERGY_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_PROPERTIES_ENERGY_H_

#include "mpqc/chemistry/qc/properties/property.h"

namespace mpqc {

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

  // clang-format off
  /**
   * KeyVal constructor
   *
   * keywords for this class
   * | KeyWord | Type | Default| Description |
   * |---------|------|--------|-------------|
   * | wfn | Wavefunction | none | the Wavefunction to use to compute energy  |
   * | precision | double | 1.0e-8 | precision for energy property |
   *
   */
  // clang-format on

  Energy(const KeyVal& kv) : Energy(get_wfn(kv), get_precision(kv)) {}

  Energy(Wavefunction* wfn_ptr,
         std::initializer_list<double> taylor_expansion_precision)
      : WavefunctionProperty<double>(wfn_ptr, taylor_expansion_precision) {}

  void evaluate() override;

private:
  Wavefunction* get_wfn(const KeyVal& kv);
  std::initializer_list<double> get_precision(const KeyVal& kv);
};

} // namespace mpqc

#endif //  MPQC4_SRC_MPQC_CHEMISTRY_QC_PROPERTIES_ENERGY_H_
