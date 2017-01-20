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
   *  inherit from Energy::Evaluator
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
   * @brief KeyVal constructor
   *
   * The KeyVal object will be queried for all keywords of the WavefunctionProperty class,
   * as well as the following keywords:
   * | KeyWord | Type | Default| Description |
   * |---------|------|--------|-------------|
   * | wfn | Wavefunction | none | the Wavefunction that will compute this  |
   * | precision | {array<real> \| real} | {none \| 1e-8} | target precision for {each \| all} derivative orders|
   * | deriv_order | int | 0 | the highest derivative order; only queried if precision is not given or is given as an array |
   */
  // clang-format on

  explicit Energy(const KeyVal& kv) : WavefunctionProperty(kv, get_precision(kv)) {}

private:
  std::vector<double> get_precision(const KeyVal& kv);
  void do_evaluate() override;
};

} // namespace mpqc

#endif //  MPQC4_SRC_MPQC_CHEMISTRY_QC_PROPERTIES_ENERGY_H_
