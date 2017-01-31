//
// Created by Chong Peng on 1/11/17.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_PROPERTIES_ENERGY_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_PROPERTIES_ENERGY_H_

#include "mpqc/chemistry/qc/properties/property.h"
#include "mpqc/math/function/optimize.h"

namespace mpqc {

/*
* to add a wavefunction property class P:
* - derive class P from WavefunctionProperty<T> and override
* P::evaluate()
* - define class P::EvaluatorBase to be used as a public base classes that
* can compute it
*/

/// Taylor expansion of the molecular energy computed by a Wavefunction
class Energy : public WavefunctionProperty<double> {
public:
  using typename WavefunctionProperty<double>::function_base_type;

  /**
   *  every class that can evaluate Energy (e.g. Wavefunction) will publicly
   *  inherit from Energy::Evaluator
   *
   *  @sa CanEvaluate
   */
  class Evaluator : public math::FunctionVisitorBase<function_base_type> {
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
   * @brief The KeyVal constructor
   * @param kv the KeyVal object to be queried
   *
   * \c kv will be queried for all keywords of the WavefunctionProperty class,
   * as well as the following keywords:
   * | KeyWord | Type | Default| Description |
   * |---------|------|--------|-------------|
   * | wfn | Wavefunction | none | the Wavefunction that will compute this  |
   *
   * This constructor overrides the default target precision to 1e-9 .
   */
  // clang-format on

  explicit Energy(const KeyVal& kv) : WavefunctionProperty(kv, 1e-9) {}

private:
  void do_evaluate() override;
};

#if 0
/// StationaryPoint finds stationary points on molecular PES.
class StationaryPoint : public Property {
 public:
  explicit StationaryPoint(const KeyVal& kv);
 private:
  std::shared_ptr<Energy> energy_;
  std::shared_ptr<math::QuasiNewtonOptimizer<double,MolecularCoordinates>> optimizer_;

  void evaluate() override;
};
#endif

} // namespace mpqc

#endif //  MPQC4_SRC_MPQC_CHEMISTRY_QC_PROPERTIES_ENERGY_H_
