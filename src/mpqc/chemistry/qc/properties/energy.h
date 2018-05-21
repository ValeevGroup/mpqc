//
// Created by Chong Peng on 1/11/17.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_PROPERTIES_ENERGY_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_PROPERTIES_ENERGY_H_

#include "mpqc/chemistry/qc/properties/property.h"
#include "mpqc/math/function/optimize.h"

namespace mpqc {

/// Taylor expansion of the molecular energy computed by a Wavefunction
class Energy : public WavefunctionProperty<double> {
public:
  using typename WavefunctionProperty<double>::function_base_type;

  /**
   *  every class that can evaluate Energy (e.g. Wavefunction) will publicly
   *  inherit from Energy::Provider
   *
   *  @sa Provides
   */
  class Provider : public math::FunctionVisitorBase<function_base_type> {
  public:
    /// @return true if \c energy can be computed.
    /// For example, if \c energy demands taylor expansion to 1st order
    /// but this wave function does not have analytic nuclear gradients,
    /// will return false.
    virtual bool can_evaluate(Energy* energy) = 0;
    /// Provider::evaluate computes the taylor expansion of the energy
    /// and uses set_value to assign the values to \c energy
    virtual void evaluate(Energy* energy) = 0;
  };

  // clang-format off
  /**
   * @brief The KeyVal constructor
   * @param kv the KeyVal object, it will be queried for all
   *        keywords of the WavefunctionProperty class. |
   *
   * @note This constructor overrides the default target precision to 1e-9 .
   */
  // clang-format on
  explicit Energy(const KeyVal& kv) : WavefunctionProperty(kv, 1e-9) {}

  /// This constructor is provided for convenient programmatic construction
  /// Request to compute energy value using Wavefunction.\c wfn to precision \c prec.
  /// @param wfn the Wavefunction object that will compute the energy value
  /// @param prec the targer precision of the energy
  Energy(std::shared_ptr<Wavefunction> wfn, double prec) :
    WavefunctionProperty(make_kv(wfn), prec) {}

  double energy() { return this->value()->value(); }
  const std::vector<double>& gradient() { return this->value()->gradient(); }
  const std::vector<double>& hessian() { return this->value()->hessian(); }

private:
  void do_evaluate() override;

  static KeyVal make_kv(std::shared_ptr<Wavefunction> wfn) {
    KeyVal kv;
    kv.assign("wfn", wfn);
    return kv;
  }

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
