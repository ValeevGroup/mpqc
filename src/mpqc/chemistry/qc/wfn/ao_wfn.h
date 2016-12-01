/*
 * ao_wfn.h
 *
 *  Created on: Aug 17, 2016
 *      Author: Drew Lewis
 */

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_WFN_AO_WFN_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_WFN_AO_WFN_H_

#include "mpqc/chemistry/qc/wfn/wfn.h"
#include "mpqc/chemistry/qc/integrals/ao_factory.h"
#include "mpqc/chemistry/qc/integrals/direct_ao_factory.h"

namespace mpqc {
namespace qc {

/// AOWavefunction is a Wavefunction with an AOFactory

/// This models wave function methods expressed in AO basis.
/// \todo elaborate AOWavefunction documentation
template<typename Tile, typename Policy>
class AOWavefunction : public Wavefunction {
 public:
  using AOIntegral = integrals::AOFactory<Tile, Policy>;
  using DirectAOIntegral = integrals::DirectAOFactory<Tile, Policy>;
  using ArrayType = typename AOIntegral::TArray;

  AOWavefunction(const KeyVal &kv) : Wavefunction(kv)
  {
    ao_factory_ = integrals::detail::construct_ao_factory<Tile, Policy>(kv);
    ao_factory_->set_orbital_basis_registry(this->wfn_world()->basis_registry());

    direct_ao_factory_ = integrals::detail::construct_direct_ao_factory<Tile,Policy>(kv);
    direct_ao_factory_->set_orbital_basis_registry(this->wfn_world()->basis_registry());

  }
  virtual ~AOWavefunction() = default;

  void compute(PropertyBase *pb) override {
    throw std::logic_error("Not Implemented!");
  }

  /// obsolete, purge the registry in AOIntegral and DirectAOIntegral
  void obsolete() override {
    ao_factory().registry().purge(wfn_world()->world());
    direct_ao_factory().registry().purge(wfn_world()->world());
    Wavefunction::obsolete();
  }

  double value() override {
    return 0.0;
  };

  /*! Return a reference to the AOFactory Library
   *
   * \note This reference can't be made const without modifying the
   * AOFactory library so that certain members are mutable.
   */
  AOIntegral &ao_factory() { return *ao_factory_; }

  /// return a reference to DirectAOFactory
  DirectAOIntegral& direct_ao_factory() { return *direct_ao_factory_;}
private:
  std::shared_ptr<AOIntegral> ao_factory_;
  std::shared_ptr<DirectAOIntegral> direct_ao_factory_;
};

}  // namespace qc
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_WFN_AO_WFN_H_
