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
#include "mpqc/chemistry/qc/integrals/periodic_ao_factory.h"

namespace mpqc {
namespace lcao {

/// AOWavefunction is a Wavefunction with an gaussian::AOFactory

/// This models wave function methods expressed in Gaussian AO basis.
/// \todo factor out the dependence on Gaussian basis into a WavefunctionPolicy
/// \todo elaborate AOWavefunction documentation
template<typename Tile, typename Policy>
class AOWavefunction : public Wavefunction {
 public:
  using AOIntegral = gaussian::AOFactory<Tile, Policy>;
  using DirectAOIntegral = gaussian::DirectAOFactory<Tile, Policy>;
  using ArrayType = typename AOIntegral::TArray;

  /**
   *  \brief The KeyVal constructor
   *
   * The KeyVal object will be queried for all keywords of the Wavefunction class,
   * as well as the following keywords:
   * | KeyWord | Type | Default| Description |
   * |---------|------|--------|-------------|
   * | \c "wfn_world:ao_factory" | integrals::AOFactory | default-constructed integrals::AOFactory | |
   * | \c "wfn_world:direct_ao_factory" | integrals::DirectAOFactory | default-constructed integrals::DirectAOFactory | |
   */
  AOWavefunction(const KeyVal &kv) : Wavefunction(kv)
  {
    ao_factory_ = gaussian::construct_ao_factory<Tile, Policy>(kv);
    ao_factory_->set_orbital_basis_registry(this->wfn_world()->basis_registry());

    direct_ao_factory_ = gaussian::construct_direct_ao_factory<Tile,Policy>(kv);
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

/// PeriodicAOWavefunction is a Wavefunction with a gaussian::PeriodicAOFactory

/// \todo factor out the dependence on Gaussian basis into a WavefunctionPolicy
/// This models wave function methods expressed in Gaussian AO basis.
/// \todo elaborate PeriodicAOWavefunction documentation
template<typename Tile, typename Policy>
class PeriodicAOWavefunction : public Wavefunction {
 public:
  using AOIntegral = gaussian::PeriodicAOFactory<Tile, Policy>;
  using ArrayType = typename AOIntegral::TArray;

  /**
   *  \brief The KeyVal constructor
   *
   * The KeyVal object will be queried for all keywords of the Wavefunction class,
   * as well as the following keywords:
   * | KeyWord | Type | Default| Description |
   * |---------|------|--------|-------------|
   * | \c "wfn_world:ao_factory" | integrals::PeriodicAOFactory | default-constructed integrals::PeriodicAOFactory | |
   */
  PeriodicAOWavefunction(const KeyVal &kv) : Wavefunction(kv)
  {
    ao_factory_ = gaussian::construct_periodic_ao_factory<Tile, Policy>(kv);
    ao_factory_->set_orbital_basis_registry(this->wfn_world()->basis_registry());
  }
  virtual ~PeriodicAOWavefunction() = default;

  void compute(PropertyBase *pb) override {
    throw std::logic_error("Not Implemented!");
  }

  /// obsolete, purge the registry in AOIntegral
  void obsolete() override {
    //ao_factory().registry().purge(wfn_world()->world());
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

private:
  std::shared_ptr<AOIntegral> ao_factory_;
};

}  // namespace lcao
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_WFN_AO_WFN_H_
