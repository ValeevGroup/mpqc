/*
 * ao_wfn.h
 *
 *  Created on: Aug 17, 2016
 *      Author: Drew Lewis
 */

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_WFN_AO_WFN_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_WFN_AO_WFN_H_

#include "mpqc/chemistry/qc/lcao/wfn/wfn.h"
#include "mpqc/chemistry/qc/lcao/factory/ao_factory.h"
#include "mpqc/chemistry/qc/lcao/factory/periodic_ao_factory.h"

namespace mpqc {
namespace lcao {

/**
 * \brief AOWavefunction is a Wavefunction with an gaussian::AOFactory
 *
 * This models wave function methods expressed in Gaussian AO basis.
 * \todo factor out the dependence on Gaussian basis into a WavefunctionPolicy
 * \todo elaborate AOWavefunction documentation
**/
template<typename Tile, typename Policy>
class AOWavefunction : public Wavefunction {
 public:
  using AOIntegral = gaussian::AOFactory<Tile, Policy>;
  using ArrayType = typename AOIntegral::TArray;

  /**
   *  \brief The KeyVal constructor
   *
   * The KeyVal object will be queried for all keywords of the Wavefunction
   * class,
   * as well as the following keywords:
   * | KeyWord | Type | Default| Description |
   * |---------|------|--------|-------------|
   * | \c "wfn_world:ao_factory" | integrals::AOFactory | default-constructed
   * integrals::AOFactory | |
   * | \c "wfn_world:direct_ao_factory" | integrals::DirectAOFactory |
   * default-constructed integrals::DirectAOFactory | |
   */
  AOWavefunction(const KeyVal &kv) : Wavefunction(kv) {
    ao_factory_ = gaussian::construct_ao_factory<Tile, Policy>(kv);
    ao_factory_->set_orbital_basis_registry(
        this->wfn_world()->basis_registry());

  }
  virtual ~AOWavefunction() = default;

  /// obsolete, purge the registry in AOIntegral and DirectAOIntegral
  void obsolete() override {
    ao_factory().obsolete();
    Wavefunction::obsolete();
  }

  /*! Return a reference to the AOFactory Library
   */
  AOIntegral &ao_factory() { return *ao_factory_; }

  /*! Return a const reference to the AOFactory Library
   */
  const AOIntegral &ao_factory() const { return *ao_factory_; }

private:
  std::shared_ptr<AOIntegral> ao_factory_;
};

/// PeriodicAOWavefunction is a Wavefunction with a gaussian::PeriodicAOFactory

/// \todo factor out the dependence on Gaussian basis into a WavefunctionPolicy
/// This models wave function methods expressed in Gaussian AO basis.
/// \todo elaborate PeriodicAOWavefunction documentation
template <typename Tile, typename Policy>
class PeriodicAOWavefunction : public Wavefunction {
 public:
  using AOIntegral = gaussian::PeriodicAOFactory<Tile, Policy>;
  using ArrayType = typename AOIntegral::TArray;
  using MatrixzVec = std::vector<Matrixz>;
  using VectorzVec = std::vector<Vectorz>;

  /**
   *  \brief The KeyVal constructor
   *
   * The KeyVal object will be queried for all keywords of the Wavefunction
   * class,
   * as well as the following keywords:
   * | KeyWord | Type | Default| Description |
   * |---------|------|--------|-------------|
   * | \c "wfn_world:ao_factory" | integrals::PeriodicAOFactory |
   * default-constructed integrals::PeriodicAOFactory | |
   */
  PeriodicAOWavefunction(const KeyVal &kv) : Wavefunction(kv) {
    ao_factory_ = gaussian::construct_periodic_ao_factory<Tile, Policy>(kv);
    ao_factory_->set_orbital_basis_registry(
        this->wfn_world()->basis_registry());
  }
  virtual ~PeriodicAOWavefunction() = default;

  virtual MatrixzVec co_coeff() = 0;
  virtual VectorzVec co_energy() = 0;
  virtual Vector3i nk() = 0;
  virtual int64_t k_size() = 0;

  /*! Return a reference to the AOFactory Library
   */
  AOIntegral &ao_factory() { return *ao_factory_; }

  /*! Return a const reference to the AOFactory Library
   */
  const AOIntegral &ao_factory() const { return *ao_factory_; }

  std::shared_ptr<AOIntegral> ao_factory_;
};

#if TA_DEFAULT_POLICY == 0
extern template class AOWavefunction<TA::TensorD, TA::DensePolicy>;
#elif TA_DEFAULT_POLICY == 1
extern template class AOWavefunction<TA::TensorD, TA::SparsePolicy>;
extern template class PeriodicAOWavefunction<TA::TensorD, TA::SparsePolicy>;
#endif

}  // namespace lcao
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_WFN_AO_WFN_H_
