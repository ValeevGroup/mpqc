//
// Created by Chong Peng on 10/6/16.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_WFN_LCAO_WFN_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_WFN_LCAO_WFN_H_

#include "mpqc/chemistry/qc/lcao/wfn/wfn.h"
#include "mpqc/util/keyval/keyval.h"
#include "mpqc/chemistry/qc/lcao/integrals/lcao_factory.h"
#include "mpqc/chemistry/qc/lcao/wfn/trange1_engine.h"
#include "mpqc/chemistry/qc/lcao/integrals/periodic_lcao_factory.h"

namespace mpqc{
namespace lcao{

/// LCAOWavefunction is a Wavefunction with an LCAOFactory

/// This models wave function methods expressed in LCAO basis (e.g. traditional
/// electron correlation methods, like MO-basis CCSD).
/// \todo elaborate LCAOWavefunction documentation
template<typename Tile, typename Policy>
class LCAOWavefunction : public Wavefunction {

public:
  using ArrayType = TA::DistArray<Tile, Policy>;
  using LCAOFactoryType = lcao::LCAOFactory<Tile,Policy>;

  /**
   *  \brief The KeyVal constructor
   *
   * The KeyVal object will be queried for all keywords of Wavefunction and LCAOFactory, and the following keywords:
   *
   * | KeyWord | Type | Default| Description |
   * |---------|------|--------|-------------|
   * | \c "frozen_core" | bool | true | if true, core electrons are not correlated |
   * | \c "charge" | int | 0 | the net charge of the molecule |
   * | \c "obs_block_size" | int | 24 | the target OBS (Orbital Basis Set) space block size |
   * | \c "occ_block_size" | int | \c "$obs_block_size" | the target block size of the occupied space |
   * | \c "unocc_block_size" | int | \c "$obs_block_size" | the target block size of the unoccupied space |
   *
   */
  LCAOWavefunction(const KeyVal &kv) : Wavefunction(kv) {
    lcao_factory_ = lcao::detail::construct_lcao_factory<Tile,Policy>(kv);

    frozen_core_ = kv.value<bool>("frozen_core",true);

    const auto net_charge = kv.value<int>("charge", 0);
    if (this->atoms()->total_atomic_number() <= net_charge)
      throw InputError("net charge cannot be greater than the sum of atomic numbers",
                       __FILE__, __LINE__, "charge");
    const auto nelectrons = this->atoms()->total_atomic_number() - net_charge;
    if (nelectrons % 2 != 0)
      throw InputError("LCAOWavefunction for now requires an even number of electrons",
                           __FILE__, __LINE__, "charge");
    ndocc_ = nelectrons / 2;
    std::size_t mo_block = kv.value<int>("obs_block_size",24);
    occ_block_ = kv.value<int>("occ_block_size",mo_block);
    unocc_block_ = kv.value<int>("unocc_block_size",mo_block);
  }

  virtual ~LCAOWavefunction() = default;

  LCAOFactoryType& lcao_factory() {
    return *lcao_factory_;
  }
  void obsolete() override {
    // obsolete factory
    lcao_factory_->obsolete();
    orbital_energy_.reset();
    trange1_engine_.reset();
    // obsolete wfn
    Wavefunction::obsolete();
  }

  const std::shared_ptr<TRange1Engine> trange1_engine() const {
    return trange1_engine_;
  }

  const std::shared_ptr<Eigen::VectorXd> orbital_energy() const {
    return orbital_energy_;
  }

  bool is_frozen_core() const {
    return frozen_core_;
  }
  /// @return # of the doubly-occupied orbitals
  size_t ndocc() const {
    return ndocc_;
  }
  size_t occ_block() const {
    return occ_block_;
  }
  size_t unocc_block() const {
    return unocc_block_;
  }

protected:
  std::shared_ptr<Eigen::VectorXd> orbital_energy_;
  std::shared_ptr<mpqc::TRange1Engine> trange1_engine_;

private:

  std::shared_ptr<LCAOFactoryType> lcao_factory_;
  bool frozen_core_;
  std::size_t ndocc_;
  std::size_t occ_block_;
  std::size_t unocc_block_;

};

/// PeriodicLCAOWavefunction is a Wavefunction with a PeriodicLCAOFactory

/// This models wave function methods expressed in LCAO basis (e.g. traditional
/// electron correlation methods, like MO-basis CCSD).
/// \todo elaborate PeriodicLCAOWavefunction documentation
template<typename Tile, typename Policy>
class PeriodicLCAOWavefunction : public Wavefunction {

public:
  using ArrayType = TA::DistArray<Tile, Policy>;
  using LCAOFactoryType = lcao::PeriodicLCAOFactory<Tile,Policy>;

  /**
   *  \brief The KeyVal constructor
   *
   * The KeyVal object will be queried for all keywords of Wavefunction and LCAOFactory, and the following keywords:
   *
   * | KeyWord | Type | Default| Description |
   * |---------|------|--------|-------------|
   * | \c "frozen_core" | bool | true | if true, core electrons are not correlated |
   * | \c "obs_block_size" | int | 24 | the target OBS (Orbital Basis Set) space block size |
   * | \c "occ_block_size" | int | \c "$obs_block_size" | the target block size of the occupied space |
   * | \c "unocc_block_size" | int | \c "$obs_block_size" | the target block size of the unoccupied space |
   *
   */
  PeriodicLCAOWavefunction(const KeyVal &kv) : Wavefunction(kv) {
    lcao_factory_ = lcao::detail::construct_periodic_lcao_factory<Tile,Policy>(kv);

    frozen_core_ = kv.value<bool>("frozen_core",true);
    std::size_t mo_block = kv.value<int>("obs_block_size",24);
    occ_block_ = kv.value<int>("occ_block_size",mo_block);
    unocc_block_ = kv.value<int>("unocc_block_size",mo_block);
  }

  virtual ~PeriodicLCAOWavefunction() = default;

  LCAOFactoryType& lcao_factory() {
    return *lcao_factory_;
  }
  void obsolete() override {
    lcao_factory_->registry().purge(wfn_world()->world());
    lcao_factory_->orbital_space().clear();
    lcao_factory_->ao_factory().registry().purge(wfn_world()->world());
    Wavefunction::obsolete();
  }

  const std::shared_ptr<TRange1Engine> trange1_engine() const {
    return trange1_engine_;
  }

  const std::shared_ptr<Eigen::VectorXd> orbital_energy() const {
    return orbital_energy_;
  }

  bool is_frozen_core() const {
    return frozen_core_;
  }
  size_t occ_block() const {
    return occ_block_;
  }
  size_t unocc_block() const {
    return unocc_block_;
  }

protected:
  std::shared_ptr<Eigen::VectorXd> orbital_energy_;
  std::shared_ptr<mpqc::TRange1Engine> trange1_engine_;

private:

  std::shared_ptr<LCAOFactoryType> lcao_factory_;
  bool frozen_core_;
  std::size_t occ_block_;
  std::size_t unocc_block_;

};

#if TA_DEFAULT_POLICY == 0
extern template class LCAOWavefunction<TA::TensorD, TA::DensePolicy>;
#elif TA_DEFAULT_POLICY == 1
extern template class LCAOWavefunction<TA::TensorD, TA::SparsePolicy>;
extern template class PeriodicLCAOWavefunction<TA::TensorZ, TA::SparsePolicy>;
#endif

}
}


#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_WFN_LCAO_WFN_H_
