//
// Created by Chong Peng on 10/6/16.
//

#ifndef MPQC_CHEMISTRY_QC_WFN_LCAO_WFN_H_
#define MPQC_CHEMISTRY_QC_WFN_LCAO_WFN_H_

#include <mpqc/chemistry/qc/wfn/wfn.h>
#include <mpqc/util/keyval/keyval.hpp>
#include <mpqc/chemistry/qc/integrals/lcao_factory.h>
#include "mpqc/chemistry/qc/wfn/trange1_engine.h"

namespace mpqc{
namespace qc{

/// LCAOWavefunction is a Wavefunction with an LCAOFactory

/// This models wave function methods expressed in LCAO basis (e.g. traditional
/// electron correlation methods, like MO-basis CCSD).
/// \todo elaborate LCAOWavefunction documentation
template<typename Tile, typename Policy>
class LCAOWavefunction : public Wavefunction {

public:
  using ArrayType = TA::DistArray<Tile, Policy>;
  using LCAOFactoryType = integrals::LCAOFactory<Tile,Policy>;

  /**
   * The KeyVal constructor uses keywords of Wavefunction and LCAOFactory, and the following keywords:
   *
   * | KeyWord | Type | Default| Description |
   * |---------|------|--------|-------------|
   * | \c "frozen_core" | bool | true | if true, core electrons are not correlated |
   * | \c "obs_block_size" | int | 24 | the target OBS (Orbital Basis Set) space block size |
   * | \c "occ_block_size" | int | \c "$obs_block_size" | the target block size of the occupied space |
   * | \c "uocc_block_size" | int | \c "$obs_block_size" | the target block size of the unoccupied space |
   *
   */
  LCAOWavefunction(const KeyVal &kv) : Wavefunction(kv) {
    lcao_factory_ = integrals::detail::construct_lcao_factory<Tile,Policy>(kv);

    frozen_core_ = kv.value<bool>("frozen_core",true);
    std::size_t mo_block = kv.value<int>("obs_block_size",24);
    occ_block_ = kv.value<int>("occ_block_size",mo_block);
    unocc_block_ = kv.value<int>("uocc_block_size",mo_block);
  }

  virtual ~LCAOWavefunction() = default;

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

}
}


#endif //MPQC_CHEMISTRY_QC_WFN_LCAO_WFN_H_
