//
// Created by Chong Peng on 10/6/16.
//

#include "lcao_wfn.h"

namespace mpqc {
namespace qc {

LCAOWavefunction::LCAOWavefunction(const KeyVal &kv) : Wavefunction(kv) {
  lcao_factory_ =
      std::make_shared<LCAOWavefunction::LCAOFactoryType>(*(this->wfn_world()), kv);

  frozen_core_ = kv.value<bool>("frozen_core",true);
  std::size_t mo_block = kv.value<int>("mo_block",24);
  occ_block_ = kv.value<int>("occ_block",mo_block);
  unocc_block_ = kv.value<int>("un_occ_block",mo_block);
}

LCAOWavefunction::LCAOFactoryType &LCAOWavefunction::lcao_factory() { return *lcao_factory_; }

void LCAOWavefunction::compute(PropertyBase *pb) {}

double LCAOWavefunction::value() { return 1.0; }

void LCAOWavefunction::obsolete() {
  lcao_factory_->registry().purge(wfn_world()->world());
  lcao_factory_->orbital_space().clear();
  lcao_factory_->atomic_integral().registry().purge(wfn_world()->world());
}

bool LCAOWavefunction::is_frozen_core() const {
  return frozen_core_;
}
size_t LCAOWavefunction::occ_block() const {
  return occ_block_;
}
size_t LCAOWavefunction::unocc_block() const {
  return unocc_block_;
}
}
}
