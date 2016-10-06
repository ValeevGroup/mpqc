//
// Created by Chong Peng on 10/6/16.
//

#include "lcao_wfn.h"

namespace mpqc {
namespace qc {

LCAOWavefunction::LCAOWavefunction(const KeyVal &kv) : Wavefunction(kv) {
  lcao_factory_ =
      std::make_shared<LCAOWavefunction::LCAOFactoryType>(*(this->wfn_world()), kv);
}

LCAOWavefunction::LCAOFactoryType &LCAOWavefunction::lcao_factory() { return *lcao_factory_; }

void LCAOWavefunction::compute(PropertyBase *pb) {}

double LCAOWavefunction::value() { return 1.0; }
}
}
