//
// Created by Eduard Valeyev on 8/5/18.
//

#include "mpqc/util/keyval/forcelink.h"
#include "mpqc/chemistry/qc/lcao/factory/wfn_world.h"

namespace mpqc{
namespace lcao{

WavefunctionWorld::WavefunctionWorld(KeyVal const &kv)
    : ::mpqc::WavefunctionWorld(kv) {
  basis_registry_ = std::make_shared<gaussian::OrbitalBasisRegistry>(kv);
}

}  // namespace lcao
}  // namespace mpqc

MPQC_CLASS_EXPORT2("LCAOWfnWorld", mpqc::lcao::WavefunctionWorld);
