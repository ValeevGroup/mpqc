//
// Created by Chong Peng on 10/6/16.
//

#include "lcao_wfn.h"

MPQC_CLASS_EXPORT_KEY2(mpqc::qc::LCAOWfn, "LCAOWfn");

namespace mpqc {
namespace qc {

LCAOWfn::LCAOWfn(const KeyVal &kv) : Wfn(kv) {
  lcao_factory_ =
      std::make_shared<LCAOWfn::LCAOIntegral>(*(this->wfn_world()), kv);
}

LCAOWfn::LCAOIntegral &LCAOWfn::lcao_factory() { return *lcao_factory_; }

void LCAOWfn::compute(PropertyBase *pb) {}

double LCAOWfn::value() { return 1.0; }
}
}
