//
// Created by Chong Peng on 10/6/16.
//

#include "lcao_wfn.h"

MPQC_CLASS_EXPORT_KEY2(mpqc::qc::LCAOWfn, "LCAOWfn");

namespace mpqc {
namespace qc {

LCAOWfn::LCAOWfn(const KeyVal &kv) : Wfn(kv) {
  lcao_factory_ =
      std::make_shared<LCAOWfn::LCAOFactoryType>(*(this->wfn_world()), kv);
}

LCAOWfn::LCAOFactoryType &LCAOWfn::lcao_factory() { return *lcao_factory_; }

void LCAOWfn::compute(PropertyBase *pb) {}

double LCAOWfn::value() { return 1.0; }
}
}
