//
// Created by Chong Peng on 1/11/17.
//

#include "mpqc/chemistry/qc/properties/energy.h"
#include "mpqc/util/keyval/forcelink.h"

MPQC_CLASS_EXPORT2("Energy", mpqc::Energy);

namespace mpqc{

void Energy::do_evaluate() {
  auto evaluator = std::dynamic_pointer_cast<Provider>(wfn());
  if (evaluator == nullptr) {
    std::ostringstream oss;
    oss << "Wavefunction " << wfn()->class_key()
        << " cannot compute Energy" << std::endl;
    throw InputError(oss.str().c_str(), __FILE__, __LINE__);
  }
  evaluator->evaluate(this);
}

#if 0
void StationaryPoint::evaluate() {
  optimizer_->value();
}
#endif

}
