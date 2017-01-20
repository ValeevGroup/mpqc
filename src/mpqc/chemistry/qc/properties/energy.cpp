//
// Created by Chong Peng on 1/11/17.
//

#include "mpqc/chemistry/qc/properties/energy.h"
#include "mpqc/util/keyval/forcelink.h"

MPQC_CLASS_EXPORT2("Energy", mpqc::Energy);

namespace mpqc{

void Energy::do_evaluate() {
  auto evaluator = dynamic_cast<Evaluator*>(wfn());
  if (evaluator == nullptr) {
    std::ostringstream oss;
    // TODO Must implement DescribedClass::key() instead of using RTTI's
    // unpredicable output
    oss << "Wavefunction " << typeid(*wfn()).name()
        << " cannot compute Energy" << std::endl;
    throw InputError(oss.str().c_str(), __FILE__, __LINE__);
  }
  evaluator->evaluate(this);
}

}
