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


std::vector<double> Energy::get_precision(const KeyVal& kv) {
  if (kv.exists("precision") && kv.count("precision") > 0) { // given an array of precisions
    return kv.value<std::vector<double>>("precision");
  }
  else {
    auto deriv_order = kv.value<size_t>("deriv_order",0);
    std::cout << "Energy::deriv_order = " << deriv_order << std::endl;
    auto precision = kv.value<double>("precision",1.0e-8);
    return std::vector<double>(deriv_order+1, precision);
  }
}

}
