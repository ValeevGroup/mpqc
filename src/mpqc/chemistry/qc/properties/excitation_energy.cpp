//
// Created by Chong Peng on 3/2/17.
//

#include "mpqc/chemistry/qc/properties/excitation_energy.h"
#include "mpqc/util/keyval/forcelink.h"

MPQC_CLASS_EXPORT2("ExcitationEnergy", mpqc::ExcitationEnergy);

namespace mpqc {

ExcitationEnergy::ExcitationEnergy(const KeyVal &kv)
    : WavefunctionProperty<std::vector<double>>(kv, 1.0e-6),
      n_roots_(kv.value<int>("n_roots", 1)),
      singlets_(kv.value<bool>("singlets", true)),
      triplets_(kv.value<bool>("triplets", false)) {
  if (n_roots_ < 1) {
    throw InputError("ExcitationEnergy: number of roots must be greater than 1",
                     __FILE__, __LINE__, "n_roots");
  }
}

unsigned int ExcitationEnergy::n_roots() const { return n_roots_; }

bool ExcitationEnergy::singlets() const { return singlets_; }

bool ExcitationEnergy::triplets() const { return triplets_; }

void ExcitationEnergy::do_evaluate() {
  auto evaluator = std::dynamic_pointer_cast<Provider>(wfn());
  if (evaluator == nullptr) {
    std::ostringstream oss;
    oss << "Wavefunction " << wfn()->class_key()
        << " cannot compute ExcitationEnergy" << std::endl;
    throw InputError(oss.str().c_str(), __FILE__, __LINE__);
  }
  evaluator->evaluate(this);
}

}  // namespace mpqc
