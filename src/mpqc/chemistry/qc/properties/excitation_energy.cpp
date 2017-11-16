//
// Created by Chong Peng on 3/2/17.
//

#include "mpqc/chemistry/qc/properties/excitation_energy.h"
#include "mpqc/util/keyval/forcelink.h"

MPQC_CLASS_EXPORT2("ExcitationEnergy", mpqc::ExcitationEnergy);

namespace mpqc {

ExcitationEnergy::ExcitationEnergy(const KeyVal &kv)
    : WavefunctionProperty<std::vector<double>>(kv, 1.0e-6),
      n_roots_(kv.value<int>("n_roots", 3)),
      n_guess_(kv.value<int>("n_guess", 2*n_roots_)),
      singlets_(kv.value<bool>("singlets", true)),
      triplets_(kv.value<bool>("triplets", false)) {
  if (n_roots_ < 1) {
    throw InputError("ExcitationEnergy: number of roots must be greater than 1",
                     __FILE__, __LINE__, "n_roots");
  }

  if (n_guess_ < n_roots_) {
    throw InputError(
        "ExcitationEnergy: number of guess must be greater than number of "
        "roots",
        __FILE__, __LINE__, "n_guess");
  }
}

std::size_t ExcitationEnergy::n_roots() const { return n_roots_; }

std::size_t ExcitationEnergy::n_guess() const { return n_guess_; }

void ExcitationEnergy::set_n_roots(unsigned int n_roots) {
  if (n_guess_ < n_roots) {
    throw ProgrammingError(
        "ExcitationEnergy: number of guess must be greater than number of "
            "roots",
        __FILE__, __LINE__, "set_n_roots");
  }
  n_roots_ = n_roots;
}

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
