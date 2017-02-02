#include "mpqc/chemistry/qc/properties/gfpole.h"
#include "mpqc/util/keyval/forcelink.h"

MPQC_CLASS_EXPORT2("GFRealPole", mpqc::GFRealPole);

namespace mpqc {

GFRealPole::GFRealPole(const KeyVal& kv)
    : WavefunctionProperty(kv, 1e-4), target_(kv.value<int>("target", -1)) {
  if (target_ == 0)
    throw std::runtime_error(
        "GFRealPole: target must be positive (for particles) or negative (for "
        "holes)");
}

int GFRealPole::target() const { return target_; }

void GFRealPole::do_evaluate() {
  auto evaluator = std::dynamic_pointer_cast<Provider>(wfn());
  if (evaluator == nullptr) {
    std::ostringstream oss;
    oss << "Wavefunction " << wfn()->class_key() << " cannot compute GFRealPole"
        << std::endl;
    throw InputError(oss.str().c_str(), __FILE__, __LINE__);
  }
  evaluator->evaluate(this);
}

}  // namespace mpqc
