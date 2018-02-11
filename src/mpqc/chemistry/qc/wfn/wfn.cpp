
#include "mpqc/chemistry/qc/wfn/wfn.h"
#include "mpqc/util/core/exception.h"

namespace mpqc{

Wavefunction::Wavefunction(const KeyVal &kv) {
  atoms_ = kv.class_ptr<Molecule>("atoms");
  if (!atoms_) {
    atoms_ = kv.class_ptr<Molecule>("molecule");
  }
  if (!atoms_)
    throw InputError("Wavefunction did not find atoms", __FILE__, __LINE__, "atoms");

  utility::Observer::register_message(atoms_.get(), [this](){
    this->obsolete();
  });
}

Wavefunction::~Wavefunction() { }

}  // namespace mpqc
