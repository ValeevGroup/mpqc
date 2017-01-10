
#include "mpqc/chemistry/qc/wfn/wfn.h"

namespace mpqc{

Wavefunction::Wavefunction(const KeyVal &kv) {
  if (kv.exists_class("atoms")) {
    atoms_ = kv.class_ptr<Molecule>("atoms");
  }
}

Wavefunction::~Wavefunction() = default;

namespace lcao{

Wavefunction::Wavefunction(const KeyVal &kv) : ::mpqc::Wavefunction(kv) {

  // first check if wfn_world is provided
  if(kv.exists("wfn_world")){
    wfn_world_ = kv.keyval("wfn_world").class_ptr<WavefunctionWorld>();
  }
  // check if wfn_world exist one level above
  else if(kv.exists("..:wfn_world")){
    wfn_world_ = kv.keyval("..:wfn_world").class_ptr<WavefunctionWorld>();
  }
  // use global (i.e. top-level) wfn_world
  else if(kv.exists("$:wfn_world")){
    wfn_world_ = kv.keyval("$:wfn_world").class_ptr<WavefunctionWorld>();
  }
  else{
    throw std::runtime_error("Wavefunction could not find wfn_world keyval! \n");
  }

}

Wavefunction::~Wavefunction() = default;

}  // namespace lcao
}  // namespace mpqc
