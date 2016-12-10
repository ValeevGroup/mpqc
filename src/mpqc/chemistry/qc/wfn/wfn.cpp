
#include "mpqc/chemistry/qc/wfn/wfn.h"

namespace mpqc{
namespace lcao{

Wavefunction::Wavefunction(const KeyVal &kv) : Energy(kv) {

  // first check if wfn_world is provided
  if(kv.exists("wfn_world")){
    wfn_world_ = kv.keyval("wfn_world").class_ptr<lcao::WavefunctionWorld>();
  }
  // check if wfn_world exist one level above
  else if(kv.exists("..:wfn_world")){
    wfn_world_ = kv.keyval("..:wfn_world").class_ptr<lcao::WavefunctionWorld>();
  }
  // use global (i.e. top-level) wfn_world
  else if(kv.exists("$:wfn_world")){
    wfn_world_ = kv.keyval("$:wfn_world").class_ptr<lcao::WavefunctionWorld>();
  }
  else{
    throw std::runtime_error("Wavefunction could not find wfn_world keyval! \n");
  }

}

Wavefunction::~Wavefunction() = default;

}  // namespace lcao
}  // namespace mpqc
