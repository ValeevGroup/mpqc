
#include <mpqc/chemistry/qc/wfn/wfn.h>

namespace mpqc{
namespace qc{


Wavefunction::Wavefunction(const KeyVal &kv) {

  // first check if wfn_world is provided
  if(kv.exists("wfn_world")){
    wfn_world_ = kv.keyval("wfn_world").class_ptr<qc::WavefunctionWorld>();
  }
  // check if wfn_world exist one level above
  else if(kv.exists("..:wfn_world")){
    wfn_world_ = kv.keyval("..:wfn_world").class_ptr<qc::WavefunctionWorld>();
  }
  // use global (i.e. top-level) wfn_world
  else if(kv.exists("$:wfn_world")){
    wfn_world_ = kv.keyval("$:wfn_world").class_ptr<qc::WavefunctionWorld>();
  }
  else{
    throw std::runtime_error("Wavefunction could not find wfn_world keyval! \n");
  }

}

Wavefunction::~Wavefunction() = default;

}
}
