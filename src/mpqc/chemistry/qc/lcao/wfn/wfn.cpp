
#include "mpqc/chemistry/qc/lcao/wfn/wfn.h"

namespace mpqc{
namespace lcao{

Wavefunction::Wavefunction(const KeyVal &kv) : ::mpqc::Wavefunction(kv) {

  // first check if wfn_world is provided
  if(kv.exists("wfn_world")){
    wfn_world_ = kv.class_ptr<WavefunctionWorld>("wfn_world");
  }
  // check if wfn_world exist one level above
  else if(kv.exists("..:wfn_world")){
    wfn_world_ = kv.class_ptr<WavefunctionWorld>("..:wfn_world");
  }
  // use global (i.e. top-level) wfn_world
  else if(kv.exists("$:wfn_world")){
    wfn_world_ = kv.class_ptr<WavefunctionWorld>("$:wfn_world");
  }
  else{
    wfn_world_ = std::make_shared<WavefunctionWorld>(kv);
    // TODO decide whether to push the newly-constructed world to kv, and where
  }

  utility::Observer::register_message(wfn_world_->atoms().get(), [this](){
    this->obsolete();
  });
}

Wavefunction::~Wavefunction() = default;

}  // namespace lcao
}  // namespace mpqc
