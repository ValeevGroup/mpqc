
#include "mpqc/chemistry/qc/wfn/wfn.h"
#include "mpqc/util/core/exception.h"

namespace mpqc{

WavefunctionWorld::WavefunctionWorld(const KeyVal& kv) :
world_(*kv.value<madness::World *>("$:world")),
atoms_(kv.class_ptr<Molecule>("atoms"))
{
  if (!atoms_)
    atoms_ = kv.class_ptr<Molecule>("molecule");
  if (!atoms_)
    throw InputError(
        "WavefunctionWorld did not find keyword \"atoms\"",
        __FILE__, __LINE__);
}

Wavefunction::Wavefunction(const KeyVal &kv) {

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

Wavefunction::~Wavefunction() { }

}  // namespace mpqc
