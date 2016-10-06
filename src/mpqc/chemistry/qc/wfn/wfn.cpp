
#include <mpqc/chemistry/qc/wfn/wfn.h>

namespace mpqc{
namespace qc{


Wfn::Wfn(const KeyVal &kv) {

  // first check if wfn_world is provided
  if(kv.exists("wfn_world")){
    wfn_world_ = kv.keyval("wfn_world").class_ptr<qc::WfnWorld>();
  }
  // check if wfn_world exist one level above
  else if(kv.exists("..:wfn_world")){
    wfn_world_ = kv.keyval("..:wfn_world").class_ptr<qc::WfnWorld>();
  }
  // use global provide wfn_world
  else if(kv.exists("$:wfn_world")){
    wfn_world_ = kv.keyval("$:wfn_world").class_ptr<qc::WfnWorld>();
  }
  else{
    throw std::runtime_error("Wfn could not find wfn_world keyval! \n");
  }

}

Wfn::~Wfn() = default;

//////void Wfn::compute(PropertyBase* pb){
////  throw std::runtime_error("Wfn::compute is abstract mehtod! \n");
//}

}
}