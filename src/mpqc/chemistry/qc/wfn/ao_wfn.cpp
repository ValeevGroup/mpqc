#include <mpqc/chemistry/qc/properties/propertybase.h>
#include <mpqc/chemistry/qc/wfn/ao_wfn.h>

//MPQC_CLASS_EXPORT_KEY2(mpqc::qc::AOWfn, "AOWfn");

namespace mpqc{
namespace qc{

AOWfn::AOWfn(KeyVal const &kv) : Wfn(kv) {}

AOWfn::~AOWfn() = default;

void AOWfn::compute(qc::PropertyBase *pb) { pb->apply(this); }

}//namespace qc
}//namespace mpqc

