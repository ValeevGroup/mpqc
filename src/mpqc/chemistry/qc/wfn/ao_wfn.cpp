#include <mpqc/chemistry/qc/properties/propertybase.h>
#include <mpqc/chemistry/qc/wfn/ao_wfn.h>

MPQC_CLASS_EXPORT_KEY(mpqc::qc::AOWfn);

using namespace mpqc;
qc::AOWfn::AOWfn(KeyVal const &kv) : Wfn(kv) {}

void qc::AOWfn::compute(qc::PropertyBase *pb) { pb->apply(this); }
