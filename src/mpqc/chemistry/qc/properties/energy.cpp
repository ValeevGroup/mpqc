#include <mpqc/chemistry/qc/properties/energy.h>
#include <mpqc/chemistry/qc/wfn/ao_wfn.h>

MPQC_CLASS_EXPORT_KEY2(mpqc::qc::Energy, "Energy");

using namespace mpqc;

qc::Energy::Energy(KeyVal const &kv){}

void qc::Energy::apply(AOWfn *aowfn){
    result_ = 1.0;
}
