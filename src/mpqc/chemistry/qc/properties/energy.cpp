#include <mpqc/chemistry/qc/properties/energy.h>
#include <mpqc/chemistry/qc/wfn/ao_wfn.h>

MPQC_CLASS_EXPORT_KEY2("Energy", mpqc::qc::Energy);

using namespace mpqc;

qc::Energy::Energy(KeyVal const &kv) : result_() {}

void qc::Energy::apply(AOWavefunction *aowfn){
    result_ = 1.0;
}
