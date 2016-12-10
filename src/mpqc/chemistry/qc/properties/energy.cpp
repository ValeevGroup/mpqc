#include "mpqc/chemistry/qc/properties/energy.h"
#include "mpqc/chemistry/qc/wfn/ao_wfn.h"

MPQC_CLASS_EXPORT_KEY2("Energy", mpqc::lcao::Energy);

using namespace mpqc;

lcao::Energy::Energy(KeyVal const &kv) : result_() {}

void lcao::Energy::apply(Wavefunction *aowfn){
    result_ = 1.0;
}
