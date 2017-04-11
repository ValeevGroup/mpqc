
#include "mpqc/chemistry/qc/lcao/scf/zrhf.h"
#include "mpqc/util/keyval/forcelink.h"

#if TA_DEFAULT_POLICY == 0

#elif TA_DEFAULT_POLICY == 1

template class mpqc::lcao::zRHF<TA::TensorD, TA::SparsePolicy>;
MPQC_CLASS_EXPORT2("zRHF", mpqc::lcao::zRHF<TA::TensorD, TA::SparsePolicy>);

template class mpqc::lcao::DFzRHF<TA::TensorD, TA::SparsePolicy>;
MPQC_CLASS_EXPORT2("DF-zRHF", mpqc::lcao::DFzRHF<TA::TensorD, TA::SparsePolicy>);

template class mpqc::lcao::FourCenterzRHF<TA::TensorD, TA::SparsePolicy>;
MPQC_CLASS_EXPORT2("FourCenter-zRHF", mpqc::lcao::FourCenterzRHF<TA::TensorD, TA::SparsePolicy>);

#endif


