
#include "mpqc/chemistry/qc/lcao/scf/zrhf.h"
#include "mpqc/util/keyval/forcelink.h"

#if TA_DEFAULT_POLICY == 0

#elif TA_DEFAULT_POLICY == 1

template class mpqc::lcao::zRHF<TA::TensorD, TA::SparsePolicy>;
MPQC_CLASS_EXPORT2("zRHF", mpqc::lcao::zRHF<TA::TensorD, TA::SparsePolicy>);

template class mpqc::lcao::DFzRHF<TA::TensorD, TA::SparsePolicy>;
MPQC_CLASS_EXPORT2("DF-zRHF",
                   mpqc::lcao::DFzRHF<TA::TensorD, TA::SparsePolicy>);

template class mpqc::lcao::FourCenterzRHF<TA::TensorD, TA::SparsePolicy>;
MPQC_CLASS_EXPORT2("FourCenter-zRHF",
                   mpqc::lcao::FourCenterzRHF<TA::TensorD, TA::SparsePolicy>);

template class mpqc::lcao::RIJCADFKzRHF<TA::TensorD, TA::SparsePolicy>;
MPQC_CLASS_EXPORT2("RIJ-CADFK-zRHF",
                   mpqc::lcao::RIJCADFKzRHF<TA::TensorD, TA::SparsePolicy>);

template class mpqc::lcao::FourCenterJCADFKzRHF<TA::TensorD, TA::SparsePolicy>;
MPQC_CLASS_EXPORT2(
    "FourCenterJ-CADFK-zRHF",
    mpqc::lcao::FourCenterJCADFKzRHF<TA::TensorD, TA::SparsePolicy>);

template class mpqc::lcao::MARIJCADFKzRHF<TA::TensorD, TA::SparsePolicy>;
MPQC_CLASS_EXPORT2("MA-RIJ-CADFK-zRHF",
                   mpqc::lcao::MARIJCADFKzRHF<TA::TensorD, TA::SparsePolicy>);

template class mpqc::lcao::MARIJFourCenterKzRHF<TA::TensorD, TA::SparsePolicy>;
MPQC_CLASS_EXPORT2(
    "MA-RIJ-FourCenterK-zRHF",
    mpqc::lcao::MARIJFourCenterKzRHF<TA::TensorD, TA::SparsePolicy>);

template class mpqc::lcao::MAFourCenterzRHF<TA::TensorD, TA::SparsePolicy>;
MPQC_CLASS_EXPORT2("MA-FourCenter-zRHF",
                   mpqc::lcao::MAFourCenterzRHF<TA::TensorD, TA::SparsePolicy>);

template class mpqc::lcao::MAFourCenterJCADFKzRHF<TA::TensorD,
                                                  TA::SparsePolicy>;
MPQC_CLASS_EXPORT2(
    "MA-FourCenterJ-CADFK-zRHF",
    mpqc::lcao::MAFourCenterJCADFKzRHF<TA::TensorD, TA::SparsePolicy>);
#endif
