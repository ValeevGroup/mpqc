#include "eom_ccsd.h"
#include "ea_eom_ccsd.h"
#include "ip_eom_ccsd.h"
#include "mpqc/util/keyval/forcelink.h"

#if TA_DEFAULT_POLICY == 0
template class mpqc::lcao::EOM_CCSD<TA::TensorD, TA::DensePolicy>;
MPQC_CLASS_EXPORT2("EOM-CCSD",
                   mpqc::lcao::EOM_CCSD<TA::TensorD, TA::DensePolicy>);
template class mpqc::lcao::IP_EOM_CCSD<TA::TensorD, TA::DensePolicy>;
MPQC_CLASS_EXPORT2("IP-EOM-CCSD",
                   mpqc::lcao::IP_EOM_CCSD<TA::TensorD, TA::DensePolicy>);
template class mpqc::lcao::EA_EOM_CCSD<TA::TensorD, TA::DensePolicy>;
MPQC_CLASS_EXPORT2("EA-EOM-CCSD",
                   mpqc::lcao::EA_EOM_CCSD<TA::TensorD, TA::DensePolicy>);
#elif TA_DEFAULT_POLICY == 1
template class mpqc::lcao::EOM_CCSD<TA::TensorD, TA::SparsePolicy>;
MPQC_CLASS_EXPORT2("EOM-CCSD",
                   mpqc::lcao::EOM_CCSD<TA::TensorD, TA::SparsePolicy>);
template class mpqc::lcao::EA_EOM_CCSD<TA::TensorD, TA::SparsePolicy>;
MPQC_CLASS_EXPORT2("EA-EOM-CCSD",
                   mpqc::lcao::EA_EOM_CCSD<TA::TensorD, TA::SparsePolicy>);
template class mpqc::lcao::IP_EOM_CCSD<TA::TensorD, TA::SparsePolicy>;
MPQC_CLASS_EXPORT2("IP-EOM-CCSD",
                   mpqc::lcao::IP_EOM_CCSD<TA::TensorD, TA::SparsePolicy>);
#endif
