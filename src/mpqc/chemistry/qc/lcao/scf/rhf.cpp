//
// Created by Chong Peng on 12/6/16.
//

#include "mpqc/chemistry/qc/lcao/scf/rhf.h"
#include "mpqc/util/keyval/forcelink.h"

#if TA_DEFAULT_POLICY == 0
template class mpqc::lcao::RHF<TA::TensorD, TA::DensePolicy>;
MPQC_CLASS_EXPORT2("RHF", mpqc::lcao::RHF<TA::TensorD, TA::DensePolicy>);

template class mpqc::lcao::RIRHF<TA::TensorD, TA::DensePolicy>;
MPQC_CLASS_EXPORT2("RI-RHF", mpqc::lcao::RIRHF<TA::TensorD, TA::DensePolicy>);

template class mpqc::lcao::DirectRHF<TA::TensorD, TA::DensePolicy>;
MPQC_CLASS_EXPORT2("Direct-RHF",
                   mpqc::lcao::DirectRHF<TA::TensorD, TA::DensePolicy>);

template class mpqc::lcao::DirectRIRHF<TA::TensorD, TA::DensePolicy>;
MPQC_CLASS_EXPORT2("Direct-RI-RHF",
                   mpqc::lcao::DirectRIRHF<TA::TensorD, TA::DensePolicy>);

template class mpqc::lcao::CadfRHF<TA::TensorD, TA::DensePolicy>;
MPQC_CLASS_EXPORT2("Cadf-RHF",
                   mpqc::lcao::CadfRHF<TA::TensorD, TA::DensePolicy>);

template class mpqc::lcao::nrCadfRHF<TA::TensorD, TA::DensePolicy>;
MPQC_CLASS_EXPORT2("nrCadf-RHF",
                   mpqc::lcao::nrCadfRHF<TA::TensorD, TA::DensePolicy>);

#elif TA_DEFAULT_POLICY == 1

template class mpqc::lcao::RHF<TA::TensorD, TA::SparsePolicy>;
MPQC_CLASS_EXPORT2("RHF", mpqc::lcao::RHF<TA::TensorD, TA::SparsePolicy>);

template class mpqc::lcao::RIRHF<TA::TensorD, TA::SparsePolicy>;
MPQC_CLASS_EXPORT2("RI-RHF", mpqc::lcao::RIRHF<TA::TensorD, TA::SparsePolicy>);

template class mpqc::lcao::DirectRHF<TA::TensorD, TA::SparsePolicy>;
MPQC_CLASS_EXPORT2("Direct-RHF",
                   mpqc::lcao::DirectRHF<TA::TensorD, TA::SparsePolicy>);

template class mpqc::lcao::DirectRIRHF<TA::TensorD, TA::SparsePolicy>;
MPQC_CLASS_EXPORT2("Direct-RI-RHF",
                   mpqc::lcao::DirectRIRHF<TA::TensorD, TA::SparsePolicy>);

template class mpqc::lcao::CadfRHF<TA::TensorD, TA::SparsePolicy>;
MPQC_CLASS_EXPORT2("Cadf-RHF",
                   mpqc::lcao::CadfRHF<TA::TensorD, TA::SparsePolicy>);

template class mpqc::lcao::nrCadfRHF<TA::TensorD, TA::SparsePolicy>;
MPQC_CLASS_EXPORT2("nrCadf-RHF",
                   mpqc::lcao::nrCadfRHF<TA::TensorD, TA::SparsePolicy>);

template class mpqc::lcao::RIJEXACTKRHF<TA::TensorD, TA::SparsePolicy>;
MPQC_CLASS_EXPORT2("RIJ-ExactK-RHF",
                   mpqc::lcao::RIJEXACTKRHF<TA::TensorD, TA::SparsePolicy>);
#endif
