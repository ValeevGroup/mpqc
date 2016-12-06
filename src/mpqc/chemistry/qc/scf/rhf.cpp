//
// Created by Chong Peng on 12/6/16.
//

#include "mpqc/chemistry/qc/scf/rhf.h"
#include "mpqc/util/keyval/forcelink.h"


template class mpqc::scf::RHF<TA::TensorD, TA::SparsePolicy>;
template class mpqc::scf::RIRHF<TA::TensorD, TA::SparsePolicy>;
template class mpqc::scf::DirectRHF<TA::TensorD, TA::SparsePolicy>;
template class mpqc::scf::DirectRIRHF<TA::TensorD, TA::SparsePolicy>;

MPQC_CLASS_EXPORT2("RHF", mpqc::scf::RHF<TA::TensorD, TA::SparsePolicy>);
MPQC_CLASS_EXPORT2("Direct-RHF", mpqc::scf::DirectRHF<TA::TensorD,TA::SparsePolicy>);
MPQC_CLASS_EXPORT2("RI-RHF", mpqc::scf::RIRHF<TA::TensorD,TA::SparsePolicy>);
MPQC_CLASS_EXPORT2("Direct-RI-RHF", mpqc::scf::DirectRIRHF<TA::TensorD,TA::SparsePolicy>);