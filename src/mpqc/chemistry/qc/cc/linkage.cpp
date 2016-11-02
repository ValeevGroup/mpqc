/*
 * linkage.cpp
 *
 *  Created on: Nov 2, 2016
 *      Author: evaleev
 */

#include "mpqc/chemistry/qc/cc/ccsd_t.h"
#include "mpqc/chemistry/qc/cc/dbccsd.h"
#include "mpqc/util/keyval/forcelink.h"

MPQC_FORCELINK_KEYVAL_CTOR(mpqc::cc::CCSD<TA::TensorD, TA::SparsePolicy>)
MPQC_CLASS_EXPORT_KEY2("CCSD", mpqc::cc::CCSD<TA::TensorD, TA::SparsePolicy>);

MPQC_FORCELINK_KEYVAL_CTOR(mpqc::cc::CCSD_T<TA::TensorD, TA::SparsePolicy>)
MPQC_CLASS_EXPORT_KEY2("CCSD(T)", mpqc::cc::CCSD_T<TA::TensorD, TA::SparsePolicy>);

MPQC_FORCELINK_KEYVAL_CTOR(mpqc::cc::DBCCSD<TA::TensorD, TA::SparsePolicy>)
MPQC_CLASS_EXPORT_KEY2("DBCCSD", mpqc::cc::DBCCSD<TA::TensorD, TA::SparsePolicy>);
