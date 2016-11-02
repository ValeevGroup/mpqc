/*
 * linkage.cpp
 *
 *  Created on: Nov 2, 2016
 *      Author: evaleev
 */

#include "mpqc/chemistry/qc/scf/rhf.h"
#include "mpqc/util/keyval/forcelink.h"

MPQC_FORCELINK_KEYVAL_CTOR(mpqc::scf::RHF)
MPQC_FORCELINK_KEYVAL_CTOR(mpqc::scf::RIRHF)
MPQC_FORCELINK_KEYVAL_CTOR(mpqc::scf::DirectRHF)
MPQC_FORCELINK_KEYVAL_CTOR(mpqc::scf::DirectRIRHF)
MPQC_CLASS_EXPORT_KEY2("RHF", mpqc::scf::RHF);
MPQC_CLASS_EXPORT_KEY2("Direct-RHF", mpqc::scf::DirectRHF);
MPQC_CLASS_EXPORT_KEY2("RI-RHF", mpqc::scf::RIRHF);
MPQC_CLASS_EXPORT_KEY2("Direct-RI-RHF", mpqc::scf::DirectRIRHF);
