/*
 * linkage.cpp
 *
 *  Created on: Nov 2, 2016
 *      Author: evaleev
 */

#include "mpqc/chemistry/qc/scf/rhf.h"
#include "mpqc/util/keyval/forcelink.h"
#include "zrhf.h"

MPQC_CLASS_EXPORT2("RHF", mpqc::lcao::RHF);
MPQC_CLASS_EXPORT2("Direct-RHF", mpqc::lcao::DirectRHF);
MPQC_CLASS_EXPORT2("RI-RHF", mpqc::lcao::RIRHF);
MPQC_CLASS_EXPORT2("Direct-RI-RHF", mpqc::lcao::DirectRIRHF);
MPQC_CLASS_EXPORT2("zRHF", mpqc::lcao::zRHF);
