/*
 * linkage.cpp
 *
 *  Created on: Nov 2, 2016
 *      Author: evaleev
 */

#include "mpqc/chemistry/qc/scf/rhf.h"
#include "mpqc/util/keyval/forcelink.h"
#include "zrhf.h"

MPQC_CLASS_EXPORT2("RHF", mpqc::scf::RHF);
MPQC_CLASS_EXPORT2("Direct-RHF", mpqc::scf::DirectRHF);
MPQC_CLASS_EXPORT2("RI-RHF", mpqc::scf::RIRHF);
MPQC_CLASS_EXPORT2("Direct-RI-RHF", mpqc::scf::DirectRIRHF);
MPQC_CLASS_EXPORT2("PRHF", mpqc::scf::zRHF);
