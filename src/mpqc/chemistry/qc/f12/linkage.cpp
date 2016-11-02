/*
 * linkage.cpp
 *
 *  Created on: Nov 2, 2016
 *      Author: evaleev
 */

#include "mpqc/chemistry/qc/f12/dbccsdf12.h"
#include "mpqc/chemistry/qc/f12/dbmp2f12.h"
#include "mpqc/chemistry/qc/f12/gf2f12.h"
#include "mpqc/util/keyval/forcelink.h"

MPQC_CLASS_EXPORT2("RMP2F12", mpqc::f12::RMP2F12);
MPQC_CLASS_EXPORT2("RI-RMP2F12", mpqc::f12::RIRMP2F12);
MPQC_CLASS_EXPORT2("RI-DBRMP2F12", mpqc::f12::RIDBRMP2F12);
MPQC_CLASS_EXPORT2("CCSD(F12)", mpqc::f12::CCSDF12<TA::TensorD>);
MPQC_CLASS_EXPORT2("DBCCSD(F12)", mpqc::f12::DBCCSDF12<TA::TensorD>);
MPQC_CLASS_EXPORT2("GF2F12", mpqc::f12::GF2F12<TA::TensorD>);

