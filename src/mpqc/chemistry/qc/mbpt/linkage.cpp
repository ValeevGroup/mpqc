/*
 * linkage.cpp
 *
 *  Created on: Nov 2, 2016
 *      Author: evaleev
 */

#include "mpqc/chemistry/qc/mbpt/mp2.h"
#include "mpqc/util/keyval/forcelink.h"

MPQC_FORCELINK_KEYVAL_CTOR(mpqc::mbpt::RMP2)
MPQC_FORCELINK_KEYVAL_CTOR(mpqc::mbpt::RIRMP2)
MPQC_CLASS_EXPORT_KEY2("RMP2", mpqc::mbpt::RMP2);
MPQC_CLASS_EXPORT_KEY2("RI-RMP2", mpqc::mbpt::RIRMP2);
