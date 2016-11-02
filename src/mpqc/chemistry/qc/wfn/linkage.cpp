/*
 * linkage.cpp
 *
 *  Created on: Nov 2, 2016
 *      Author: evaleev
 */

#include "mpqc/util/keyval/forcelink.h"
#include "mpqc/chemistry/qc/wfn/wfn_world.h"

MPQC_FORCELINK_KEYVAL_CTOR(mpqc::qc::WavefunctionWorld)
MPQC_CLASS_EXPORT_KEY2("WfnWorld", mpqc::qc::WavefunctionWorld);

