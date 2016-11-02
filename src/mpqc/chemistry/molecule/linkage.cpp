/*
 * linkage.cpp
 *
 *  Created on: Nov 2, 2016
 *      Author: evaleev
 */

#include "mpqc/chemistry/molecule/molecule.h"
#include "mpqc/util/keyval/forcelink.h"

MPQC_CLASS_EXPORT_KEY2("Molecule", mpqc::Molecule);
MPQC_FORCELINK_KEYVAL_CTOR(mpqc::Molecule)
