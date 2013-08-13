/*
 * dmm_scf.cc
 *
 *  Created on: Jul 22, 2013
 *      Author: drewlewis
 */

#include "dmm_scf.h"
#include "k-means.h"

using namespace sc;
using namespace dmm;

DMMHF::DMMHF(const Ref<KeyVal> &keyval) :
{

}

DMMHF::~DMMHF(){}

DMMHF::print(std::ostream&o) const {
    o << "Density Matrix Minimization Closed Shell Hartree-Fock :" << endl;
}

