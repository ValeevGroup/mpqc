//
// Created by Chong Peng on 1/12/17.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_PROPERTIES_LINKAGE_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_PROPERTIES_LINKAGE_H_

#include "mpqc/chemistry/qc/properties/energy.h"
#include "mpqc/util/keyval/forcelink.h"

namespace mpqc{

class Energy;
mpqc::detail::ForceLink<Energy> fl_energy;
template<> class mpqc::MolecularFiniteDifferenceDerivative<1,double>;
mpqc::detail::ForceLink<mpqc::MolecularFiniteDifferenceDerivative<1,double>> fl_fdgrad;

}

#endif //MPQC4_SRC_MPQC_CHEMISTRY_QC_PROPERTIES_LINKAGE_H_
