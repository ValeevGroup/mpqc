//
// Created by Chong Peng on 1/12/17.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_PROPERTIES_LINKAGE_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_PROPERTIES_LINKAGE_H_

#include "mpqc/util/keyval/forcelink.h"

namespace mpqc{

class Energy;
detail::ForceLink<Energy> fl_energy;
template <size_t Order, typename Value> class MolecularFiniteDifferenceDerivative;
detail::ForceLink<MolecularFiniteDifferenceDerivative<1,double>> fl_fdgrad;

}

#endif //MPQC4_SRC_MPQC_CHEMISTRY_QC_PROPERTIES_LINKAGE_H_
