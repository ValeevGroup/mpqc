//
// Created by Chong Peng on 10/5/16.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_WFN_LINKAGE_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_WFN_LINKAGE_H_

#include "mpqc/util/keyval/forcelink.h"
#include "mpqc/chemistry/molecule/linkage.h"

namespace mpqc {
namespace lcao {

class WavefunctionWorld;
mpqc::detail::ForceLink<WavefunctionWorld> fl1;

}
}


#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_WFN_LINKAGE_H_
