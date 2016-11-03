//
// Created by Chong Peng on 10/5/16.
//

#ifndef SRC_MPQC_CHEMISTRY_QC_WFN_LINKAGE_H_
#define SRC_MPQC_CHEMISTRY_QC_WFN_LINKAGE_H_

#include "mpqc/util/keyval/forcelink.h"
#include "mpqc/chemistry/molecule/linkage.h"

namespace mpqc {
namespace qc {

class WavefunctionWorld;
mpqc::detail::ForceLink<WavefunctionWorld> fl1;

}
}


#endif  // SRC_MPQC_CHEMISTRY_QC_WFN_LINKAGE_H_
