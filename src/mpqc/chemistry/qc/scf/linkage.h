//
// Created by Chong Peng on 10/5/16.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_LINKAGE_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_LINKAGE_H_

#include "mpqc/chemistry/qc/wfn/linkage.h"
#include "mpqc/util/keyval/forcelink.h"

namespace mpqc {
namespace scf {

class RHF;
class RIRHF;
class DirectRHF;
class DirectRIRHF;
class PRHF;
mpqc::detail::ForceLink<RHF> fl1;
mpqc::detail::ForceLink<RIRHF> fl2;
mpqc::detail::ForceLink<DirectRHF> fl3;
mpqc::detail::ForceLink<DirectRIRHF> fl4;
mpqc::detail::ForceLink<PRHF> fl5;
}
}

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_LINKAGE_H_
