//
// Created by Chong Peng on 10/6/16.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_MBPT_LINKAGE_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_MBPT_LINKAGE_H_

#include "mpqc/chemistry/qc/scf/linkage.h"
#include "mpqc/util/keyval/forcelink.h"

namespace mpqc {
namespace lcao {
class RMP2;
class RIRMP2;
class DBRMP2;
class RIDBRMP2;

namespace mbpt {
mpqc::detail::ForceLink<RMP2> fl1;
mpqc::detail::ForceLink<RIRMP2> fl2;
mpqc::detail::ForceLink<DBRMP2> fl3;
mpqc::detail::ForceLink<RIDBRMP2> fl4;
}  // namespace mbpt
}  // namespace lcao
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_MBPT_LINKAGE_H_
