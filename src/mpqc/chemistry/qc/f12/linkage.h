//
// Created by Chong Peng on 10/10/16.
//

#ifndef SRC_MPQC_CHEMISTRY_QC_F12_LINKAGE_H_
#define SRC_MPQC_CHEMISTRY_QC_F12_LINKAGE_H_

#include "mpqc/chemistry/qc/cc/linkage.h"
#include "mpqc/util/keyval/forcelink.h"

namespace mpqc {
namespace f12 {
class RMP2F12;
class RIRMP2F12;
class RIDBRMP2F12;
mpqc::detail::ForceLink<RMP2F12> fl1;
mpqc::detail::ForceLink<RIRMP2F12> fl2;
mpqc::detail::ForceLink<RIDBRMP2F12> fl3;

template <typename Tile>
class CCSDF12;
template <typename Tile>
class DBCCSDF12;
mpqc::detail::ForceLink<CCSDF12<TA::TensorD>> fl4;
mpqc::detail::ForceLink<DBCCSDF12<TA::TensorD>> fl5;

template <typename Tile>
class GF2F12;
mpqc::detail::ForceLink<GF2F12<TA::TensorD>> fl6;
}
}

#endif  // SRC_MPQC_CHEMISTRY_QC_F12_LINKAGE_H_
