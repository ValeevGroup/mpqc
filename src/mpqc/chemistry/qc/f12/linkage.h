//
// Created by Chong Peng on 10/10/16.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_F12_LINKAGE_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_F12_LINKAGE_H_

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
class CCSD_F12;
mpqc::detail::ForceLink<CCSD_F12<TA::TensorD>> fl4;

template <typename Tile>
class DBCCSD_F12;
mpqc::detail::ForceLink<DBCCSD_F12<TA::TensorD>> fl5;

template <typename Tile>
class CCSD_T_F12;
mpqc::detail::ForceLink<CCSD_T_F12<TA::TensorD>> fl6;

template <typename Tile>
class GF2F12;
mpqc::detail::ForceLink<GF2F12<TA::TensorD>> fl7;
}
}

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_F12_LINKAGE_H_
