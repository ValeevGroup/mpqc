//
// Created by Chong Peng on 10/6/16.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_MBPT_LINKAGE_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_MBPT_LINKAGE_H_

#include "mpqc/chemistry/qc/scf/linkage.h"
#include "mpqc/util/keyval/forcelink.h"

namespace mpqc {
namespace mbpt {
template<typename Tile, typename Policy>
class RMP2;

template<typename Tile, typename Policy>
class RIRMP2;

template<typename Tile, typename Policy>
class DBRMP2;

template<typename Tile, typename Policy>
class RIDBRMP2;

mpqc::detail::ForceLink<RMP2<TA::TensorD, TA::SparsePolicy>> fl1;
mpqc::detail::ForceLink<RIRMP2<TA::TensorD, TA::SparsePolicy>> fl2;
mpqc::detail::ForceLink<DBRMP2<TA::TensorD, TA::SparsePolicy>> fl3;
mpqc::detail::ForceLink<RIDBRMP2<TA::TensorD, TA::SparsePolicy>> fl4;
}
}

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_MBPT_LINKAGE_H_
