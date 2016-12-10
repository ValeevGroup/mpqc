//
// Created by Chong Peng on 10/31/16.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_CC_LINKAGE_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_CC_LINKAGE_H_

#include "mpqc/chemistry/qc/scf/linkage.h"
#include "mpqc/util/keyval/forcelink.h"

namespace mpqc {
namespace lcao {

template <typename Tile, typename Policy>
class CCSD;
template <typename Tile, typename Policy>
class CCSD_T;
template <typename Tile, typename Policy>
class DBCCSD;

namespace cc {
mpqc::detail::ForceLink<CCSD<TA::TensorD, TA::SparsePolicy>> fl1;
mpqc::detail::ForceLink<CCSD_T<TA::TensorD, TA::SparsePolicy>> fl2;
mpqc::detail::ForceLink<DBCCSD<TA::TensorD, TA::SparsePolicy>> fl3;
}  // namespace cc
}  // namespace lcao
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_CC_LINKAGE_H_
