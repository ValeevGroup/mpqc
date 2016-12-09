//
// Created by Chong Peng on 10/31/16.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_CC_LINKAGE_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_CC_LINKAGE_H_

#include "mpqc/chemistry/qc/scf/linkage.h"
#include "mpqc/util/keyval/forcelink.h"
#include "mpqc/mpqc_config.h"

namespace mpqc {
namespace cc {

template <typename Tile, typename Policy>
class CCSD;
template <typename Tile, typename Policy>
class CCSD_T;
template <typename Tile, typename Policy>
class DBCCSD;

#if TA_DEFAULT_POLICY == 0
mpqc::detail::ForceLink<CCSD<TA::TensorD, TA::DensePolicy>> fl1;
mpqc::detail::ForceLink<CCSD_T<TA::TensorD, TA::DensePolicy>> fl2;
mpqc::detail::ForceLink<DBCCSD<TA::TensorD, TA::DensePolicy>> fl3;
#elif TA_DEFAULT_POLICY == 1
mpqc::detail::ForceLink<CCSD<TA::TensorD, TA::SparsePolicy>> fl1;
mpqc::detail::ForceLink<CCSD_T<TA::TensorD, TA::SparsePolicy>> fl2;
mpqc::detail::ForceLink<DBCCSD<TA::TensorD, TA::SparsePolicy>> fl3;
#endif
}
}

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_CC_LINKAGE_H_
