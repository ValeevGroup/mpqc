//
// Created by Chong Peng on 10/5/16.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_LINKAGE_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_LINKAGE_H_

#include "mpqc/chemistry/qc/wfn/linkage.h"
#include "mpqc/util/keyval/forcelink.h"

namespace mpqc {
namespace scf {

template <typename Tile, typename Policy>
class RHF;

template <typename Tile, typename Policy>
class RIRHF;

template <typename Tile, typename Policy>
class DirectRHF;

template <typename Tile, typename Policy>
class DirectRIRHF;

class zRHF;

mpqc::detail::ForceLink<RHF<TA::TensorD, TA::SparsePolicy>> fl1;
mpqc::detail::ForceLink<RIRHF<TA::TensorD, TA::SparsePolicy>> fl2;
mpqc::detail::ForceLink<DirectRHF<TA::TensorD, TA::SparsePolicy>> fl3;
mpqc::detail::ForceLink<DirectRIRHF<TA::TensorD, TA::SparsePolicy>> fl4;
mpqc::detail::ForceLink<zRHF> fl5;
}
}

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_LINKAGE_H_
