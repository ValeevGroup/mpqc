//
// Created by Chong Peng on 10/5/16.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_LINKAGE_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_LINKAGE_H_

#include "mpqc/chemistry/qc/lcao/wfn/linkage.h"
#include "mpqc/util/keyval/forcelink.h"
#include "mpqc/mpqc_config.h"

namespace mpqc {
namespace lcao {

template <typename Tile, typename Policy>
class RHF;

template <typename Tile, typename Policy>
class RIRHF;

template <typename Tile, typename Policy>
class DirectRHF;

template <typename Tile, typename Policy>
class DirectRIRHF;

template <typename Tile, typename Policy>
class zRHF;

class DFzRHF;

namespace scf {
#if TA_DEFAULT_POLICY == 0
mpqc::detail::ForceLink<RHF<TA::TensorD, TA::DensePolicy>> fl1;
mpqc::detail::ForceLink<RIRHF<TA::TensorD, TA::DensePolicy>> fl2;
mpqc::detail::ForceLink<DirectRHF<TA::TensorD, TA::DensePolicy>> fl3;
mpqc::detail::ForceLink<DirectRIRHF<TA::TensorD, TA::DensePolicy>> fl4;
#elif TA_DEFAULT_POLICY == 1
mpqc::detail::ForceLink<RHF<TA::TensorD, TA::SparsePolicy>> fl1;
mpqc::detail::ForceLink<RIRHF<TA::TensorD, TA::SparsePolicy>> fl2;
mpqc::detail::ForceLink<DirectRHF<TA::TensorD, TA::SparsePolicy>> fl3;
mpqc::detail::ForceLink<DirectRIRHF<TA::TensorD, TA::SparsePolicy>> fl4;
mpqc::detail::ForceLink<zRHF<TA::TensorD, TA::SparsePolicy>> fl5;
mpqc::detail::ForceLink<DFzRHF> fl6;
#endif
}

} // namespace lcao
} // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_LINKAGE_H_
