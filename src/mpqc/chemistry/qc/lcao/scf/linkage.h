//
// Created by Chong Peng on 10/5/16.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_LINKAGE_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_LINKAGE_H_

#include "mpqc/chemistry/qc/lcao/factory/linkage.h"
#include "mpqc/mpqc_config.h"
#include "mpqc/util/keyval/forcelink.h"

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

template <typename Tile, typename Policy>
class DFzRHF;

namespace scf {

template <typename Tile, typename Policy>
class FosterBoysLocalizer;

template <typename Tile, typename Policy>
class RRQRLocalizer;

#if TA_DEFAULT_POLICY == 0
mpqc::detail::ForceLink<RHF<TA::TensorD, TA::DensePolicy>> fl1;
mpqc::detail::ForceLink<RIRHF<TA::TensorD, TA::DensePolicy>> fl2;
mpqc::detail::ForceLink<DirectRHF<TA::TensorD, TA::DensePolicy>> fl3;
mpqc::detail::ForceLink<DirectRIRHF<TA::TensorD, TA::DensePolicy>> fl4;
mpqc::detail::ForceLink<FosterBoysLocalizer<TA::TensorD, TA::DensePolicy>> fl7;
mpqc::detail::ForceLink<RRQRLocalizer<TA::TensorD, TA::DensePolicy>> fl8;
#elif TA_DEFAULT_POLICY == 1
mpqc::detail::ForceLink<RHF<TA::TensorD, TA::SparsePolicy>> fl1;
mpqc::detail::ForceLink<RIRHF<TA::TensorD, TA::SparsePolicy>> fl2;
mpqc::detail::ForceLink<DirectRHF<TA::TensorD, TA::SparsePolicy>> fl3;
mpqc::detail::ForceLink<DirectRIRHF<TA::TensorD, TA::SparsePolicy>> fl4;
mpqc::detail::ForceLink<zRHF<TA::TensorD, TA::SparsePolicy>> fl5;
mpqc::detail::ForceLink<DFzRHF<TA::TensorD, TA::SparsePolicy>> fl6;
mpqc::detail::ForceLink<FosterBoysLocalizer<TA::TensorD, TA::SparsePolicy>> fl7;
mpqc::detail::ForceLink<RRQRLocalizer<TA::TensorD, TA::SparsePolicy>> fl8;
#endif
}

}  // namespace lcao
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_LINKAGE_H_
