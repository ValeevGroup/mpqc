//
// Created by Chong Peng on 10/31/16.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_CC_LINKAGE_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_CC_LINKAGE_H_

#include "mpqc/chemistry/qc/lcao/scf/linkage.h"
#include "mpqc/util/keyval/forcelink.h"
#include "mpqc/mpqc_config.h"

namespace mpqc {
namespace lcao {

template <typename Tile, typename Policy>
class CCSD;
template <typename Tile, typename Policy>
class CCSD_T;
template <typename Tile, typename Policy>
class GammaPointCCSD;
template <typename Tile, typename Policy>
class EOM_CCSD;
template <typename Tile, typename Policy>
class CCSDT1;
template <typename Tile, typename Policy>
class CCSDT1b;
template <typename Tile, typename Policy>
class CCSDT2;
template <typename Tile, typename Policy>
class CCSDT3;
template <typename Tile, typename Policy>
class CCSDT4;

namespace cc {
#if TA_DEFAULT_POLICY == 0
mpqc::detail::ForceLink<CCSD<TA::TensorD, TA::DensePolicy>> fl1;
mpqc::detail::ForceLink<CCSD_T<TA::TensorD, TA::DensePolicy>> fl2;
mpqc::detail::ForceLink<EOM_CCSD<TA::TensorD,TA::DensePolicy>> fl5;
mpqc::detail::ForceLink<CCSDT1<TA::TensorD,TA::DensePolicy>> fl6;
mpqc::detail::ForceLink<CCSDT1b<TA::TensorD,TA::DensePolicy>> fl7;
mpqc::detail::ForceLink<CCSDT2<TA::TensorD,TA::DensePolicy>> fl8;
mpqc::detail::ForceLink<CCSDT3<TA::TensorD,TA::DensePolicy>> fl9;
mpqc::detail::ForceLink<CCSDT4<TA::TensorD,TA::DensePolicy>> fl10;
#elif TA_DEFAULT_POLICY == 1
mpqc::detail::ForceLink<CCSD<TA::TensorD, TA::SparsePolicy>> fl1;
mpqc::detail::ForceLink<CCSD_T<TA::TensorD, TA::SparsePolicy>> fl2;
mpqc::detail::ForceLink<GammaPointCCSD<TA::TensorD, TA::SparsePolicy>> fl4;
mpqc::detail::ForceLink<EOM_CCSD<TA::TensorD,TA::SparsePolicy>> fl3;
mpqc::detail::ForceLink<CCSDT1<TA::TensorD,TA::SparsePolicy>> fl5;
mpqc::detail::ForceLink<CCSDT1b<TA::TensorD,TA::SparsePolicy>> fl6;
mpqc::detail::ForceLink<CCSDT2<TA::TensorD,TA::SparsePolicy>> fl7;
mpqc::detail::ForceLink<CCSDT3<TA::TensorD,TA::SparsePolicy>> fl8;
mpqc::detail::ForceLink<CCSDT4<TA::TensorD,TA::SparsePolicy>> fl9;
#endif
}  // namespace
}  // namespace lcao
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_CC_LINKAGE_H_
