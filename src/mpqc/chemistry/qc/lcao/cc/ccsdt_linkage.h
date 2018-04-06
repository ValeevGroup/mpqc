//
// Created by Chong Peng on 10/31/16.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_CC_CCSDT_LINKAGE_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_CC_CCSDT_LINKAGE_H_

#include "mpqc/chemistry/qc/lcao/scf/linkage.h"
#include "mpqc/mpqc_config.h"
#include "mpqc/util/keyval/forcelink.h"

namespace mpqc {
namespace lcao {

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
template <typename Tile, typename Policy>
class CCSDT;
template <typename Tile, typename Policy>
class CC3;

namespace cc {
#if TA_DEFAULT_POLICY == 0
mpqc::detail::ForceLink<CCSDT1<TA::TensorD, TA::DensePolicy>> fl6;
mpqc::detail::ForceLink<CCSDT1b<TA::TensorD, TA::DensePolicy>> fl7;
mpqc::detail::ForceLink<CCSDT2<TA::TensorD, TA::DensePolicy>> fl8;
mpqc::detail::ForceLink<CCSDT3<TA::TensorD, TA::DensePolicy>> fl9;
mpqc::detail::ForceLink<CCSDT4<TA::TensorD, TA::DensePolicy>> fl10;
mpqc::detail::ForceLink<CCSDT<TA::TensorD, TA::DensePolicy>> fl11;
mpqc::detail::ForceLink<CC3<TA::TensorD, TA::DensePolicy>> fl12;
#elif TA_DEFAULT_POLICY == 1
mpqc::detail::ForceLink<CCSDT1<TA::TensorD, TA::SparsePolicy>> fl5;
mpqc::detail::ForceLink<CCSDT1b<TA::TensorD, TA::SparsePolicy>> fl6;
mpqc::detail::ForceLink<CCSDT2<TA::TensorD, TA::SparsePolicy>> fl7;
mpqc::detail::ForceLink<CCSDT3<TA::TensorD, TA::SparsePolicy>> fl8;
mpqc::detail::ForceLink<CCSDT4<TA::TensorD, TA::SparsePolicy>> fl9;
mpqc::detail::ForceLink<CCSDT<TA::TensorD, TA::SparsePolicy>> fl10;
mpqc::detail::ForceLink<CC3<TA::TensorD, TA::SparsePolicy>> fl11;
#endif
}  // namespace cc
}  // namespace lcao
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_CC_CCSDT_LINKAGE_H_
