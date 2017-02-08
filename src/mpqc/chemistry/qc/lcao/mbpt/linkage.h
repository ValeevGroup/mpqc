//
// Created by Chong Peng on 10/6/16.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_MBPT_LINKAGE_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_MBPT_LINKAGE_H_

#include "mpqc/chemistry/qc/lcao/scf/linkage.h"
#include "mpqc/util/keyval/forcelink.h"
#include "mpqc/mpqc_config.h"

namespace mpqc {
namespace lcao {
template<typename Tile, typename Policy>
class RMP2;

template<typename Tile, typename Policy>
class RIRMP2;

template<typename Tile, typename Policy>
class DBRMP2;

template<typename Tile, typename Policy>
class RIDBRMP2;

//template<typename Tile, typename Policy>
//class GammaPointMP2;
template<typename Tile, typename Policy>
class GammaPointMP2;

namespace mbpt {
#if TA_DEFAULT_POLICY == 0
mpqc::detail::ForceLink<RMP2<TA::TensorD, TA::DensePolicy>> fl1;
mpqc::detail::ForceLink<RIRMP2<TA::TensorD, TA::DensePolicy>> fl2;
mpqc::detail::ForceLink<DBRMP2<TA::TensorD, TA::DensePolicy>> fl3;
mpqc::detail::ForceLink<RIDBRMP2<TA::TensorD, TA::DensePolicy>> fl4;
mpqc::detail::ForceLink<GammaPointMP2<TA::TensorD, TA::DensePolicy>> fl5;

#elif TA_DEFAULT_POLICY == 1
mpqc::detail::ForceLink<RMP2<TA::TensorD, TA::SparsePolicy>> fl1;
mpqc::detail::ForceLink<RIRMP2<TA::TensorD, TA::SparsePolicy>> fl2;
mpqc::detail::ForceLink<DBRMP2<TA::TensorD, TA::SparsePolicy>> fl3;
mpqc::detail::ForceLink<RIDBRMP2<TA::TensorD, TA::SparsePolicy>> fl4;
mpqc::detail::ForceLink<GammaPointMP2<TA::TensorD, TA::SparsePolicy>> fl5;

#endif
}

} // namespace lcao
} // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_MBPT_LINKAGE_H_
