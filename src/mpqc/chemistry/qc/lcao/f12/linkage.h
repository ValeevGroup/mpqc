//
// Created by Chong Peng on 10/10/16.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_F12_LINKAGE_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_F12_LINKAGE_H_

#include "mpqc/chemistry/qc/lcao/cc/linkage.h"
#include "mpqc/chemistry/qc/lcao/mbpt/linkage.h"
#include "mpqc/util/keyval/forcelink.h"
#include "mpqc/mpqc_config.h"

namespace mpqc {
namespace lcao {

template <typename Tile>
class RMP2F12;

template <typename Tile>
class RIRMP2F12;

template <typename Tile>
class RIDBRMP2F12;

template <typename Tile>
class CCSD_F12;

template <typename Tile>
class DBCCSD_F12;

template <typename Tile>
class CCSD_T_F12;

template <typename Tile>
class GF2F12;

namespace f12{
#if TA_DEFAULT_POLICY == 1
mpqc::detail::ForceLink<RMP2F12<TA::TensorD>> fl1;
mpqc::detail::ForceLink<RIRMP2F12<TA::TensorD>> fl2;
//mpqc::detail::ForceLink<RIDBRMP2F12<TA::TensorD>> fl3;
mpqc::detail::ForceLink<CCSD_F12<TA::TensorD>> fl4;
//mpqc::detail::ForceLink<DBCCSD_F12<TA::TensorD>> fl5;
mpqc::detail::ForceLink<CCSD_T_F12<TA::TensorD>> fl6;
//mpqc::detail::ForceLink<GF2F12<TA::TensorD>> fl7;
#endif
}  // namespace
}  // namespace lcao
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_F12_LINKAGE_H_
