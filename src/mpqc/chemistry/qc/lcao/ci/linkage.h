//
// Created by Chong Peng on 3/2/17.
//

#ifndef SRC_MPQC_CHEMISTRY_QC_LCAO_CI_LINKAGE_H_
#define SRC_MPQC_CHEMISTRY_QC_LCAO_CI_LINKAGE_H_

#include <tiledarray.h>
#include "mpqc/chemistry/qc/lcao/scf/linkage.h"
#include "mpqc/util/keyval/forcelink.h"
#include "mpqc/mpqc_config.h"

namespace mpqc{
namespace lcao{

template <typename Tile, typename Policy>
class CIS;

namespace ci{
#if TA_DEFAULT_POLICY == 0
mpqc::detail::ForceLink<CIS<TA::TensorD, TA::DensePolicy>> fl1;
#elif TA_DEFAULT_POLICY == 1
mpqc::detail::ForceLink<CIS<TA::TensorD, TA::SparsePolicy>> fl1;
#endif

} //namespace ci

} // namespace lcao
} // namespace mpqc

#endif //SRC_MPQC_CHEMISTRY_QC_LCAO_CI_LINKAGE_H_
