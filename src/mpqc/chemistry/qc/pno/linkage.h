/*
 *  Created on: Jan 6, 2017
 *      Author: jinmei
 */

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_PNO_CCSD_PNO_LINKAGE_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_PNO_CCSD_PNO_LINKAGE_H_

#include "mpqc/chemistry/qc/scf/linkage.h"
#include "mpqc/util/keyval/forcelink.h"
#include "mpqc/mpqc_config.h"

namespace mpqc {
namespace lcao {

template <typename Tile, typename Policy>
class CCSD_PNO;

namespace pno {
#if TA_DEFAULT_POLICY == 0
mpqc::detail::ForceLink<CCSD_PNO<TA::TensorD, TA::DensePolicy>> fl1;
#elif TA_DEFAULT_POLICY == 1
mpqc::detail::ForceLink<CCSD_PNO<TA::TensorD, TA::SparsePolicy>> fl1;
#endif
}  // namespace
}  // namespace lcao
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_PNO_CCSD_PNO_LINKAGE_H_
