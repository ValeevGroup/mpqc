//
// Created by Chong Peng on 10/31/16.
//

#ifndef MPQC_CHEMISTRY_QC_CC_LINKAGE_H
#define MPQC_CHEMISTRY_QC_CC_LINKAGE_H

//#include <mpqc/util/keyval/forcelink.h>
#include "ccsd_t.h"
#include "dbccsd.h"

namespace mpqc{
namespace cc{


mpqc::detail::ForceLink<CCSD<TA::TensorD, TA::SparsePolicy>> fl1;
mpqc::detail::ForceLink<CCSD_T<TA::TensorD,TA::SparsePolicy>> fl2;
mpqc::detail::ForceLink<DBCCSD<TA::TensorD,TA::SparsePolicy>> fl3;

}
}

#endif //MPQC_CHEMISTRY_QC_CC_LINKAGE_H
