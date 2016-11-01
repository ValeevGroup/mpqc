//
// Created by Chong Peng on 10/31/16.
//

#ifndef MPQC_CHEMISTRY_QC_CC_LINKAGE_H
#define MPQC_CHEMISTRY_QC_CC_LINKAGE_H

//#include <mpqc/util/keyval/forcelink.h>
#include "ccsd_t.h"

namespace mpqc{
namespace cc{


mpqc::detail::ForceLink<CCSD<TA::TensorD>> fl1;
mpqc::detail::ForceLink<CCSD_T> fl2;

}
}

#endif //MPQC_CHEMISTRY_QC_CC_LINKAGE_H
