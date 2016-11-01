//
// Created by Chong Peng on 10/10/16.
//

#ifndef MPQC_CHEMISTRY_QC_F12_LINKAGE_H
#define MPQC_CHEMISTRY_QC_F12_LINKAGE_H

#include "dbmp2f12.h"
#include "dbccsdf12.h"
#include <mpqc/util/keyval/forcelink.h>

namespace mpqc{
namespace f12{

mpqc::detail::ForceLink<RMP2F12> fl1;
mpqc::detail::ForceLink<RIRMP2F12> fl2;
mpqc::detail::ForceLink<RIDBRMP2F12> fl3;

mpqc::detail::ForceLink<CCSDF12<TA::TensorD>> fl4;
mpqc::detail::ForceLink<DBCCSDF12<TA::TensorD>> fl5;
}
}

#endif //MPQC_CHEMISTRY_QC_F12_LINKAGE_H
