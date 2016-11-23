#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_PHF_LINKAGE_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_PHF_LINKAGE_H_

#include "mpqc/chemistry/qc/wfn/linkage.h"
#include "mpqc/util/keyval/forcelink.h"

namespace mpqc {
namespace phf {

class PHF;
mpqc::detail::ForceLink<PHF> fl1;

}
}

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_PHF_LINKAGE_H_

