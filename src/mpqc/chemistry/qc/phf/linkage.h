#ifndef MPQC_CHEMISTRY_QC_PHF_LINKAGE_H
#define MPQC_CHEMISTRY_QC_PHF_LINKAGE_H

#include "mpqc/chemistry/qc/wfn/linkage.h"
#include "mpqc/util/keyval/forcelink.h"

namespace mpqc {
namespace phf {

class PHF;
mpqc::detail::ForceLink<PHF> fl1;

}
}

#endif  // MPQC_CHEMISTRY_QC_PHF_LINKAGE_H

