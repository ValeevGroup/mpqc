
#ifndef SRC_MPQC_CHEMISTRY_MOLECULE_LINKAGE_H_
#define SRC_MPQC_CHEMISTRY_MOLECULE_LINKAGE_H_

#include "mpqc/util/keyval/forcelink.h"

namespace mpqc {
class Molecule;
mpqc::detail::ForceLink<Molecule> fl1;
namespace molecule {
class PeriodicSystem;
mpqc::detail::ForceLink<PeriodicSystem> fl2;
}
}

#endif  // SRC_MPQC_CHEMISTRY_MOLECULE_LINKAGE_H_
