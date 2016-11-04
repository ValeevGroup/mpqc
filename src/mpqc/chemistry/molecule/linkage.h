
#ifndef SRC_MPQC_CHEMISTRY_MOLECULE_LINKAGE_H_
#define SRC_MPQC_CHEMISTRY_MOLECULE_LINKAGE_H_


#include <mpqc/util/keyval/forcelink.h>

namespace mpqc {

class Molecule;
class UnitCell;

mpqc::detail::ForceLink<Molecule> fl1;
mpqc::detail::ForceLink<UnitCell> fl2;

}  // namespace mpqc

#endif  // SRC_MPQC_CHEMISTRY_MOLECULE_LINKAGE_H_
