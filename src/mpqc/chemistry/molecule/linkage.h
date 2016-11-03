
#ifndef SRC_MPQC_CHEMISTRY_MOLECULE_LINKAGE_H_
#define SRC_MPQC_CHEMISTRY_MOLECULE_LINKAGE_H_


#include <mpqc/util/keyval/forcelink.h>

namespace mpqc {

class Molecule;
class UnitCell;

mpqc::detail::ForceLink<mpqc::Molecule> fl1;
mpqc::detail::ForceLink<mpqc::UnitCell> fl2;

}

#endif  // SRC_MPQC_CHEMISTRY_MOLECULE_LINKAGE_H_
