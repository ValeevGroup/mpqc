#ifndef SRC_MPQC_CHEMISTRY_MOLECULE_LINKAGE_H_
#define SRC_MPQC_CHEMISTRY_MOLECULE_LINKAGE_H_

#include <mpqc/chemistry/molecule/molecule.h>
#include<mpqc/chemistry/molecule/periodic_system.h>
#include <mpqc/util/keyval/forcelink.h>

namespace mpqc {

mpqc::detail::ForceLink<mpqc::Molecule> fl1;
mpqc::detail::ForceLink<mpqc::molecule::PeriodicSystem> fl2;

}

#endif // MPQC_MOLECULE_CLUSTER_H
