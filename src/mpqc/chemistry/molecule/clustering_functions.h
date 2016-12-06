
#ifndef MPQC4_SRC_MPQC_CHEMISTRY_MOLECULE_CLUSTERING_FUNCTIONS_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_MOLECULE_CLUSTERING_FUNCTIONS_H_

#include "mpqc/chemistry/molecule/molecule_fwd.h"

#include <vector>

namespace mpqc {

Molecule attach_hydrogens_and_kmeans(std::vector<AtomBasedClusterable> const &,
                                     int64_t);

Molecule kmeans(std::vector<AtomBasedClusterable> const &, int64_t);

}  // namespace tcc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_MOLECULE_CLUSTERING_FUNCTIONS_H_
