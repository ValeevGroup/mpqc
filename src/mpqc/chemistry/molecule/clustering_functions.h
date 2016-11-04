
#ifndef MPQC_CLUSTERING_FUNCTIONS_H
#define MPQC_CLUSTERING_FUNCTIONS_H

#include "mpqc/chemistry/molecule/molecule_fwd.h"

#include <vector>

namespace mpqc {

Molecule attach_hydrogens_and_kmeans(std::vector<AtomBasedClusterable> const &,
                                     int64_t);

Molecule kmeans(std::vector<AtomBasedClusterable> const &, int64_t);

}  // namespace tcc

#endif  // MPQC_CLUSTERING_FUNCTIONS_H
