#pragma once
#ifndef MPQC_MOLECULE_ATTACHHYDROGENS_H
#define MPQC_MOLECULE_ATTACHHYDROGENS_H

#include "molecule_fwd.h"
#include <vector>

namespace mpqc {
namespace molecule {
namespace clustering {

/// Functor that takes clusterables and returns a vector of clusters which have
/// attached all hydrogens to their nearest heavy atom.  If there are only
/// hydrogens then return each hydrogen as its own cluster.  If there are no
/// hydrogens then return each heavy atom as its own cluster.
class attach_hydrogens {
  public:
    /// take a vector of clusterables and return a vector of clusters.
    /// clusterables are copied because they will be sorted.
    std::vector<Cluster> operator()(std::vector<Clusterable> clusterables);
};

} // namespace clustering
} // namespace molecule
} // namespace mpqc


#endif /* end of include guard: MPQC_MOLECULE_ATTACHHYDROGENS_H */
