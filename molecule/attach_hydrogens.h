#pragma once
#ifndef TCC_MOLECULE_ATTACHHYDROGENS_H
#define TCC_MOLECULE_ATTACHHYDROGENS_H

#include "molecule_fwd.h"
#include <vector>

namespace tcc {
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
} // namespace tcc


#endif /* end of include guard: TCC_MOLECULE_ATTACHHYDROGENS_H */
