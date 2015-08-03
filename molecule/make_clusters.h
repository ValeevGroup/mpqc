#pragma once
#ifndef TCC_MOLECULE_MAKECLUSTERS_H
#define TCC_MOLECULE_MAKECLUSTERS_H

#include "molecule_fwd.h"
#include <memory>
#include <vector>

namespace tcc {
namespace molecule {

std::vector<std::shared_ptr<Cluster>>
attach_hydrogens_kmeans(Molecule const &, unsigned long);

std::vector<std::shared_ptr<Cluster>>
kmeans(Molecule const &, unsigned long);

} // namespace molecule
} // namespace tcc

#endif // TCC_MOLECULE_MAKECLUSTERS_H
