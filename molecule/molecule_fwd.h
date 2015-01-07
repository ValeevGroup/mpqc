#pragma once
#ifndef TCC_MOLECULE_MOLECULE_FWD_H
#define TCC_MOLECULE_MOLECULE_FWD_H

#include "../include/eigen.h"

namespace tcc {
namespace molecule {

using position_t = Eigen::Vector3d;

class Atom;
class Cluster;
class Clusterable;
class Molecule;

namespace clustering {
class kmeans;
}

} // namespace molecule
} // namespace tcc


#endif // TCC_MOLECULE_MOLECULE_FWD_H
