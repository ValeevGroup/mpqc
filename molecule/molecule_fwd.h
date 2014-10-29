#pragma once
#ifndef MOLECULE_FWD_H
#define MOLECULE_FWD_H

#include "../include/eigen.h"

namespace tcc {
namespace molecule {

class Atom;
class Cluster;
class Clusterable;
class Molecule;
using position_t = Eigen::Vector3d;

namespace clustering {
class kmeans;
}

} // namespace molecule
} // namespace tcc


#endif // MOLECULE_FWD_H
