#pragma once
#ifndef MPQC_MOLECULE_MOLECULE_FWD_H
#define MPQC_MOLECULE_MOLECULE_FWD_H

#include "../common/typedefs.h"
#include "../include/eigen.h"

namespace mpqc {
namespace molecule {

class Atom;
class Cluster;
class Clusterable;
class AtomBasedClusterable;
class Molecule;

namespace clustering {
class kmeans;
}

} // namespace molecule
} // namespace mpqc


#endif // MPQC_MOLECULE_MOLECULE_FWD_H
