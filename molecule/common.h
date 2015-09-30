#pragma once
#ifndef MPQC_MOLECULE_COMMON_H
#define MPQC_MOLECULE_COMMON_H

#include "cluster_concept.h"
#include "molecule_fwd.h"
#include <vector>

namespace mpqc {
namespace molecule {

/// Function takes mpqc::molecule::Atom vector and converts it to a vector 
/// of libint atoms. 
std::vector<libint2::Atom> to_libint_atom(std::vector<Atom> const &);

} // namespace molecule
} // namespace mpqc

#endif // MPQC_MOLECULE_COMMON_H
