#pragma once
#ifndef TCC_MOLECULE_COMMON_H
#define TCC_MOLECULE_COMMON_H

#include "cluster_concept.h"
#include "molecule_fwd.h"
#include <vector>

// Fwd declare Libint2 Atom
namespace libint2 {
    struct Atom;
};

namespace tcc {
namespace molecule {
/**
 * center_of_mass returns an position vector which is the center of mass
 * of the input types with centers.
 */
position_t center_of_mass(const std::vector<Clusterable> cs, double mass);

inline double diff_squaredNorm(position_t const &a, position_t const &b) {
    double out = 0.0;
    for (unsigned int i = 0; i < 3; ++i) {
        double temp = a[i] - b[i];
        out += temp * temp;
    }
    return out;
}

/// Function takes tcc::molecule::Atom vector and converts it to a vector 
/// of libint atoms. 
std::vector<libint2::Atom> to_libint_atom(std::vector<Atom> const &);

} // namespace molecule
} // namespace tcc

#endif // TCC_MOLECULE_COMMON_H
