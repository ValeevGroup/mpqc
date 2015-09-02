#include "cluster.h"
#include "common.h"

#include <functional>
#include <cmath>
#include <iostream>

namespace tcc {
namespace molecule {

void Cluster::compute_com() { center_ = center_of_mass(elements_, mass_); }

double Cluster::sum_distances_from_center() const {
    auto reduce_r = [&](double d, const Clusterable &c) {
        return d + std::sqrt(diff_squaredNorm(c.center(), center_));
    };

    return std::accumulate(elements_.begin(), elements_.end(), 0.0, reduce_r);
}

std::ostream & operator<<(std::ostream &os, Cluster const &c){
    auto atoms = collapse_to_atoms(c);

    // xyz format specifies number of atoms on first line.
    os << atoms.size() << std::endl;
    os << std::endl; // comment line
    for(auto const & atom : atoms){
        os << atom.xyz_string() << std::endl;
    }

    return os;
}


} // namespace molecule
} // namespace tcc
