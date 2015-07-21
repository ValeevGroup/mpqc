#include "../include/tbb.h"
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

    using iter_t = decltype(elements_.begin());
    return tbb::parallel_reduce(
        tbb::blocked_range<iter_t>(elements_.begin(), elements_.end()), 0.0,
        [&](const tbb::blocked_range<iter_t> &r, double d) {
            return std::accumulate(r.begin(), r.end(), d, reduce_r);
        },
        std::plus<double>());
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
