#include "atom_based_cluster.h"
#include "atom.h"

#include "common.h"

#include <iostream>

namespace mpqc {
namespace molecule {

void AtomBasedCluster::update_center(){
    center_ = Vec3D{0,0,0};

    for(auto const &elem : elements_){
        center_ += elem.center();
    }

    center_ /= elements_.size();
}

double AtomBasedCluster::sum_distances_from_center() const {
    assert(false); // I am not sure this is correct.
    // auto reduce_r = [&](double d, const Clusterable &c) {
    //     return d + std::sqrt(diff_squaredNorm(c.center(), center_));
    // };

    // return std::accumulate(elements_.begin(), elements_.end(), 0.0, reduce_r);
}

std::vector<Atom> AtomBasedCluster::atoms() const {
    std::vector<Atom> atoms;
    for(auto const &elem : elements_){
        std::vector<Atom> temp_atoms = elem.atoms();
        atoms.insert(atoms.end(), temp_atoms.begin(), temp_atoms.end());
    }
    return atoms;
}

std::ostream & operator<<(std::ostream &os, AtomBasedCluster const &c){
    const auto end = c.end();
    const auto last = end - 1;
    os << "AtomBasedCluster: {";
    os << "Center: " << c.center().transpose() << ", elements: {";
    for(auto i = c.begin(); i != end; ++i){
        if(i != last){
            i->print(os) << ", ";
        } else {
            i->print(os) << "}}";
        }
    }

    return os;
}


} // namespace molecule
} // namespace mpqc
