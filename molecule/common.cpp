#include "common.h"
#include "../include/libint.h"
#include "cluster_concept.h"

#include <algorithm>
#include <vector>


namespace tcc {
namespace molecule {

position_t center_of_mass(const std::vector<Clusterable> cs, double mass) {

    position_t com = {0,0,0};
    for(auto const &c : cs){
        com += c.center() * c.mass();
    }
    return com/mass;
}

std::vector<libint2::Atom> to_libint_atom(std::vector<Atom> const &atoms){
    std::vector<libint2::Atom> l_atoms; 
    l_atoms.reserve(atoms.size());

    for(auto const &atom : atoms){
        const auto &c = atom.center();

        libint2::Atom l_atom; 
        l_atom.atomic_number = atom.charge();
        l_atom.x = c[0];
        l_atom.y = c[1];
        l_atom.z = c[2];

        l_atoms.push_back(std::move(l_atom));
    }

    return l_atoms; 
}

} // namespace molecule
} // namespace tcc