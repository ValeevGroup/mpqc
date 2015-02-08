#include "common.h"
#include "../include/tbb.h"
#include "../include/libint.h"
#include "cluster_concept.h"

#include <algorithm>
#include <vector>


namespace tcc {
namespace molecule {

position_t center_of_mass(const std::vector<Clusterable> cs, double mass) {

    return tbb::parallel_reduce(
               tbb::blocked_range<unsigned long>(0, cs.size()),
               position_t({0, 0, 0}),
               [&](const tbb::blocked_range<unsigned long> & r, position_t p)
                   ->position_t {
                   auto i = r.begin();
                   const auto end = r.end();
                   for (; i != end; ++i) {
                       p += cs[i].center() * cs[i].mass();
                   }
                   return p;
               },
               std::plus<position_t>()) / mass;
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
