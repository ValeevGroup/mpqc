#include <mpqc/chemistry/molecule/atom_based_cluster.h>
#include <mpqc/chemistry/molecule/common.h>

namespace mpqc {
namespace molecule {

void AtomBasedCluster::update_cluster() {
    mass_ = sum_mass(elements_);
    charge_ = sum_charge(elements_);
    com_ = center_of_mass(elements_);
}

std::vector<Atom> AtomBasedCluster::atoms() const {
    std::vector<Atom> atoms;
    for (auto const &elem : elements_) {
        std::vector<Atom> temp_atoms = elem.atoms();
        atoms.insert(atoms.end(), temp_atoms.begin(), temp_atoms.end());
    }
    return atoms;
}

std::vector<Atom>
collapse_to_atoms(AtomBasedCluster const &abc) {
    return abc.atoms();
}

std::ostream &operator<<(std::ostream &os, AtomBasedCluster const &c) {
    const auto end = c.end();
    const auto last = end - 1;
    os << "AtomBasedCluster: {";
    os << "C. Of Mass: " << c.com().transpose() << ", elements: {";
    for (auto i = c.begin(); i != end; ++i) {
        if (i != last) {
            i->print(os) << ", ";
        } else {
            i->print(os) << "}}";
        }
    }

    return os;
}


} // namespace molecule
} // namespace mpqc
