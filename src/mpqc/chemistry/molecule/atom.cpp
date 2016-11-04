#include "mpqc/chemistry/molecule/atom.h"

#include <map>
#include <iostream>

namespace mpqc {

static std::map<int, std::string> atom_names = {{1, "H"},
                                                {2, "He"},
                                                {3, "Li"},
                                                {4, "Be"},
                                                {5, "B"},
                                                {6, "C"},
                                                {7, "N"},
                                                {8, "O"},
                                                {9, "F"},
                                                {10, "Ne"},
                                                {11, "Na"},
                                                {12, "Mg"}};

// Taken from https://en.wikipedia.org/wiki/Atomic_units on 07/08/15
double bohr_to_ang = 0.52917721092;

std::string Atom::xyz_string(bool convert_to_angstroms) const {
    std::string name = atom_names[atomic_number_];
    name += ' ';

    Vector3d center = center_;
    if (convert_to_angstroms) {
        center *= bohr_to_ang;
    }

    name += (std::to_string(center[0]) + ' ');
    name += (std::to_string(center[1]) + ' ');
    name += std::to_string(center[2]);

    return name;
}

std::ostream &operator<<(std::ostream &os, Atom const &a) {
    os << "Atom: {Z: " << a.charge() << ", mass: " << a.mass()
       << ", pos: " << a.center().transpose() << "}";
    return os;
}

} // namespace mpqc
