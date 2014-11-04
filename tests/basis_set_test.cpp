#include "../basis/basis_set.h"
#include "../basis/atom_basisset.h"
#include "../molecule/atom.h"
#include "../molecule/cluster.h"
#include "gtest.h"
#include "../include/libint.h"

#include <iostream>

TEST(BasisSet, FileConstructor) {
    tcc::basis::BasisSet bs("ano_basis_G94.txt");
}

TEST(BasisSet, CreateLibintShell) {
    tcc::basis::BasisSet bs("ano_basis_G94.txt");
    tcc::molecule::Atom a{{0, 0, 1}, 1, 1};
    auto shells = bs.atom_basis(a);
    for (auto const &sh : shells) {
        std::cout << sh << std::endl;
    }
}
