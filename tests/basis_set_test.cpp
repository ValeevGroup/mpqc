#include "../basis/basis_set.h"
#include "../basis/atom_basisset.h"
#include "../molecule/atom.h"
#include "../molecule/cluster.h"
#include "../basis/cluster_shells.h"
#include "gtest.h"
#include "../include/libint.h"

#include <iostream>

TEST(BasisSet, FileConstructor) {
    tcc::basis::BasisSet bs("3-21G_basis_G94.txt");
    auto abs = bs.atom_basis_set();

    EXPECT_FALSE(abs[0].is_SP(0));
    EXPECT_FALSE(abs[0].shell(0).spherical());
    EXPECT_EQ(abs[0].shell(0).contraction_length(), 1) << abs[0].shell(0);
    EXPECT_DOUBLE_EQ(abs[0].exponent(0, 0), 5.4471780);
    EXPECT_DOUBLE_EQ(abs[0].coeff(0, 0), 0.1562850);
}

TEST(BasisSet, CreateLibintShellShells) {
    tcc::basis::BasisSet bs("3-21G_basis_G94.txt");
    tcc::molecule::Atom a{{0, 0, 1}, 1, 1};
    auto shells = bs.atom_basis(a);
    EXPECT_FALSE(shells[0].contr[0].pure);
    EXPECT_TRUE(shells[0].ncontr() == 1);
    EXPECT_DOUBLE_EQ(shells[0].alpha[0], 5.4471780);
    EXPECT_DOUBLE_EQ(shells[0].contr[0].coeff[0], 0.1562850);
    EXPECT_DOUBLE_EQ(shells[0].O[0], 0.0);
    EXPECT_DOUBLE_EQ(shells[0].O[1], 0.0);
    EXPECT_DOUBLE_EQ(shells[0].O[2], 1.0);
}

TEST(BasisSet, CreateBinnedLibintShells) {
    tcc::basis::BasisSet bs("3-21G_basis_G94.txt");
    tcc::molecule::Atom h{{0, 0, 0}, 1, 1};
    tcc::molecule::Atom f{{0, 0, 1}, 19, 9};
    tcc::molecule::Atom h2{{0, 0, 2}, 1, 1};
    tcc::molecule::Atom f2{{0, 0, 3}, 19, 9};

    auto clust1 = std::make_shared<tcc::molecule::Cluster>();
    auto clust2 = std::make_shared<tcc::molecule::Cluster>();

    clust1->add_clusterable(std::move(h));
    clust1->add_clusterable(std::move(f));
    clust2->add_clusterable(std::move(h2));
    clust2->add_clusterable(std::move(f2));

    auto clustered_shells = bs.create_basis({clust1, clust2});
    //TODO test clustered shells
}
