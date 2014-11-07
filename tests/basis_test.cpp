#include "../basis/basis_set.h"
#include "../basis/basis.h"
#include "../basis/atom_basisset.h"
#include "../molecule/atom.h"
#include "../molecule/cluster.h"
#include "../basis/cluster_shells.h"
#include "gtest.h"
#include "../include/libint.h"
#include "../include/tiledarray.h"

#include <iostream>

TEST(Basis, ClusterShellsConstructor) {

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

    tcc::basis::Basis base { bs.create_basis({clust1, clust2}) };
}

TEST(Basis, TrangeRange1Creation) {

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

    tcc::basis::Basis base { bs.create_basis({clust1, clust2}) };
    std::cout << base.create_trange1() << std::endl;
    std::cout << base.create_flattend_trange1() << std::endl;
}
