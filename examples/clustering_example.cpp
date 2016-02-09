#include "../molecule/atom_based_cluster_concept.h"
#include "../molecule/atom_based_cluster.h"

#include "../molecule/molecule.h"
#include "../molecule/clustering_functions.h"

#include "../clustering/kmeans.h"

#include "../include/libint.h"

using namespace mpqc;

int main(int argc, char **argv) {
    mol::AtomBasedCluster h201(mol::Atom({0, -1, 1}, 1, 1),
                               mol::Atom({0, 1, 1}, 1, 1),
                               mol::Atom({0, 0, 0}, 16, 8));
    mol::AtomBasedCluster h202(mol::Atom({2, -1, 1}, 1, 1),
                               mol::Atom({2, 1, 1}, 1, 1),
                               mol::Atom({2, 0, 0}, 16, 8));
    mol::AtomBasedCluster h203(mol::Atom({4, -1, 1}, 1, 1),
                               mol::Atom({4, 1, 1}, 1, 1),
                               mol::Atom({4, 0, 0}, 16, 8));

    std::vector<mol::AtomBasedClusterable> initial_clusters{
          std::move(h201), std::move(h202), std::move(h203)};

    mol::Molecule mol(std::move(initial_clusters));

    std::cout << "Molecule initial = \n" << mol << std::endl;

    clustering::Kmeans kmeans;

    auto clusters
          = kmeans.cluster<mol::AtomBasedCluster>(mol.clusterables(), 2);

    std::cout << "\nClusters from kmeans\n";
    for (auto const &c : clusters) {
        std::cout << c << "\n\n";
    }

    auto mol2 = mol::attach_hydrogens_and_kmeans(mol.clusterables(), 2);

    std::cout << "\nMolecule Clusters\n";
    std::cout << mol2 << std::endl;

    return 0;
}
