#include "../molecule/atom_based_cluster.h"

using namespace mpqc;

int main(int argc, char **argv) {
    Atom a1({0, 0, 0}, 1, 1);
    Atom a2({0, 0, 1}, 1, 1);
    Atom a3({0, 0, 2}, 1, 1);
    Atom a4({0, 0, 3}, 2, 1);

    std::cout << "Atom 1 = " << a1 << std::endl;
    std::cout << "Atom 2 = " << a2 << std::endl;
    std::cout << "Atom 3 = " << a3 << std::endl;
    std::cout << "Atom 4 = " << a4 << std::endl;

    // Add clusterable construction
    AtomBasedCluster c12;
    c12.add_clusterable(std::move(a1));
    c12.add_clusterable(std::move(a2));

    // Varadic template construction
    AtomBasedCluster c34(std::move(a3), std::move(a4));

    c12.update_cluster();
    c34.update_cluster();

    // We can print clusters
    std::cout << "\nCluster 12 = " << c12 << std::endl;
    std::cout << "Cluster 34 = " << c34 << std::endl;

    // We can compute centers for clusters

    std::cout << "\nCluster 12 center = " << c12.com().transpose()
              << std::endl;
    std::cout << "Cluster 34 center = " << c34.com().transpose()
              << std::endl;

    // Get masses
    std::cout << "\nCluster 12 mass = " << c12.mass() << std::endl;
    std::cout << "Cluster 34 mass = " << c34.mass() << std::endl;

    // Get charges
    std::cout << "\nCluster 12 charge = " << c12.charge() << std::endl;
    std::cout << "Cluster 34 charge = " << c34.charge() << std::endl;

    // Get atoms
    std::cout << "\nCluster 12 atoms = \n";
    for (auto const &atom : c12.atoms()) {
        std::cout << atom << std::endl;
    }
    std::cout << "\nCluster 34 atoms = \n";
    for (auto const &atom : c34.atoms()) {
        std::cout << atom << std::endl;
    }

    // We can nest clusters
    AtomBasedCluster nc1234(std::move(c12), std::move(c34));
    nc1234.update_cluster();
    std::cout << "\nNested Clusters 1 2 3 4 = " << nc1234 << std::endl;

    // Finally get its center
    std::cout << "\nNested Clusters 1 2 3 4 center = "
              << nc1234.com().transpose() << std::endl;

    // Get masses
    std::cout << "Nested Cluster 1 2 3 4 mass = " << nc1234.mass() << std::endl;
    std::cout << "Nested Cluster 1 2 3 4 charge = " << nc1234.charge() << std::endl;

    std::cout << "\nNested Cluster 1 2 3 4 atoms = \n";
    for (auto const &atom : nc1234.atoms()) {
        std::cout << atom << std::endl;
    }

    return 0;
}
