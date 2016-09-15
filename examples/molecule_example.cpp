#include "../molecule/atom_based_cluster.h"

#include "../molecule/molecule.h"

using namespace mpqc;
using namespace molecule;

int main(int argc, char **argv) {
    AtomBasedCluster h20f(Atom({0, -1, 1}, 1, 1), Atom({0, 1, 1}, 1, 1),
                          Atom({0, 0, 0}, 16, 8));
    AtomBasedCluster h20s(Atom({2, 1, 3}, 1, 1), Atom({2, 3, 3}, 1, 1),
                          Atom({2, 2, 2}, 16, 8));

    std::vector<AtomBasedClusterable> clusters{std::move(h20f),
                                               std::move(h20s)};

    std::cout << "Input atoms: " << std::endl;
    for (auto const &c : clusters) {
        std::cout << "\tCluster Center: " << center(c).transpose() << std::endl;
        for (auto const &a : c.atoms()) {
            std::cout << "\t\t" << a << std::endl;
        }
    }
    std::cout << std::endl;

    Molecule mol(std::move(clusters));
    std::cout << "Molecule center of mass " << mol.com().transpose()
              << std::endl;
    std::cout << "Atoms as seen by molecule:\n";
    for (const auto &a : mol.atoms()) {
        std::cout << a << " ";
        std::cout << "distance from COM: " << (a.center() - mol.com()).norm()
                  << std::endl;
    }
    std::cout << std::endl;

    // Molecule methods
    std::cout << "Molecule nclusters = " << mol.nclusters() << std::endl;
    std::cout << "Molecule charge = " << mol.charge() << std::endl;
    std::cout << "Molecule mass = " << mol.mass() << std::endl;

    // Molecule nuclear repulsion
    std::cout << "Molecule repulsion = " << mol.nuclear_repulsion()
              << std::endl;
    std::cout << "Molecule core_electrons = " << mol.core_electrons()
              << std::endl;

    return 0;
}
