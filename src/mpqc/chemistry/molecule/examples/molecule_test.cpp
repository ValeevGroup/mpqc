#include <chrono>
#include <iostream>
#include <vector>

#include "mpqc/chemistry/molecule/atom.h"
#include "mpqc/chemistry/molecule/molecule.h"

using namespace mpqc;
using namespace mpqc::molecule;

Molecule read_xyz(std::ifstream &f) {
    // Get number of atoms.
    unsigned long natoms = 0;
    f >> natoms;

    std::string line;
    std::vector<Clusterable> clusterables;
    while (std::getline(f, line)) {
        if (!line.empty()) {
            std::stringstream ss(line);
            std::string atom = "";
            double x = 0.0;
            double y = 0.0;
            double z = 0.0;
            ss >> atom;
            ss >> x;
            ss >> y;
            ss >> z;
            if (atom == "H") {
                clusterables.emplace_back(molecule::Atom({x, y, z}, 1, 1));
            } else if (atom == "O") {
                clusterables.emplace_back(molecule::Atom({x, y, z}, 16, 8));
            }
        }
    }
    return Molecule{std::move(clusterables)};
}

int main(int argc, char **argv) {
    std::string mol_file = "";
    int nclusters = 0;
    if(argc == 3){
        mol_file = argv[1];
        nclusters = std::stoi(argv[2]);
    } else {
        std::cout << "Need input file and/or number of clusters\n";
        return 0;
    }
    
    std::ifstream molecule_file(mol_file);
    auto mol = read_xyz(molecule_file);
    molecule_file.close();

    // Making clusters
    auto mc0 = std::chrono::high_resolution_clock::now();
    auto clusters = mol.attach_H_and_kmeans(nclusters);
    auto mc1 = std::chrono::high_resolution_clock::now();
    double mc_alloc =
        std::chrono::duration_cast<std::chrono::duration<double>>(mc1 - mc0)
            .count();
    std::cout << "cluster allocing time = " << mc_alloc << std::endl;

    auto i = 0;
    for(auto &c : clusters){
        auto atoms = collapse_to_atoms(c);
        std::cout << "Cluster " << i << " has " << atoms.size() << " atoms\n";
        ++i;
    }

    return 0;
}
