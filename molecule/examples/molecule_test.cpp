#include <iostream>
#include <vector>
#include <chrono>

#include "../Atom.h"
#include "../../include/tbb.h"
#include "../cluster_concept.h"
#include "../molecule.h"
#include "../cluster_collapse.h"

#include "../cluster.h"
#include "../clustering_functions.h"

#include "../../include/libint.h"
#include <fstream>

using namespace tcc;
using namespace tcc::molecule;

molecule::Molecule read_xyz(std::ifstream &f) {
    // Get number of atoms.
    unsigned long natoms = 0;
    f >> natoms;

    std::string line;
    std::vector<molecule::Clusterable> clusterables;
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
    return molecule::Molecule{std::move(clusterables)};
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
    tbb::tick_count mc0 = tbb::tick_count::now();
    auto clusters = mol.attach_H_and_kmeans(nclusters);
    tbb::tick_count mc1 = tbb::tick_count::now();
    double mc_alloc = (mc1 - mc0).seconds();
    std::cout << "cluster allocing time = " << mc_alloc << std::endl;

    auto i = 0;
    for(auto &c : clusters){
        auto atoms = collapse_to_atoms(c);
        std::cout << "Cluster " << i << " has " << atoms.size() << " atoms\n";
        ++i;
    }
            

    /*
    using iter_t = decltype(clusters.begin());
    double sum = tbb::parallel_reduce(
        tbb::blocked_range<iter_t>(clusters.begin(), clusters.end()), 0.0,
        [](const tbb::blocked_range<iter_t> & r, double d)->double {
            return std::accumulate(r.begin(), r.end(), d,
                                   [](double d, const Clusterable
                                      &b) { return d + b.center().norm(); });
        },
        std::plus<double>());

    std::cout << "Sum of distances = " << sum << std::endl;
    */

    return 0;
}
