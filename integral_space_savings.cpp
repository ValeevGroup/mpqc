#include <memory>
#include <fstream>
#include <algorithm>

#include "include/tbb.h"
#include "include/libint.h"
#include "include/tiledarray.h"
#include "include/btas.h"

#include "molecule/atom.h"
#include "molecule/cluster.h"
#include "molecule/molecule.h"
#include "molecule/clustering_functions.h"

#include "basis/atom_basisset.h"
#include "basis/basis_set.h"
#include "basis/cluster_shells.h"
#include "basis/basis.h"

#include "integrals/ta_compute_functors.h"
#include "integrals/integral_engine_pool.h"
#include "integrals/task_integrals.h"
#include "integrals/ta_sparse_space_calculator.h"

#include "purification/purification_devel.h"
#include "purification/sqrt_inv.h"

using namespace tcc;

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

int main(int argc, char *argv[]) {
    auto &world = madness::initialize(argc, argv);
    std::string mol_file = "";
    std::string basis_file = "";
    std::string df_basis_file = "";
    int nclusters = 0;
    if (argc == 5) {
        mol_file = argv[1];
        basis_file = argv[2];
        df_basis_file = argv[3];
        nclusters = std::stoi(argv[4]);
    } else {
        std::cout << "input is $./program mol_file basis_file df_basis_file "
                     "nclusters \n";
        return 0;
    }

    std::ifstream molecule_file(mol_file);
    auto mol = read_xyz(molecule_file);

    basis::BasisSet bs{basis_file};
    basis::BasisSet df_bs{df_basis_file};

    if (world.rank() == 0) {
        std::cout << mol.nelements() << " atoms  with " << nclusters
                  << " clusters" << std::endl;
    }


    std::vector<std::shared_ptr<molecule::Cluster>> clusters;
    clusters.reserve(nclusters);

    for (auto &&cluster : mol.attach_H_and_kmeans(nclusters)) {
        clusters.push_back(
            std::make_shared<molecule::Cluster>(std::move(cluster)));
    }
    if (world.rank() == 0) {
        auto i = 0;
        for (auto const &cluster : clusters) {
            auto atoms = collapse_to_atoms(*cluster);
            std::cout << "Cluster " << i << " has " << atoms.size()
                      << " atoms\n";
            ++i;
        }
    }

    basis::Basis basis{bs.create_basis(clusters)};
    basis::Basis df_basis{df_bs.create_basis(clusters)};

    auto max_nprim = 0ul;
    auto max_am = 0ul;

    for (auto const &cluster : basis.cluster_shells()) {
        auto const &shell_vec = cluster.flattened_shells();

        auto temp_nprim = std::max_element(shell_vec.begin(), shell_vec.end(),
                                           [](libint2::Shell const &s,
                                              libint2::Shell const &t) {
                                               return t.nprim() > s.nprim();
                                           })->nprim();

        auto temp_am = cluster.max_am();

        max_nprim = (temp_nprim > max_nprim) ? temp_nprim : max_nprim;
        max_am = (temp_am > max_am) ? temp_am : max_am;
    }
    for (auto const &cluster : df_basis.cluster_shells()) {
        auto const &shell_vec = cluster.flattened_shells();

        auto temp_nprim = std::max_element(shell_vec.begin(), shell_vec.end(),
                                           [](libint2::Shell const &s,
                                              libint2::Shell const &t) {
                                               return t.nprim() > s.nprim();
                                           })->nprim();

        auto temp_am = cluster.max_am();

        max_nprim = (temp_nprim > max_nprim) ? temp_nprim : max_nprim;
        max_am = (temp_am > max_am) ? temp_am : max_am;
    }

    libint2::init();
    if (world.rank() == 0) {
        std::cout << "Computing 3 center integrals now" << std::endl;
    }

    libint2::TwoBodyEngine<libint2::Coulomb> eri{max_nprim,
                                                 static_cast<int>(max_am)};

    auto eri_pool = integrals::make_pool(std::move(eri));
    integrals::space_calculation::Compute_Storage(
        world, eri_pool, integrals::compute_functors::TaTileFunctor<double>{},
        df_basis, basis, basis);

    world.gop.fence();
    if (world.rank() == 0) {
        std::cout << "Finished 3 center integrals now" << std::endl;
    }

    madness::finalize();
    return 0;
}
