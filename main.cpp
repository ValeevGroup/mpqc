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
#include "integrals/sparse_task_integrals.h"

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
    tbb::task_scheduler_init(1);
    basis::BasisSet bs{"3-21G_basis_G94.txt"};

    std::vector<molecule::Clusterable> clusterables;
    clusterables.emplace_back(molecule::Atom{{0, -1, 1}, 1, 1});
    clusterables.emplace_back(molecule::Atom{{0, 1, 1}, 1, 1});
    clusterables.emplace_back(molecule::Atom{{0, 0, 0}, 16, 8});

    std::ifstream molecule_file("mol.xyz");
    auto mol = read_xyz(molecule_file);
    molecule::Molecule mol_test{std::move(clusterables)};

    int nclusters = std::max(1, static_cast<int>(mol.nelements() / 9));
    std::cout << mol.nelements() << " atoms  with " << nclusters << " clusters"
              << std::endl;

    auto cluster_func = molecule::clustering::kmeans{127};
    std::vector<std::shared_ptr<molecule::Cluster>> clusters;
    clusters.reserve(nclusters);

    for (auto &&cluster : mol.cluster_molecule(cluster_func, nclusters)) {
        clusters.push_back(
            std::make_shared<molecule::Cluster>(std::move(cluster)));
    }

    basis::Basis basis{bs.create_basis(clusters)};

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

    libint2::init();
    libint2::OneBodyEngine overlap{libint2::OneBodyEngine::overlap, max_nprim,
                                   static_cast<int>(max_am)};

    auto overlap_pool = integrals::make_pool(std::move(overlap));

    auto sparse_S = integrals::SparseIntegrals(
        world, overlap_pool, basis,
        integrals::compute_functors::TaTileFunctor<double>{});
    world.gop.fence();

    if (world.rank() == 0) {
        std::cout << "Computing 4 center integrals now" << std::endl;
    }

    libint2::TwoBodyEngine<libint2::Coulomb> eri{max_nprim,
                                                 static_cast<int>(max_am)};

    auto eri_pool = integrals::make_pool(std::move(eri));
    auto eri4 = integrals::Integrals(world, eri_pool, basis);

    world.gop.fence();
    if (world.rank() == 0) {
        std::cout << "Finished 4 center integrals now" << std::endl;
    }


#if 0
    auto a = pure::inverse_sqrt(sparse_S);
    //    decltype(sparse_S) ident;
    //    ident("i,j") = a("i,k") * sparse_S("k,l") * a("l,j");

    world.gop.fence();
#endif
    madness::finalize();
    return 0;
}
