#include <memory>
#include <fstream>
#include <algorithm>

#include "include/tbb.h"
#include "include/libint.h"
#include "include/tiledarray.h"
#include "include/btas.h"

#include "utility/make_array.h"

#include "molecule/atom.h"
#include "molecule/cluster.h"
#include "molecule/molecule.h"
#include "molecule/clustering_functions.h"

#include "basis/atom_basisset.h"
#include "basis/basis_set.h"
#include "basis/cluster_shells.h"
#include "basis/basis.h"

#include "integrals/btas_to_ta_tensor.h"
#include "integrals/integral_engine_pool.h"
#include "integrals/sparse_task_integrals.h"

#include "purification/purification_devel.h"
#include "purification/sqrt_inv.h"

using namespace tcc;

molecule::Molecule read_xyz(std::ifstream &f) {
    // Get number of atoms.
    unsigned long natoms = 0;
    f >> natoms;

    constexpr auto ang_to_bohr = 1.0 / 0.52917721092;

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
            x *= ang_to_bohr;
            y *= ang_to_bohr;
            z *= ang_to_bohr;
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
        std::cout << mol.nelements() << " elements with " << nclusters
                  << " clusters" << std::endl;
    }

    std::vector<std::shared_ptr<molecule::Cluster>> clusters;
    clusters.reserve(nclusters);

    for (auto &&cluster : mol.attach_H_and_kmeans(nclusters)) {
        clusters.push_back(
            std::make_shared<molecule::Cluster>(std::move(cluster)));
    }

    basis::Basis basis{bs.create_basis(clusters)};
    basis::Basis df_basis{df_bs.create_basis(clusters)};

    auto max_am = std::max(basis.max_am(), df_basis.max_am());
    auto max_nprim = std::max(basis.max_nprim(), df_basis.max_nprim());

    libint2::init();
    libint2::TwoBodyEngine<libint2::Coulomb> eri{max_nprim,
                                                 static_cast<int>(max_am)};

    auto eri_pool = integrals::make_pool(std::move(eri));
    auto eri2 = integrals::BlockSparseIntegrals(
        world, eri_pool, utility::make_array(df_basis, df_basis),
        integrals::compute_functors::BtasToTaTensor{});

    auto sqrt_inv = pure::inverse_sqrt(eri2);
    {
        decltype(sqrt_inv) ident;
        ident("i,j") = sqrt_inv("i,k") * eri2("k,l") * sqrt_inv("l,j");
        const auto id_norm = double{ident("a,b").norm()};
        const auto real_norm = std::sqrt(ident.trange().elements().size()[0]);
        std::cout << "Ident norm = " << id_norm << " correct = " << real_norm
                  << std::endl;
    }

    /*
    auto eri3 = integrals::BlockSparseIntegrals(
        world, eri_pool, utility::make_array(df_basis, basis, basis),
        integrals::compute_functors::BtasToTaTensor{});
    eri3("P,u,v") = sqrt_inv("P,X")*eri3("X,u,v");
    */

    world.gop.fence();
    if (world.rank() == 0) {
        std::cout << "Finished integrals." << std::endl;
    }
    libint2::cleanup();
    madness::finalize();
    return 0;
}
