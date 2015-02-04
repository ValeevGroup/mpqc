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
#include "integrals/btas_to_low_rank_tensor.h"
#include "integrals/ta_tensor_to_low_rank_tensor.h"
#include "integrals/integral_engine_pool.h"
#include "integrals/sparse_task_integrals.h"

#include "purification/purification_devel.h"
#include "purification/sqrt_inv.h"

using namespace tcc;

struct dummy_converter {

    template <typename T>
    TiledArray::Tensor<double> operator()(T const &) {
        return TiledArray::Tensor<double>{};
    }
};

int main(int argc, char *argv[]) {
    auto &world = madness::initialize(argc, argv);
    std::string mol_file = "";
    std::string basis_name = "";
    std::string df_basis_name = "";
    int nclusters = 0;
    int debug = 0;
    if (argc >= 5) {
        mol_file = argv[1];
        basis_name = argv[2];
        df_basis_name = argv[3];
        nclusters = std::stoi(argv[4]);
    } else {
        std::cout << "input is $./program mol_file basis_file df_basis_file "
                     "nclusters \n";
        return 0;
    }
    if (argc == 6) {
        debug = std::stoi(argv[5]);
    }

    auto mol = molecule::read_xyz(mol_file);

    basis::BasisSet bs{basis_name};
    basis::BasisSet df_bs{df_basis_name};

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
    world.gop.fence();
    if (world.rank() == 0) {
        std::cout << "Computing Integrals\n";
    }

    libint2::TwoBodyEngine<libint2::Coulomb> eri{max_nprim,
                                                 static_cast<int>(max_am)};


    auto eri_pool = integrals::make_pool(std::move(eri));
    /*
    auto eri3 = integrals::BlockSparseIntegrals(
        world, eri_pool, utility::make_array(df_basis, basis, basis),
        integrals::compute_functors::BtasToLowRankTensor{1e-8});

    world.gop.fence();

    if (world.rank() == 0) {
        std::cout << "Finished 3 center integrals." << std::endl;
    }
    */

    auto eri2 = integrals::BlockSparseIntegrals(
        world, eri_pool, utility::make_array(df_basis, df_basis),
        integrals::compute_functors::BtasToTaTensor{});

    world.gop.fence();

    if (world.rank() == 0) {
        std::cout << "Finished Integrals\n" << std::endl;
    }

    auto eri2_inv_sqrt = pure::inverse_sqrt(eri2);
    auto eri2_inv_sqrt_low_rank = TiledArray::conversion::to_new_tile_type(
        eri2_inv_sqrt, integrals::compute_functors::TaToLowRankTensor<2>{1e-8});

    world.gop.fence();

    if (world.rank() == 0) {
        std::cout << "Finished inverting and converting 2 center intgrals"
                     "\nStarting space counting" << std::endl;
    }
    double full_storage = 0.0;
    double low_rank_storage = 0.0;
    for (auto it = eri2_inv_sqrt_low_rank.begin();
         it != eri2_inv_sqrt_low_rank.end(); ++it) {
        const auto tensor = it->get();
        if (tensor.isFull()) {
            low_rank_storage += tensor.tile().ftile().size();
            full_storage += tensor.tile().ftile().size();
        } else {
            std::cout << "There was a low rank tile." << std::endl;
            low_rank_storage += tensor.tile().lrtile().matrixL().size();
            low_rank_storage += tensor.tile().lrtile().matrixR().size();
            full_storage += tensor.tile().lrtile().matrixL().rows()
                            * tensor.tile().lrtile().matrixR().cols();
        }
    }

    full_storage *= sizeof(double) * 1e-9;
    low_rank_storage *= sizeof(double) * 1e-9;
    world.gop.fence();

    if (world.rank() == 0) {
        std::cout << "WARNING FOLLOWING CURRENTLY ONLY VALID ON ONE "
                     "NODE!!!!!!!!!!!!!!!" << std::endl;
        std::cout << "Full storage = " << full_storage
                  << " GB low rank storage = " << low_rank_storage << " GB"
                  << std::endl;
    }

    /*
    decltype(eri3) Xab{eri3.get_world(), eri3.trange(), eri3.get_shape()};
    Xab("X,a,b") = eri2_inv_sqrt_low_rank("X,P") * eri3("P,a,b");
    */

    world.gop.fence();
    if (world.rank() == 0) {
        std::cout << "Finished contracting 2 and 3 center ints" << std::endl;
    }

    world.gop.fence();
    libint2::cleanup();
    madness::finalize();
    return 0;
}
