#include <memory>
#include <fstream>
#include <algorithm>

#include "../common/namespaces.h"
#include "../common/typedefs.h"

#include "../include/tbb.h"
#include "../include/libint.h"
#include "../include/tiledarray.h"
#include "../include/btas.h"

#include "../utility/make_array.h"
#include "../utility/parallel_print.h"
#include "../utility/array_storage.h"
#include "../utility/ta_helpers.h"
#include "../utility/time.h"

#include "../molecule/atom.h"
#include "../molecule/cluster.h"
#include "../molecule/molecule.h"
#include "../molecule/clustering_functions.h"
#include "../molecule/make_clusters.h"

#include "../integrals/ta_tensor_to_low_rank_tensor.h"

#include "../basis/atom_basisset.h"
#include "../basis/basis_set.h"
#include "../basis/cluster_shells.h"
#include "../basis/basis.h"

#include "../integrals/btas_to_ta_tensor.h"
#include "../integrals/integral_engine_pool.h"
#include "../integrals/sparse_task_integrals.h"
#include "../integrals/make_engine.h"

#include "../purification/purification_devel.h"
#include "../purification/sqrt_inv.h"

using namespace tcc;

template <typename T, unsigned int DIM, typename TileType, typename Policy>
void print_size_info(TiledArray::Array<T, DIM, TileType, Policy> const &a,
                     std::string name) {
    utility::print_par(a.get_world(), "Printing size information for ", name,
                       "\n");

    auto data = utility::array_storage(a);

    utility::print_par(a.get_world(), "\tFull   = ", data[0], " GB\n",
                       "\tSparse = ", data[1], " GB\n", "\tLow Rank = ",
                       data[2], " GB\n");
}

int main(int argc, char *argv[]) {
    auto &world = madness::initialize(argc, argv);
    std::string mol_file = "";
    std::string basis_name = "";
    int nclusters = 0;
    double threshold = 1e-11;
    double low_rank_threshold = 1e-7;
    bool use_columb_metric = true;
    if (argc >= 4) {
        mol_file = argv[1];
        basis_name = argv[2];
        nclusters = std::stoi(argv[3]);
    } else {
        std::cout << "input is $./program mol_file basis_file df_basis_file "
                     "nclusters \n";
        return 0;
    }
    if (argc == 5) {
        use_columb_metric = std::stoi(argv[4]);
    }
    if (argc == 6) {
        threshold = std::stod(argv[5]);
    }
    TiledArray::SparseShape<float>::threshold(threshold);
    utility::print_par(world, "Sparse threshold is ",
                       TiledArray::SparseShape<float>::threshold(), "\n");
    utility::print_par(world, "Low Rank threshold is ", low_rank_threshold,
                       "\n");


    auto mol = molecule::read_xyz(mol_file);
    auto repulsion_energy = mol.nuclear_repulsion();

    utility::print_par(world, "Computing ", mol.nelements(), " elements with ",
                       nclusters, " clusters. Nuclear repulsion energy = ",
                       repulsion_energy, "\n");

    auto clusters = molecule::attach_hydrogens_kmeans(mol, nclusters);

    basis::BasisSet bs{basis_name};
    basis::Basis basis{bs.create_basis(clusters)};

    libint2::init();
    utility::print_par(world, "Starting integrals\n");

    TAArray<2, TA::SparsePolicy> eri2;

    if (use_columb_metric) {
        utility::print_par(world, "Using Coulomb Metric\n");
        auto eri_pool = integrals::make_pool(integrals::make_2body(basis));

        eri2 = integrals::BlockSparseIntegrals(
            world, eri_pool, utility::make_array(basis, basis),
            integrals::compute_functors::BtasToTaTensor{});
        world.gop.fence();
    } else {
        utility::print_par(world, "Using Kinetic Metric\n");
        auto eri_pool
            = integrals::make_pool(integrals::make_1body("kinetic", basis));

        eri2 = integrals::BlockSparseIntegrals(
            world, eri_pool, utility::make_array(basis, basis),
            integrals::compute_functors::BtasToTaTensor{});
        world.gop.fence();
    }

    { // Compute low rank information for Eri2
        auto eri2_lr = TA::to_new_tile_type(
            eri2, integrals::compute_functors::TaToLowRankTensor<2>{
                      low_rank_threshold});
        print_size_info(eri2_lr, "Eri2");
    }

    if (world.rank() == 0) {
        std::cout << "Finished Integrals\nTesting inverse sqrt." << std::endl;
    }
    {
        auto inv_sqrt_timer
            = tcc_time::make_timer([&]() { return pure::inverse_sqrt(eri2); });

        auto sqrt_inv = inv_sqrt_timer.apply();
        auto inv_lr = TA::to_new_tile_type(
            sqrt_inv, integrals::compute_functors::TaToLowRankTensor<2>{
                          low_rank_threshold});

        utility::print_par(world, "V^{-1/2} took ", inv_sqrt_timer.time(),
                           " s\n");
        print_size_info(inv_lr, "V^{-1/2}");

        // Check accuracy
        decltype(sqrt_inv) inv;
        inv("i,j") = sqrt_inv("i,k") * sqrt_inv("k,j");
        inv.truncate();
        auto inv_full_lr = TA::to_new_tile_type(
            inv, integrals::compute_functors::TaToLowRankTensor<2>{
                     low_rank_threshold});
        print_size_info(inv_full_lr, "V^{-1}");
        // Reuse inv to compute approximate identity.
        inv("i,j") = inv("i,k") * eri2("k,j");
        inv.truncate();

        auto ident = pure::create_diagonal_matrix(inv, 1.0);

        auto f_norm_diff = utility::array_fnorm_diff(inv, ident);
        utility::print_par(world, "Idenity-(S*S^{-1}) F norm difference = ",
                           f_norm_diff, "\n");
    }

    world.gop.fence();
    libint2::cleanup();
    madness::finalize();
    return 0;
}
