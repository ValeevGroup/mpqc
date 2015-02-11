#include <memory>
#include <fstream>
#include <algorithm>

#include "include/tbb.h"
#include "include/libint.h"
#include "include/tiledarray.h"
#include "include/btas.h"

#include "utility/make_array.h"
#include "utility/parallel_print.h"

#include "molecule/atom.h"
#include "molecule/cluster.h"
#include "molecule/molecule.h"
#include "molecule/clustering_functions.h"
#include "molecule/make_clusters.h"

#include "basis/atom_basisset.h"
#include "basis/basis_set.h"
#include "basis/cluster_shells.h"
#include "basis/basis.h"

#include "integrals/btas_to_ta_tensor.h"
#include "integrals/btas_to_low_rank_tensor.h"
#include "integrals/make_engine.h"
#include "integrals/ta_tensor_to_low_rank_tensor.h"
#include "integrals/integral_engine_pool.h"
#include "integrals/sparse_task_integrals.h"

#include "purification/sqrt_inv.h"
#include "purification/purification_devel.h"

using namespace tcc;
namespace ints = integrals;

void debug_par(madness::World &world, volatile int debug) {
    if (0 != debug) {
        char hostname[256];
        gethostname(hostname, sizeof(hostname));
        printf("PID %d on %s ready for attach\n", getpid(), hostname);
        fflush(stdout);
        if (world.rank() == 0) {
            while (0 != debug) sleep(5);
        }
    }
    world.gop.fence();
}

int main(int argc, char *argv[]) {
    auto &world = madness::initialize(argc, argv);
    std::string mol_file = "";
    std::string basis_name = "";
    std::string df_basis_name = "";
    int nclusters = 0;
    volatile int debug = 0;
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
    debug_par(world, debug);

    auto mol = molecule::read_xyz(mol_file);

    utility::print_par(world, "Computing ", mol.nelements(), " elements with ",
                       nclusters, " clusters \n");

    auto clusters = molecule::attach_hydrogens_kmeans(mol, nclusters);

    basis::BasisSet bs{basis_name};
    basis::BasisSet df_bs{df_basis_name};

    basis::Basis basis{bs.create_basis(clusters)};
    basis::Basis df_basis{df_bs.create_basis(clusters)};


    libint2::init();
    utility::print_par(world, "Computing overlap\n");
    auto overlap_pool = ints::make_pool(ints::make_1body("overlap", basis));
    auto overlap = integrals::BlockSparseIntegrals(
        world, overlap_pool, utility::make_array(basis, basis),
        integrals::compute_functors::BtasToTaTensor{});

    utility::print_par(world, "Computing eri2\n");
    auto eri_pool = ints::make_pool(ints::make_2body(basis, df_basis));
    auto eri2 = integrals::BlockSparseIntegrals(
        world, eri_pool, utility::make_array(df_basis, df_basis),
        integrals::compute_functors::BtasToTaTensor{});
    std::cout << overlap << std::endl;

    utility::print_par(world, "Computing overlap inverse\n");
    auto overlap_inv_sqrt = pure::inverse_sqrt(overlap);

    utility::print_par(world, "Computing kinetic\n");
    auto kinetic_pool = ints::make_pool(ints::make_1body("kinetic", basis));
    auto T = integrals::BlockSparseIntegrals(
        world, kinetic_pool, utility::make_array(basis, basis),
        integrals::compute_functors::BtasToTaTensor{});

    utility::print_par(world, "Computing nuclear\n");
    auto nuclear_pool = ints::make_pool(ints::make_1body("nuclear", basis));
    auto V = integrals::BlockSparseIntegrals(
        world, nuclear_pool, utility::make_array(basis, basis),
        integrals::compute_functors::BtasToTaTensor{});

    utility::print_par(world, "Computing Hcore\n");
    decltype(V) H;
    H("i,j") = T("i,j") + V("i,j");

    utility::print_par(world, "Computing Density\n");
    auto purifier = pure::make_orthogonal_tr_reset_pure(overlap_inv_sqrt);
    auto D = purifier(H, 10);

    /*
    auto eri2_inv_sqrt_low_rank = TiledArray::conversion::to_new_tile_type(
        eri2, integrals::compute_functors::TaToLowRankTensor<2>{1e-8});

    if (world.rank() == 0) {
        std::cout << "Eri2 inverse square root storage: ";
    }
    auto eri2_inv_sqrt_storage = array_storage(eri2_inv_sqrt_low_rank);

    auto eri3 = integrals::BlockSparseIntegrals(
        world, eri_pool, utility::make_array(df_basis, basis, basis),
        integrals::compute_functors::BtasToTaTensor{});

    if (world.rank() == 0) {
        std::cout << "Eri3 TA::Tensor sparse storage: ";
    }
    auto eri3_storage = array_storage(eri3);

    auto eri3_low_rank = TiledArray::conversion::to_new_tile_type(
        eri3, integrals::compute_functors::TaToLowRankTensor<2>{1e-8});

    if (world.rank() == 0) {
        std::cout << "Eri3 low rank tile storage: ";
    }
    auto eri3_low_rank_storage = array_storage(eri3_low_rank);

    double eri3_diff = array_diff(eri3, eri3_low_rank);
    if (world.rank() == 0) {
        std::cout << "The difference between initial eri3 arrays was "
                  << eri3_diff << std::endl;
    }

    world.gop.fence();

    decltype(eri3_low_rank) Xab;
    Xab("X,a,b") = eri2_inv_sqrt_low_rank("X,P") * eri3_low_rank("P,a,b");
    Xab.truncate();

    world.gop.fence();
    for (auto it = Xab.begin(); it != Xab.end(); ++it) {
        it->get().compress();
    }
    if (world.rank() == 0) {
        std::cout << "Eri3 sqrt inverse storage: ";
    }
    auto eri3_contract_storage = array_storage(Xab);

    eri3("X,a,b") = eri2("X,P") * eri3("P,a,b");
    eri3.truncate();

    double diff = array_diff(eri3, Xab);
    if (world.rank() == 0) {
        std::cout << "The difference between contracted eri3 arrays was "
                  << diff << std::endl;
    }

    */

    world.gop.fence();
    libint2::cleanup();
    madness::finalize();
    return 0;
}
