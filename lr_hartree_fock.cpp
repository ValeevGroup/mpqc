#include <memory>
#include <fstream>
#include <algorithm>
#include <iomanip>

#include "include/tbb.h"
#include "include/libint.h"
#include "include/tiledarray.h"
#include "include/btas.h"

#include "utility/make_array.h"
#include "utility/parallel_print.h"
#include "utility/array_storage.h"
#include "utility/function_timer.h"

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
#include "integrals/dense_task_integrals.h"

#include "purification/sqrt_inv.h"
#include "purification/purification_devel.h"

#include <chrono>

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

template <typename Pool, typename Basis, unsigned long DIM, typename Fn>
TiledArray::Array<double, DIM, typename Fn::TileType, TiledArray::SparsePolicy>
time_and_print_block_sparse(madness::World &world, Pool &&pool,
                            std::array<Basis, DIM> bs, Fn fn,
                            std::string name) {

    utility::print_par(world, name, "\n");
    auto timer = utility::make_timer([=, &world]() {
        return integrals::BlockSparseIntegrals(world, pool, bs, fn);
    });
    auto result = timer.apply();
    utility::print_par(world, name, " time = ", timer.time(), "\n");
    print_size_info(result, name);
    return result;
}

int main(int argc, char *argv[]) {
    auto &world = madness::initialize(argc, argv);
    std::string mol_file = "";
    std::string basis_name = "";
    std::string df_basis_name = "";
    int nclusters = 0;
    double threshold = 1e-7;
    double low_rank_threshold = 1e-9;
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
        threshold = std::stod(argv[5]);
    }
    if (argc == 7) {
        debug = std::stoi(argv[6]);
    }
    debug_par(world, debug);
    TiledArray::SparseShape<float>::threshold(threshold);
    utility::print_par(world, "Sparse threshold is ",
                       TiledArray::SparseShape<float>::threshold(), "\n");

    auto mol = molecule::read_xyz(mol_file);
    auto charge = 0;
    auto occupation = mol.occupation(charge);
    auto repulsion_energy = mol.nuclear_repulsion();

    utility::print_par(world, "Computing ", mol.nelements(), " elements with ",
                       nclusters, " clusters. Nuclear repulsion energy = ",
                       repulsion_energy, "\n");

    auto clusters = molecule::attach_hydrogens_kmeans(mol, nclusters);

    basis::BasisSet bs{basis_name};
    basis::BasisSet df_bs{df_basis_name};

    basis::Basis basis{bs.create_basis(clusters)};
    basis::Basis df_basis{df_bs.create_basis(clusters)};


    libint2::init();

    // Compute overlap.
    auto overlap_pool = ints::make_pool(ints::make_1body("overlap", basis));
    auto S = time_and_print_block_sparse(
        world, overlap_pool, utility::make_array(basis, basis),
        integrals::compute_functors::BtasToTaTensor{}, "Overlap");

    // Invert overlap
    utility::print_par(world, "\nComputing overlap inverse\n");
    auto overlap_inv_sqrt = pure::inverse_sqrt(S);
    print_size_info(S, "Overlap inverse sqrt");

    // Compute T
    auto kinetic_pool = ints::make_pool(ints::make_1body("kinetic", basis));
    utility::print_par(world, "\n");
    auto T = time_and_print_block_sparse(
        world, kinetic_pool, utility::make_array(basis, basis),
        integrals::compute_functors::BtasToTaTensor{}, "Kinetic");

    // Compute V
    auto nuclear_pool = ints::make_pool(ints::make_1body("nuclear", basis));
    utility::print_par(world, "\n");
    auto V = time_and_print_block_sparse(
        world, nuclear_pool, utility::make_array(basis, basis),
        integrals::compute_functors::BtasToTaTensor{}, "Nuclear");

    // Compute Hcore
    utility::print_par(world, "\nComputing Hcore\n");
    decltype(V) H;
    H("i,j") = T("i,j") + V("i,j");
    print_size_info(H, "Hcore");

    // Compute intial density
    utility::print_par(world, "\nComputing Density\n");
    auto purifier = pure::make_orthogonal_tr_reset_pure(overlap_inv_sqrt);
    auto D = purifier(H, occupation);
    print_size_info(D, "D initial");

    // Begin Two electron integrals section.
    auto eri_pool = ints::make_pool(ints::make_2body(basis, df_basis));

    // Computing Eri2
    utility::print_par(world, "\n");
    auto eri2 = time_and_print_block_sparse(
        world, eri_pool, utility::make_array(df_basis, df_basis),
        integrals::compute_functors::BtasToTaTensor{}, "Eri2 Integrals");

    // Computing the sqrt inverse of Eri2
    utility::print_par(world, "\nComputing eri2 sqrt Inverse\n");
    auto inv_timer
        = utility::make_timer([&]() { return pure::inverse_sqrt(eri2); });
    auto eri2_sqrt_inv = inv_timer.apply();
    utility::print_par(world, "Eri2 inverse computation time = ",
                       inv_timer.time(), "\n");
    print_size_info(eri2_sqrt_inv, "Eri2 sqrt inverse");

    /*
     * Start using Low rank arrays
     */

    // Convert 1/sqrt(V) into a low rank tensor.
    auto eri2_inv_lr = TiledArray::conversion::to_new_tile_type(
        eri2_sqrt_inv,
        integrals::compute_functors::TaToLowRankTensor<2>{low_rank_threshold});

    // Compute center integrals
    utility::print_par(world, "\n");
    auto Xab = time_and_print_block_sparse(
        world, eri_pool, utility::make_array(df_basis, basis, basis),
        integrals::compute_functors::BtasToLowRankTensor{low_rank_threshold},
        "Eri3 integrals");

    // Form symetric order 3 tensor
    utility::print_par(world, "\nForming the symmetric three center product\n");
    auto Xabt = utility::make_timer(
        [&]() { Xab("X,a,b") = eri2_inv_lr("X,P") * Xab("P,a,b"); });
    Xabt.apply();
    Xab.truncate();
    utility::print_par(world, "Eri3X contraction time = ", Xabt.time(), "\n");
    print_size_info(Xab, "Eri3 * 1/sqrt(V)");

    // Convert D to low rank
    utility::print_par(world, "\nCreating Exchange Temp.");
    auto D_lr = TiledArray::conversion::to_new_tile_type(
        D,
        integrals::compute_functors::TaToLowRankTensor<2>{low_rank_threshold});

    // Compute the exchange intermediate tensor
    decltype(Xab) Exch;
    Exch("i,X,a") = Xab("X,i,b") * D_lr("b,a");
    print_size_info(Exch, "Exchange temp.");


    /*
    decltype(D) J, K, F;
    J("i,j") = (Xab("X,a,b") * D("a,b")) * Xab("X,i,j");
    K("i,j") = (Xab("X,i,b") * D("b,a")) * Xab("X,a,j");
    F("i,j") = H("i,j") + 2 * J("i,j") - K("i,j");

    D = purifier(F, occupation);
    auto energy = D("i,j").dot(F("i,j") + H("i,j"), world).get();

    utility::print_par(world, "\nStarting SCF iterations\n");
    auto diis = TiledArray::DIIS<decltype(D)>{3, 7};
    auto iter = 1;
    decltype(F) Ferror;
    auto error = 1.0;
    const auto volume = double(F.trange().elements().volume());
    while (error >= 1e-13 && iter <= 35) {
        auto t0 = std::chrono::high_resolution_clock::now();
        J("i,j") = (Xab("X,a,b") * D("a,b")) * Xab("X,i,j");
        J.truncate();
        auto k0 = std::chrono::high_resolution_clock::now();
        K("i,j") = (Xab("X,i,b") * D("b,a")) * Xab("X,a,j");
        auto k1 = std::chrono::high_resolution_clock::now();
        F("i,j") = H("i,j") + 2 * J("i,j") - K("i,j");
        Ferror("i,j") = F("i,k") * D("k,l") * S("l,j")
                        - S("i,k") * D("k,l") * F("l,j");
        Ferror.truncate();
        error = Ferror("i,j").norm().get() / volume;
        diis.extrapolate(F, Ferror);

        D = purifier(F, occupation);
        energy = D("i,j").dot(F("i,j") + H("i,j"), world).get();
        world.gop.fence();
        auto t1 = std::chrono::high_resolution_clock::now();
        auto time = std::chrono::duration_cast<std::chrono::duration<double>>(
                        t1 - t0).count();
        auto ktime = std::chrono::duration_cast<std::chrono::duration<double>>(
                         k1 - k0).count();
        utility::print_par(world, "Iteration: ", iter++, " has energy ",
                           std::setprecision(11), energy, " with error ", error,
                           " in ", time, " s with K time ", ktime, "\n");
    }


    utility::print_par(world, "Final energy = ", std::setprecision(11),
                       energy + repulsion_energy, "\n");

   */
    world.gop.fence();
    libint2::cleanup();
    madness::finalize();
    return 0;
}
