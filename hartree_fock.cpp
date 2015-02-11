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

template <typename Array>
void print_size_info(Array const &a, std::string name) {
    utility::print_par(a.get_world(), "Printing size information for ", name,
                       "\n");
    auto data = utility::array_storage(a);
    utility::print_par(a.get_world(), "\tFull   = ", data[0], " GB\n",
                       "\tSparse = ", data[1], " GB\n", "\tLow =    ", data[2],
                       " GB\n");
}


int main(int argc, char *argv[]) {
    auto &world = madness::initialize(argc, argv);
    std::string mol_file = "";
    std::string basis_name = "";
    std::string df_basis_name = "";
    int nclusters = 0;
    double threshold = 1e-7;
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

    utility::print_par(world, "\nComputing overlap\n");
    auto overlap_pool = ints::make_pool(ints::make_1body("overlap", basis));
    auto S = integrals::BlockSparseIntegrals(
        world, overlap_pool, utility::make_array(basis, basis),
        integrals::compute_functors::BtasToTaTensor{});
    print_size_info(S, "Overlap");

    utility::print_par(world, "\nComputing overlap inverse\n");
    auto overlap_inv_sqrt = pure::inverse_sqrt(S);
    print_size_info(S, "Overlap inverse sqrt");

    utility::print_par(world, "\nComputing kinetic\n");
    auto kinetic_pool = ints::make_pool(ints::make_1body("kinetic", basis));
    auto T = integrals::BlockSparseIntegrals(
        world, kinetic_pool, utility::make_array(basis, basis),
        integrals::compute_functors::BtasToTaTensor{});

    utility::print_par(world, "\nComputing nuclear\n");
    auto nuclear_pool = ints::make_pool(ints::make_1body("nuclear", basis));
    auto V = integrals::BlockSparseIntegrals(
        world, nuclear_pool, utility::make_array(basis, basis),
        integrals::compute_functors::BtasToTaTensor{});

    utility::print_par(world, "Computing Hcore\n");
    decltype(V) H;
    H("i,j") = T("i,j") + V("i,j");
    print_size_info(H, "Hcore");

    utility::print_par(world, "\nComputing Density\n");
    auto purifier = pure::make_orthogonal_tr_reset_pure(overlap_inv_sqrt);
    auto D = purifier(H, occupation);
    print_size_info(D, "D initial");

    utility::print_par(world, "\nComputing eri2\n");
    auto eri_pool = ints::make_pool(ints::make_2body(basis, df_basis));
    auto eri2 = integrals::BlockSparseIntegrals(
        world, eri_pool, utility::make_array(df_basis, df_basis),
        integrals::compute_functors::BtasToTaTensor{});
    print_size_info(eri2, "Eri2");

    utility::print_par(world, "\nComputing eri2 sqrt Inverse\n");
    auto eri2_sqrt_inv = pure::inverse_sqrt(eri2);
    print_size_info(eri2_sqrt_inv, "Eri2 sqrt inverse");

    utility::print_par(world, "\nComputing eri3\n");
    auto Xab = integrals::BlockSparseIntegrals(
        world, eri_pool, utility::make_array(df_basis, basis, basis),
        integrals::compute_functors::BtasToTaTensor{});
    print_size_info(Xab, "Eri3");

    utility::print_par(world, "\nForming the symmetric three center product\n");
    Xab("X,a,b") = eri2_sqrt_inv("X,P") * Xab("P,a,b");
    Xab.truncate();
    print_size_info(Xab, "Eri3 * 1/sqrt(V)");

    utility::print_par(world, "\nCreating intial F matrix\n");
    decltype(D) J, K, F;
    decltype(Xab) Exch;
    J("i,j") = (Xab("X,a,b") * D("b,a")) * Xab("X,i,j");
    Exch("i,a,X") = Xab("X,i,b") * D("b,a");
    Exch.truncate();
    print_size_info(Exch, "Exchange intermediate");
    K("i,j") = Exch("i,a,X") * Xab("X,a,j");
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
        J("i,j") = (Xab("X,a,b") * D("b,a")) * Xab("X,i,j");
        J.truncate();
        Exch("i,a,X") = Xab("X,i,b") * D("b,a");
        Exch.truncate();
        K("i,j") = Exch("i,a,X") * Xab("X,a,j");
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
        utility::print_par(world, "Iteration: ", iter++, " has energy ",
                           std::setprecision(11), energy, " with error ", error,
                           " in ", time, " s\n");
    }


    utility::print_par(world, "Final energy = ", std::setprecision(11),
                       energy + repulsion_energy, "\n");

    print_size_info(Exch, "Final exchange temp");

    world.gop.fence();
    libint2::cleanup();
    madness::finalize();
    return 0;
}
