#include <memory>
#include <fstream>
#include <algorithm>
#include <iomanip>

#include "include/tbb.h"
#include "include/libint.h"
#include "include/tiledarray.h"
#include "include/btas.h"

#include <TiledArray/algebra/diis.h>

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

    utility::print_par(world, "Computing overlap\n");
    auto overlap_pool = ints::make_pool(ints::make_1body("overlap", basis));
    auto S = integrals::BlockSparseIntegrals(
        world, overlap_pool, utility::make_array(basis, basis),
        integrals::compute_functors::BtasToTaTensor{});

    utility::print_par(world, "Computing overlap inverse\n");
    auto overlap_inv_sqrt = pure::inverse_sqrt(S);

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
    auto D = purifier(H, occupation);

    utility::print_par(world, "Computing eri2\n");
    auto eri_pool = ints::make_pool(ints::make_2body(basis, df_basis));
    auto eri2 = integrals::BlockSparseIntegrals(
        world, eri_pool, utility::make_array(df_basis, df_basis),
        integrals::compute_functors::BtasToTaTensor{});

    utility::print_par(world, "Computing eri2 sqrt Inverse\n");
    auto eri2_sqrt_inv = pure::inverse_sqrt(eri2);

    utility::print_par(world, "Computing eri3\n");
    auto Xab = integrals::BlockSparseIntegrals(
        world, eri_pool, utility::make_array(df_basis, basis, basis),
        integrals::compute_functors::BtasToTaTensor{});

    utility::print_par(world, "Forming the symmetric three center product\n");
    Xab("X,a,b") = eri2_sqrt_inv("X,P") * Xab("P,a,b");

    decltype(D) J, K, F;
    decltype(Xab) Exch;
    J("i,j") = (Xab("X,a,b") * D("b,a")) * Xab("X,i,j");
    // Makes order 4 temp, currently only for testing purposes.
    Exch("i,a,X") = Xab("X,i,b") * D("b,a");
    K("i,j") = Exch("i,a,X") * Xab("X,a,j");
    F("i,j") = H("i,j") + 2 * J("i,j") - K("i,j");

    D = purifier(F, occupation);
    auto energy = D("i,j").dot(F("i,j") + H("i,j"), world).get();

    auto diis = TiledArray::DIIS<decltype(D)>{3, 7};
    auto iter = 1;
    decltype(F) Ferror;
    auto error = 1.0;
    const auto volume = double(F.trange().elements().volume());
    while (error >= 1e-9 && iter <= 100) {
        auto t0 = std::chrono::high_resolution_clock::now();
        J("i,j") = (Xab("X,a,b") * D("b,a")) * Xab("X,i,j");
        Exch("i,a,X") = Xab("X,i,b") * D("b,a");
        K("i,j") = Exch("i,a,X") * Xab("X,a,j");
        F("i,j") = H("i,j") + 2 * J("i,j") - K("i,j");
        Ferror("i,j") = F("i,k") * D("k,l") * S("l,j")
                        - S("i,k") * D("k,l") * F("l,j");
        error = Ferror("i,j").abs_max().get() / volume;
        diis.extrapolate(F, Ferror);

        D("i,j") = purifier(F, occupation)("i,j");
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


    world.gop.fence();
    libint2::cleanup();
    madness::finalize();
    return 0;
}
