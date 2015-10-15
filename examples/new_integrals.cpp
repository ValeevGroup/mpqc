#include "../common/namespaces.h"
#include "../common/typedefs.h"

#include "../include/libint.h"
#include "../include/tiledarray.h"

#include "../utility/make_array.h"

#include "../clustering/kmeans.h"

#include "../molecule/atom.h"
#include "../molecule/atom_based_cluster.h"
#include "../molecule/molecule.h"
#include "../molecule/clustering_functions.h"
#include "../molecule/make_clusters.h"

#include "../basis/atom_basisset.h"
#include "../basis/basis_set.h"
#include "../basis/basis.h"

#include "../integrals/integrals.h"

#include "../utility/time.h"
#include "../utility/array_storage.h"

#include <memory>

using namespace mpqc;

int main(int argc, char *argv[]) {

    auto &world = madness::initialize(argc, argv);
    std::string mol_file = "";
    std::string basis_name = "";
    int nclusters = 0;
    double threshold = 1e-10;
    if (argc == 4) {
        mol_file = argv[1];
        basis_name = argv[2];
        nclusters = std::stoi(argv[3]);
    } else {
        std::cout << "input is $./program mol_file basis_file nclusters ";
        return 0;
    }
    TiledArray::SparseShape<float>::threshold(threshold);

    auto mol = molecule::read_xyz(mol_file);
    auto clustered_mol
          = molecule::attach_hydrogens_and_kmeans(mol.clusterables(),
                                                  nclusters);

    basis::BasisSet bs(basis_name);
    basis::Basis basis(bs.get_cluster_shells(clustered_mol));

    libint2::init();

    auto eri_pool
          = tcc::integrals::make_pool(tcc::integrals::make_2body(basis));

    auto ta_pass_through =
          [](TA::TensorD &&ten) { return TA::TensorD(std::move(ten)); };

    { // Overlap ints
        if (world.rank() == 0) {
            std::cout << "Overlap ints\n";
        }
        world.gop.fence();
        auto t0 = tcc_time::now();
        auto overlap_pool
              = tints::make_pool(tints::make_1body("overlap", basis));
        auto S = mpqc_ints::sparse_integrals(
              world, overlap_pool, tcc::utility::make_array(basis, basis),
              ta_pass_through);
        auto S_norm = S("i,j").norm(world).get();
        world.gop.fence();
        auto t1 = tcc_time::now();

        if (world.rank() == 0) {
            std::cout << "\tnorm of ints was " << S_norm << std::endl;
            std::cout << "\tin " << tcc_time::duration_in_s(t0, t1)
                      << " seconds" << std::endl;
        }
    }

    { // Two electron Three center
        if (world.rank() == 0) {
            std::cout << "\nTwo E Three Center direct unscreened ints\n";
        }
        world.gop.fence();
        auto t0 = tcc_time::now();

        auto eri_pool
              = tcc::integrals::make_pool(tcc::integrals::make_2body(basis));

        auto eri3 = mpqc_ints::direct_sparse_integrals(
              world, eri_pool, tcc::utility::make_array(basis, basis, basis),
              ta_pass_through);

        auto eri3_norm = eri3("i,j,k").norm(world).get();
        world.gop.fence();
        auto t1 = tcc_time::now();

        if (world.rank() == 0) {
            std::cout << "\tnorm of ints was " << eri3_norm << std::endl;
            std::cout << "\tin " << tcc_time::duration_in_s(t0, t1)
                      << " seconds" << std::endl;
        }
    }

    { // Two electron Three center
        if (world.rank() == 0) {
            std::cout << "Two E Three Center Schwarz direct Screened ints\n";
        }
        world.gop.fence();
        auto t0 = tcc_time::now();

        auto eri_pool
              = tcc::integrals::make_pool(tcc::integrals::make_2body(basis));

        auto eri3 = mpqc_ints::sparse_integrals<integrals::init_schwarz_screen>(
              world, eri_pool, tcc::utility::make_array(basis, basis, basis),
              ta_pass_through);

        auto eri3_norm = eri3("i,j,k").norm(world).get();
        world.gop.fence();
        auto t1 = tcc_time::now();

        if (world.rank() == 0) {
            std::cout << "\tnorm of ints was " << eri3_norm << std::endl;
            std::cout << "\tin " << tcc_time::duration_in_s(t0, t1)
                      << " seconds" << std::endl;
        }
    }

    { // Two electron Four center
        if (world.rank() == 0) {
            std::cout << "\nTwo E four Center Schwarz Screened direct ints\n";
        }
        world.gop.fence();
        auto t0 = tcc_time::now();

        auto eri_pool
              = tcc::integrals::make_pool(tcc::integrals::make_2body(basis));

        auto eri4 = mpqc_ints::
              direct_sparse_integrals<integrals::init_schwarz_screen>(
                    world, eri_pool,
                    tcc::utility::make_array(basis, basis, basis, basis),
                    ta_pass_through);

        auto eri4_norm = eri4("i,j,k,l").norm(world).get();
        world.gop.fence();
        auto t1 = tcc_time::now();

        if (world.rank() == 0) {
            std::cout << "\tnorm of ints was " << eri4_norm << std::endl;
            std::cout << "\tin " << tcc_time::duration_in_s(t0, t1)
                      << " seconds" << std::endl;
        }
    }

    { // Two electron Four center
        if (world.rank() == 0) {
            std::cout << "Two E four Center QQR Screened direct ints\n";
        }
        world.gop.fence();
        auto t0 = tcc_time::now();

        auto eri_pool
              = tcc::integrals::make_pool(tcc::integrals::make_2body(basis));

        // Try looser thresh
        integrals::QQR::well_sep_threshold(0.1);
        auto eri4 = mpqc_ints::
              direct_sparse_integrals<integrals::init_qqr_screen>(
                    world, eri_pool,
                    tcc::utility::make_array(basis, basis, basis, basis),
                    ta_pass_through);

        auto eri4_norm = eri4("i,j,k,l").norm(world).get();
        world.gop.fence();
        auto t1 = tcc_time::now();

        if (world.rank() == 0) {
            std::cout << "\tnorm of ints was " << eri4_norm << std::endl;
            std::cout << "\tin " << tcc_time::duration_in_s(t0, t1)
                      << " seconds" << std::endl;
        }
    }

    { // Two electron Four center
        if (world.rank() == 0) {
            std::cout << "Two E four Center direct unscreened ints\n";
        }
        world.gop.fence();
        auto t0 = tcc_time::now();

        auto eri_pool
              = tcc::integrals::make_pool(tcc::integrals::make_2body(basis));

        auto eri4 = mpqc_ints::direct_sparse_integrals(
              world, eri_pool,
              tcc::utility::make_array(basis, basis, basis, basis),
              ta_pass_through);

        auto eri4_norm = eri4("i,j,k,l").norm(world).get();
        world.gop.fence();
        auto t1 = tcc_time::now();

        if (world.rank() == 0) {
            std::cout << "\tnorm of ints was " << eri4_norm << std::endl;
            std::cout << "\tin " << tcc_time::duration_in_s(t0, t1)
                      << " seconds" << std::endl;
        }
    }

    libint2::cleanup();
    madness::finalize();
    return 0;
}
