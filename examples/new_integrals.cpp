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

#include "../integrals/integral_engine_pool.h"
#include "../integrals/make_engine.h"

#include "../integrals/task_integrals.h"
#include "../integrals/screened_task_integrals.h"
#include "../integrals/direct_task_integrals.h"

#include "../tensor/tcc_tile.h"
#include "../tensor/decomposed_tensor.h"
#include "../tensor/decomposed_tensor_algebra.h"

#include "../utility/time.h"
#include "../utility/array_storage.h"

#include <memory>

using namespace mpqc;

int main(int argc, char *argv[]) {
    auto &world = madness::initialize(argc, argv);
    std::string mol_file = "";
    std::string basis_name = "";
    int nclusters = 0;
    double threshold = 1e-11;
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
    auto mol2 = molecule::attach_hydrogens_and_kmeans(mol.clusterables(),
                                                      nclusters);

    basis::BasisSet bs(basis_name);
    basis::Basis basis(bs.get_cluster_shells(mol2));

    libint2::init();

    auto eri_pool
          = tcc::integrals::make_pool(tcc::integrals::make_2body(basis));

    auto ta_pass_through =
          [](TA::TensorD &&ten) { return TA::TensorD(std::move(ten)); };

    { // Over lap ints
        try {
            if (world.rank() == 0) {
                std::cout << "Overlap ints\n";
            }
            world.gop.fence();
            auto t0 = tcc_time::now();
            auto overlap_pool
                  = tints::make_pool(tints::make_1body("overlap", basis));
            auto S = mpqc_ints::TaskInts(world, overlap_pool,
                                         tcc::utility::make_array(basis, basis),
                                         ta_pass_through);
            auto S_norm = S("i,j").norm(world).get();
            world.gop.fence();
            auto t1 = tcc_time::now();

            if (world.rank() == 0) {
                std::cout << "\tnorm of ints was " << S_norm << std::endl;
                std::cout << "\tin " << tcc_time::duration_in_s(t0, t1)
                          << " seconds" << std::endl;
            }
        } catch (TA::Exception &e) {
            std::cout << "Caught exception: " << e.what() << std::endl;
            std::cout << "\nInput Molecule has " << mol.nclusters()
                      << " atoms\n";
            std::cout << "\nClustered Molecule has " << mol2.nclusters()
                      << " clusters\n";
            for (auto const &c : mol2) {
                std::cout << c << std::endl;
                for (auto const &a : c.atoms()) {
                    std::cout << "\t" << a << std::endl;
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
            std::terminate();
        }
    }

    { // Two electron two center
        if (world.rank() == 0) {
            std::cout << "Two E two Center ints\n";
        }
        world.gop.fence();
        auto t0 = tcc_time::now();

        auto eri_pool
              = tcc::integrals::make_pool(tcc::integrals::make_2body(basis));
        auto eri2 = mpqc_ints::TaskInts(world, eri_pool,
                                        tcc::utility::make_array(basis, basis),
                                        ta_pass_through);

        auto eri2_norm = eri2("i,j").norm(world).get();
        world.gop.fence();
        auto t1 = tcc_time::now();

        if (world.rank() == 0) {
            std::cout << "\tnorm of ints was " << eri2_norm << std::endl;
            std::cout << "\tin " << tcc_time::duration_in_s(t0, t1)
                      << " seconds" << std::endl;
        }
    }

    { // Two electron three center
        if (world.rank() == 0) {
            std::cout << "\nTwo E three Center ints\n";
        }
        world.gop.fence();
        auto t0 = tcc_time::now();

        auto eri_pool = tcc::integrals::make_pool(
              tcc::integrals::make_2body(basis, basis));
        auto eri3 = mpqc_ints::TaskInts(
              world, eri_pool, tcc::utility::make_array(basis, basis, basis),
              ta_pass_through);

        auto eri3_norm = eri3("x,i,j").norm(world).get();
        world.gop.fence();
        auto t1 = tcc_time::now();

        if (world.rank() == 0) {
            std::cout << "\tnorm of ints was " << eri3_norm << std::endl;
            std::cout << "\tin " << tcc_time::duration_in_s(t0, t1)
                      << " seconds" << std::endl;
        }
    }
    {
        if (world.rank() == 0) {
            std::cout << "Two E three Center ints screened\n";
        }
        world.gop.fence();
        auto t0 = tcc_time::now();

        auto eri_pool = tcc::integrals::make_pool(
              tcc::integrals::make_2body(basis, basis));
        auto eri3 = mpqc_ints::ScreenedTaskInts(
              world, eri_pool, tcc::utility::make_array(basis, basis, basis),
              ta_pass_through);

        auto eri3_norm = eri3("x,i,j").norm(world).get();
        world.gop.fence();
        auto t1 = tcc_time::now();

        if (world.rank() == 0) {
            std::cout << "\tnorm of ints was " << eri3_norm << std::endl;
            std::cout << "\tin " << tcc_time::duration_in_s(t0, t1)
                      << " seconds" << std::endl;
        }
    }
    { // Two electron three center
        if (world.rank() == 0) {
            std::cout << "\nTwo E three Center ints Direct\n";
        }
        world.gop.fence();
        auto t0 = tcc_time::now();

        auto eri_pool = tcc::integrals::make_pool(
              tcc::integrals::make_2body(basis, basis));
        auto eri3 = mpqc_ints::direct_task_ints_unscreeened(
              world, eri_pool, tcc::utility::make_array(basis, basis, basis),
              ta_pass_through);

        auto eri3_norm = eri3("x,i,j").norm(world).get();
        world.gop.fence();
        auto t1 = tcc_time::now();

        if (world.rank() == 0) {
            std::cout << "\tnorm of ints was " << eri3_norm << std::endl;
            std::cout << "\tin " << tcc_time::duration_in_s(t0, t1)
                      << " seconds" << std::endl;
        }
    }


    auto clr_tensor = [](TA::TensorD &&t) {
        auto range = t.range();

        auto dense = tcc::tensor::DecomposedTensor<double>(1e-7, std::move(t));
        auto test = tcc::tensor::algebra::two_way_decomposition(dense);

        if (!test.empty()) {
            dense = std::move(test);
        }

        return tcc::tensor::Tile<decltype(dense)>(std::move(range),
                                                  std::move(dense));
    };

    if (world.rank() == 0) {
        std::cout << "\n\nTesting CLR times\n";
    }
    { // Two electron three center CLR
        if (world.rank() == 0) {
            std::cout << "Two E three Center ints CLR no screen\n";
        }
        world.gop.fence();
        auto t0 = tcc_time::now();

        auto eri_pool = tcc::integrals::make_pool(
              tcc::integrals::make_2body(basis, basis));
        auto eri3 = mpqc_ints::TaskInts(
              world, eri_pool, tcc::utility::make_array(basis, basis, basis),
              clr_tensor);

        auto eri3_norm = eri3("x,i,j").norm(world).get();
        world.gop.fence();
        auto t1 = tcc_time::now();

        tcc::utility::print_size_info(eri3, "CLR no screening");
        if (world.rank() == 0) {
            std::cout << "\n\tnorm of ints was " << eri3_norm << std::endl;
            std::cout << "\tin " << tcc_time::duration_in_s(t0, t1)
                      << " seconds" << std::endl;
        }
    }
    {
        if (world.rank() == 0) {
            std::cout << "\nTwo E three Center ints CLR screened\n";
        }
        world.gop.fence();
        auto t0 = tcc_time::now();

        auto eri_pool = tcc::integrals::make_pool(
              tcc::integrals::make_2body(basis, basis));
        auto eri3 = mpqc_ints::ScreenedTaskInts(
              world, eri_pool, tcc::utility::make_array(basis, basis, basis),
              clr_tensor);

        auto eri3_norm = eri3("x,i,j").norm(world).get();
        world.gop.fence();
        auto t1 = tcc_time::now();

        tcc::utility::print_size_info(eri3, "CLR with screening");
        if (world.rank() == 0) {
            std::cout << "\n\tnorm of ints was " << eri3_norm << std::endl;
            std::cout << "\tin " << tcc_time::duration_in_s(t0, t1)
                      << " seconds" << std::endl;
        }
    }

    libint2::cleanup();
    madness::finalize();
    return 0;
}
