#include "../common/namespaces.h"
#include "../common/typedefs.h"

#include "../include/libint.h"
#include "../include/tiledarray.h"
#include "../include/btas.h"

#include "../utility/make_array.h"
#include "../utility/time.h"

#include "../molecule/atom.h"
#include "../molecule/cluster.h"
#include "../molecule/molecule.h"
#include "../molecule/clustering_functions.h"
#include "../molecule/make_clusters.h"

#include "../basis/atom_basisset.h"
#include "../basis/basis_set.h"
#include "../basis/cluster_shells.h"
#include "../basis/basis.h"

#include "../integrals/btas_to_ta_tensor.h"
#include "../integrals/integral_engine_pool.h"
#include "../integrals/sparse_task_integrals.h"
#include "../integrals/make_engine.h"

#include "../tensor/tcc_tile.h"
#include "../tensor/decomposed_tensor.h"
#include "../tensor/decomposed_tensor_algebra.h"
#include "../tensor/decomposed_tensor_unary.h"

#include <memory>

// #include "../ta_routines/sqrt_inv.h"
// #include "../ta_routines/inverse.h"
// #include <madness/world/array_addons.h>

using namespace tcc;

int main(int argc, char *argv[]) {
    auto &world = madness::initialize(argc, argv);
    std::string mol_file = "";
    std::string basis_name = "";
    int nclusters = 0;
    double threshold = 1e-13;
    if (argc >= 4) {
        mol_file = argv[1];
        basis_name = argv[2];
        nclusters = std::stoi(argv[3]);
    } else {
        std::cout << "input is $./program mol_file basis_file nclusters ";
        return 0;
    }
    TiledArray::SparseShape<float>::threshold(threshold);

    auto mol = molecule::read_xyz(mol_file);
    auto clusters = molecule::attach_hydrogens_kmeans(mol, nclusters);

    std::streambuf *cout_sbuf = std::cout.rdbuf(); // Silence libint printing.
    std::ofstream fout("/dev/null");
    std::cout.rdbuf(fout.rdbuf());
    basis::BasisSet bs{basis_name};
    basis::Basis basis{bs.create_basis(clusters)};
    std::cout.rdbuf(cout_sbuf);

    libint2::init();

    auto btas_to_ta = tints::compute_functors::BtasToTaTensor();

    { // Over lap ints
        if (world.rank() == 0) {
            std::cout << "Overlap ints\n";
        }
        world.gop.fence();
        auto t0 = tcc_time::now();
        auto overlap_pool
              = tints::make_pool(tints::make_1body("overlap", basis));
        auto S = tints::BlockSparseIntegrals(world, overlap_pool,
                                             utility::make_array(basis, basis),
                                             btas_to_ta);
        auto S_norm = S("i,j").norm(world).get();
        world.gop.fence();
        auto t1 = tcc_time::now();

        if (world.rank() == 0) {
            std::cout << "\tnorm of ints was " << S_norm << std::endl;
            std::cout << "\tin " << tcc_time::duration_in_s(t0, t1)
                      << " seconds" << std::endl;
        }
    }

    { // Two electron two center
        if (world.rank() == 0) {
            std::cout << "Two E two Center ints\n";
        }
        world.gop.fence();
        auto t0 = tcc_time::now();

        auto eri_pool = integrals::make_pool(integrals::make_2body(basis));
        auto eri2 = integrals::BlockSparseIntegrals(
              world, eri_pool, utility::make_array(basis, basis),
              tints::compute_functors::BtasToTaTensor());

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
            std::cout << "Two E three Center ints\n";
        }
        world.gop.fence();
        auto t0 = tcc_time::now();

        auto eri_pool = integrals::make_pool(integrals::make_2body(basis, basis));
        auto eri3 = integrals::BlockSparseIntegrals(
              world, eri_pool, utility::make_array(basis, basis, basis),
              tints::compute_functors::BtasToTaTensor());

        auto eri3_norm = eri3("x,i,j").norm(world).get();
        world.gop.fence();
        auto t1 = tcc_time::now();

        if (world.rank() == 0) {
            std::cout << "\tnorm of ints was " << eri3_norm << std::endl;
            std::cout << "\tin " << tcc_time::duration_in_s(t0, t1)
                      << " seconds" << std::endl;
        }
    }

    libint2::cleanup();
    madness::finalize();
    return 0;
}
