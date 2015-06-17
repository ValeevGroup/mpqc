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

#include "../density/purification_devel.h"
#include "../density/sqrt_inv.h"

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
    bool use_columb_metric = true;
    if (argc >= 4) {
        mol_file = argv[1];
        basis_name = argv[2];
        nclusters = std::stoi(argv[3]);
    } else {
        std::cout << "input is $./program mol_file basis_file nclusters "
                     "use_coulomb_metric (Optioinal set to 0 to use overlap "
                     "matrix) sparse_threshold (Optional, default is 1e-11)\n";

        return 0;
    }
    if (argc >= 5) {
        utility::print_par(world, "User Selected Metric\n");
        use_columb_metric = std::stoi(argv[4]);
    }
    if (argc == 6) {
        threshold = std::stod(argv[5]);
    }
    TiledArray::SparseShape<float>::threshold(threshold);
    utility::print_par(world, "Sparse threshold is ",
                       TiledArray::SparseShape<float>::threshold(), "\n");

    auto mol = molecule::read_xyz(mol_file);
    auto repulsion_energy = mol.nuclear_repulsion();

    utility::print_par(world, "Computing ", mol.nelements(), " elements with ",
                       nclusters, " clusters.\n");

    auto clusters = molecule::attach_hydrogens_kmeans(mol, nclusters);

    std::streambuf *cout_sbuf = std::cout.rdbuf(); // Silence libint printing.
    std::ofstream fout("/dev/null");
    std::cout.rdbuf(fout.rdbuf());
    basis::BasisSet bs{basis_name};
    basis::Basis basis{bs.create_basis(clusters)};
    std::cout.rdbuf(cout_sbuf);

    libint2::init();
    utility::print_par(world, "Starting integrals ");

    TAArray<2, TA::SparsePolicy> eri2;

    if (use_columb_metric) {
        utility::print_par(world, "Using Coulomb Metric\n");
        auto eri_pool = integrals::make_pool(integrals::make_2body(basis));

        eri2 = integrals::BlockSparseIntegrals(
              world, eri_pool, utility::make_array(basis, basis),
              integrals::compute_functors::BtasToTaTensor{});
        world.gop.fence();
    } else {
        utility::print_par(world, "Using Overlap Metric\n");
        auto eri_pool
              = integrals::make_pool(integrals::make_1body("overlap", basis));

        eri2 = integrals::BlockSparseIntegrals(
              world, eri_pool, utility::make_array(basis, basis),
              integrals::compute_functors::BtasToTaTensor{});
        world.gop.fence();
    }

    auto volume = eri2.elements().volume();
    { // Compute low rank information for Eri2
        eri2.truncate();
        utility::print_par(world, "\nM Tensor volume = ", volume, "\n");
        utility::print_par(world, "M row and col dim = ", std::sqrt(volume), "\n");
        print_size_info(eri2, "M");
        if (world.rank() == 0) {
            std::cout << "\tM sparsity percent = "
                      << eri2.get_shape().sparsity() << "\n";
        }
        auto shape = eri2.get_shape().data();
        if (world.rank() == 0) {
            std::ofstream tiles_file("tiles.txt");
            tiles_file << "tile: (rows, cols) volume scaled_norm\n";
            auto volume = eri2.trange().tiles().volume();
            for (auto i = 0; i < volume; ++i) {
                    auto range = eri2.trange().make_tile_range(i);
                    tiles_file << i << ": ("
                              << range.size()[0] << " "
                              << range.size()[1] << ") "
                              << range.volume() << " "
                              << shape[i]
                              << "\n";
            }
        }
    }

    if (world.rank() == 0) {
        std::cout << "\nStarting inverse sqrt." << std::endl;
    }
    {
        auto inv_sqrt_timer = tcc_time::make_timer(
              [&]() { return pure::inverse_sqrt(eri2); });

        auto sqrt_inv = inv_sqrt_timer.apply();

        utility::print_par(world, "M^{-1/2} took in total ",
                           inv_sqrt_timer.time(), " s\n");
        print_size_info(sqrt_inv, "M^{-1/2}");
        if (world.rank() == 0) {
            std::cout << "\tM^{-1/2} sparsity percent = "
                      << sqrt_inv.get_shape().sparsity() << "\n";
        }

        // Check accuracy
        decltype(sqrt_inv) inv;
        inv("i,j") = sqrt_inv("i,k") * sqrt_inv("k,j");
        inv.truncate();
        // Reuse inv to compute approximate identity.
        inv("i,j") = inv("i,k") * eri2("k,j");
        inv.truncate();

        auto ident = pure::create_diagonal_matrix(inv, 1.0);

        auto f_norm_diff = utility::array_fnorm_diff(inv, ident);
        utility::print_par(world, "\nIdenity-(M*M^{-1}) F norm difference = ",
                           f_norm_diff, "\n");
        utility::print_par(world,
                           "Idenity-(M*M^{-1}) F norm difference / volume = ",
                           f_norm_diff / volume, "\n");
    }

    world.gop.fence();
    libint2::cleanup();
    madness::finalize();
    return 0;
}
