#include <memory>
#include <fstream>
#include <algorithm>

#include "../common/namespaces.h"
#include "../common/typedefs.h"

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

#include "../basis/atom_basisset.h"
#include "../basis/basis_set.h"
#include "../basis/cluster_shells.h"
#include "../basis/basis.h"

#include "../integrals/btas_to_ta_tensor.h"
#include "../integrals/integral_engine_pool.h"
#include "../integrals/sparse_task_integrals.h"
#include "../integrals/make_engine.h"

/* #include "../scf/purification/purification_devel.h" */
#include "../ta_routines/sqrt_inv.h"
#include "../ta_routines/inverse.h"
#include <madness/world/array_addons.h>

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
    double threshold = 1e-14;
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
        utility::print_par(world, "M row and col dim = ", std::sqrt(volume),
                           "\n");
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
                tiles_file << i << ": (" << range.extent()[0] << " "
                           << range.extent()[1] << ") " << range.volume() << " "
                           << shape[i] << "\n";
            }
        }
    }
    {
        using Matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                                     Eigen::RowMajor>;
        int64_t dim = std::sqrt(eri2.trange().elements().volume());
        Matrix Print_mat = Matrix::Zero(dim, dim);
        for (auto it : eri2) {
            auto tile = it.get();
            auto const &range = tile.range();

            auto const &extent = range.extent();
            auto const &lb = range.lobound();
            auto const &ub = range.upbound();

            Matrix block(extent[0], extent[1]);

            for (auto i = 0; i < block.size(); ++i) {
                *(block.data() + i) = *(tile.data() + i);
            }

            Print_mat.block(lb[0], lb[1], extent[0], extent[1]) = block;
        }

        std::ofstream mat_file("dense_matrix_file.dat", std::ios::binary);
        auto rows = Print_mat.rows();
        auto cols = Print_mat.cols();
        mat_file.write((char*) Print_mat.data(), rows * cols * sizeof(typename Matrix::Scalar));
        mat_file.close();
    }

    if (world.rank() == 0) {
        std::cout << "\nStarting inverse sqrt." << std::endl;
    }
    if (false) {
        auto inv_sqrt_timer = tcc_time::make_timer(
              [&]() { return tcc::pure::inverse_sqrt(eri2); });

        auto sqrt_inv = inv_sqrt_timer.apply();

        utility::print_par(world, "M^{-1/2} took in total ",
                           inv_sqrt_timer.time(), " s\n\n");
    }

    if (world.rank() == 0) {
        std::cout << "\nStarting LR inverse sqrt." << std::endl;
    }
    {
        auto to_decomp_with_decompose = [=](TA::Tensor<double> const &t) {
            auto range = t.range();

            auto const extent = range.extent();
            const auto i = extent[0];
            const auto j = extent[1];
            auto local_range = TA::Range{i, j};

            auto tensor = TA::Tensor<double>(local_range, t.data());
            auto dense
                  = tensor::DecomposedTensor<double>(1e-6, std::move(tensor));

            auto test = tensor::algebra::two_way_decomposition(dense);
            if (!test.empty()) {
                dense = std::move(test);
            }

            return tensor::Tile<tensor::DecomposedTensor<double>>(
                  range, std::move(dense));
        };

        auto eri2_lr = TA::to_new_tile_type(eri2, to_decomp_with_decompose);
        {
            using Matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                                         Eigen::RowMajor>;
            int64_t dim = std::sqrt(eri2_lr.trange().elements().volume());
            Matrix Print_mat = Matrix::Zero(dim, dim);
            for (auto it : eri2_lr) {
                auto tile_wrapper = it.get();
                auto const &range = tile_wrapper.range();
                auto tile = tile_wrapper.tile();

                auto const &extent = range.extent();
                auto const &lb = range.lobound();
                auto const &ub = range.upbound();

                Matrix block = Matrix::Zero(extent[0], extent[1]);

                if (tile.ndecomp() == 1) {
                    for (auto i = 0; i < block.size(); ++i) {
                        *(block.data() + i) = *(tile.tensors()[0].data() + i);
                    }
                } else {
                    auto const &t1 = tile.tensors()[0];
                    auto e1 = t1.range().extent();
                    Eigen::Map<const Matrix> map(t1.data(), e1[0], e1[1]);
                    block.block(0, 0, e1[0], e1[1]) = map;

                    auto const &t2 = tile.tensors()[1];
                    auto e2 = t2.range().extent();
                    Eigen::Map<const Matrix> map2(t2.data(), e2[0], e2[1]);
                    block.block(0, 0, e2[0], e2[1]) = map2;
                }

                Print_mat.block(lb[0], lb[1], extent[0], extent[1]) = block;
            }

            std::ofstream mat_file("lr_matrix_file.dat", std::ios::binary);
            auto rows = Print_mat.rows();
            auto cols = Print_mat.cols();
            mat_file.write((char*) Print_mat.data(), rows * cols * sizeof(typename Matrix::Scalar));
            mat_file.close();
        }
        //   auto inv_sqrt_timer = tcc_time::make_timer(
        //         [&]() { return tcc::pure::inverse_sqrt(eri2_lr); });

        //   auto sqrt_inv = inv_sqrt_timer.apply();

        //   utility::print_par(world, "M^{-1/2} took in total ",
        //                      inv_sqrt_timer.time(), " s\n");
    }

    world.gop.fence();
    libint2::cleanup();
    madness::finalize();
    return 0;
}
