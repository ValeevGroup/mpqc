#include <memory>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <chrono>

#include "../include/tbb.h"
#include "../include/libint.h"
#include "../include/tiledarray.h"
#include "../include/btas.h"

#include "../utility/make_array.h"
#include "../utility/parallel_print.h"
#include "../utility/array_storage.h"
#include "../utility/time.h"
#include "../utility/ta_helpers.h"

#include "../tensor/conversions/tile_pimpl_to_ta_tensor.h"

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
#include "../integrals/btas_to_low_rank_tensor.h"
#include "../integrals/make_engine.h"
#include "../integrals/ta_tensor_to_low_rank_tensor.h"
#include "../integrals/integral_engine_pool.h"
#include "../integrals/sparse_task_integrals.h"
#include "../integrals/dense_task_integrals.h"

#include "../purification/sqrt_inv.h"
#include "../purification/purification_devel.h"


using namespace tcc;
namespace ints = integrals;

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
    volatile int debug = 0;
    if (argc >= 3) {
        mol_file = argv[1];
        basis_name = argv[2];
        nclusters = std::stoi(argv[3]);
    } else {
        std::cout << "input is $./program mol_file basis_file df_basis_file "
                     "nclusters \n";
        return 0;
    }
    if (argc == 5) {
        threshold = std::stod(argv[4]);
    }
    if (argc == 6) {
        debug = std::stoi(argv[5]);
    }
    TiledArray::SparseShape<float>::threshold(threshold);
    utility::print_par(world, "Sparse threshold is ",
                       TiledArray::SparseShape<float>::threshold(), "\n");

    auto mol = molecule::read_xyz(mol_file);

    utility::print_par(world, "Computing ", mol.nelements(), " elements with ",
                       nclusters, " clusters.\n");

    auto clusters = molecule::attach_hydrogens_kmeans(mol, nclusters);

    basis::BasisSet bs{basis_name};
    basis::Basis basis{bs.create_basis(clusters)};


    libint2::init();

    // Begin Two electron integrals section.
    auto eri_pool
        = ints::make_pool(ints::make_2body(basis, basis, basis, basis));

    class btas_to_four_k {
      private:
        double low_rank_threshold;

      public:
        using TileType = tensor::TilePimpl<double>;

        btas_to_four_k(double l) : low_rank_threshold(l) {}

        tensor::TilePimpl<double>
        operator()(tensor::ShallowTensor<4> const bt) {
            auto const &extent = bt.tensor().extent();
            auto ij = extent[0] * extent[1];
            auto prange = btas::permute(bt.tensor().range(), {0, 3, 2, 1});
            remove_ref_const_t<decltype(bt.tensor())> view
                = btas::make_view(prange, bt.tensor().storage());
            RowMatrixXd Tile(ij, ij);
            for (auto i = 0; i < Tile.size(); ++i) {
                auto elem = view[i];
                *(Tile.data() + i) = elem;
            }

            if (Tile.lpNorm<2>()
                >= TiledArray::SparseShape<float>::threshold()) {
                RowMatrixXd L, R;
                bool is_full
                    = algebra::Decompose_Matrix(Tile, L, R, low_rank_threshold);
                if (!is_full) {
                    tensor::TileVariant<double> tile_variant{
                        tensor::LowRankTile<double>{std::move(L),
                                                    std::move(R)}};

                    return tensor::TilePimpl<double>{bt.range(),
                                                     std::move(tile_variant),
                                                     low_rank_threshold};
                }

                tensor::TileVariant<double> tile_variant{
                    tensor::FullRankTile<double>{std::move(Tile)}};

                return tensor::TilePimpl<double>{bt.range(),
                                                 std::move(tile_variant)};
            } else { // For now if tile is super sparse return a 1x1 zero tile.
                Tile = RowMatrixXd::Zero(1, 1);
                tensor::TileVariant<double> tile_variant{
                    tensor::FullRankTile<double>{std::move(Tile)}};

                return tensor::TilePimpl<double>{bt.range(),
                                                 std::move(tile_variant)};
            }
        }
    };

    // Compute center integrals
    utility::print_par(world, "\n");

    auto EriK = ints::BlockSparseIntegrals(
        world, eri_pool, utility::make_array(basis, basis, basis, basis),
        btas_to_four_k{low_rank_threshold});

    print_size_info(EriK, "Low Rank Exchange Tensor");


    world.gop.fence();
    libint2::cleanup();
    madness::finalize();
    return 0;
}
