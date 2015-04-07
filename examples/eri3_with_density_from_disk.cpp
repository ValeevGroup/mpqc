#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/array.hpp>

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
#include "../tensor/tcc_tile.h"
#include "../tensor/decomposed_tensor.h"
#include "../tensor/decomposed_tensor_nonintrusive_interface.h"

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

namespace boost {
namespace serialization {

template <typename Archive>
void serialize(Archive &ar, RowMatrixXd &m, const unsigned int version) {
    unsigned int rows = 0;
    unsigned int cols = 0;
    ar &rows;
    ar &cols;

    m.resize(rows, cols);

    ar &boost::serialization::make_array(m.data(), m.size());
}


} // serialization
} // boost

RowMatrixXd read_density_from_file(std::string const &file_name) {
    RowMatrixXd D;
    std::ifstream dfile(file_name.c_str());
    if (dfile.good()) {
        boost::archive::binary_iarchive ia(dfile, std::ios::binary);
        ia >> D;
        dfile.close();
    } else {
        throw;
    }

    return D;
}

int main(int argc, char **argv) {
    auto &world = madness::initialize(argc, argv);
    std::string mol_file = "";
    std::string basis_name = "";
    std::string df_basis_name = "";
    std::string density_file = "";
    int bs_nclusters = 0;
    int dfbs_nclusters = 0;
    if (argc >= 7) {
        mol_file = argv[1];
        density_file = argv[2];
        basis_name = argv[3];
        df_basis_name = argv[4];
        bs_nclusters = std::stoi(argv[5]);
        dfbs_nclusters = std::stoi(argv[6]);
    } else if (argc >= 10 || argc < 7) {
        std::cout << "input is $./program mol_file density_matrix_file "
                     "basis_file df_basis_file "
                     "bs_clusters dfbs_clusters low_rank_threshhold(optional) "
                     "debug(optional)\n";
        return 0;
    }

    double threshold = (argc == 8) ? std::stod(argv[7]) : 1e-11;
    auto low_rank_threshold = (argc == 9) ? std::stod(argv[8]) : 1e-7;
    volatile auto debug = (argc == 10) ? std::stoi(argv[9]) : 0;

    if (world.rank() == 0) {
        std::cout << "Mol file is " << mol_file << std::endl;
        std::cout << "D file is " << density_file << std::endl;
        std::cout << "basis is " << basis_name << std::endl;
        std::cout << "df basis is " << df_basis_name << std::endl;
        std::cout << "Using " << bs_nclusters << " bs clusters" << std::endl;
        std::cout << "Using " << dfbs_nclusters << " dfbs clusters"
                  << std::endl;
        std::cout << "low rank threshhold is " << low_rank_threshold
                  << std::endl;
    }

    TiledArray::SparseShape<float>::threshold(threshold);
    utility::print_par(world, "Sparse threshold is ",
                       TiledArray::SparseShape<float>::threshold(), "\n");

    auto mol = molecule::read_xyz(mol_file);

    auto bs_clusters = molecule::attach_hydrogens_kmeans(mol, bs_nclusters);
    auto dfbs_clusters = molecule::attach_hydrogens_kmeans(mol, dfbs_nclusters);

    basis::BasisSet bs{basis_name};
    basis::BasisSet df_bs{df_basis_name};

    std::streambuf *cout_sbuf = std::cout.rdbuf(); // Silence libint printing.
    std::ofstream fout("/dev/null");
    std::cout.rdbuf(fout.rdbuf());
    basis::Basis basis{bs.create_basis(bs_clusters)};
    basis::Basis df_basis{df_bs.create_basis(dfbs_clusters)};
    std::cout.rdbuf(cout_sbuf);

    auto tr1 = basis.create_flattend_trange1();
    auto tr = TA::TiledRange{tr1, tr1};

    RowMatrixXd D_eig;
    if (world.rank() == 0) {
        D_eig = read_density_from_file(density_file);
    }
    world.gop.fence();

    TA::Tensor<float> tile_norms(tr.tiles(), 0.0);

    if (world.rank() == 0) {
        for (auto i = 0; i < tile_norms.size(); ++i) {
            auto range = tr.make_tile_range(i);
            auto const &size = range.size();
            auto const &start = range.start();

            tile_norms[i] = D_eig.block(start[0], start[1], size[0], size[1])
                                  .lpNorm<2>();
        }
    }

    TA::SparseShape<float> shape(world, tile_norms, tr);

    auto D_TA = TA::Array<double, 2, TA::Tensor<double>, TA::SparsePolicy>(
          world, tr, shape);

    auto const &pmap = D_TA.get_pmap();

    if (world.rank() == 0) {
        for (auto i = 0; i < pmap->size(); ++i) {
            if (!D_TA.is_zero(i)) {
                auto range = tr.make_tile_range(i);
                auto const &start = range.start();
                auto const &size = range.size();
                auto tile = TA::Tensor<double>(range);
                auto tile_map = TA::eigen_map(tile, size[0], size[1]);
                tile_map = D_eig.block(start[0], start[1], size[0], size[1]);
                D_TA.set(i, tile);
            }
        }
    }

    {
        auto D_lr = TA::to_new_tile_type(
              D_TA, integrals::compute_functors::TaToLowRankTensor<2>(
                          low_rank_threshold));
        utility::print_par(world, "\n");
        utility::print_size_info(D_lr, "D");
    }

    libint2::init();
    auto eri_pool = ints::make_pool(ints::make_2body(basis, df_basis));

    utility::print_par(world, "\n");
    auto basis_array = utility::make_array(df_basis, basis, basis);
    auto Xab
          = BlockSparseIntegrals(world, eri_pool, basis_array,
                                 integrals::compute_functors::BtasToTaTensor{});

    auto t_ta0 = std::chrono::high_resolution_clock::now();
    decltype(Xab) Xak;
    Xak("X, a, k") = Xab("X,a,b") * D_TA("b,k");
    auto t_ta1 = std::chrono::high_resolution_clock::now();
    auto time = std::chrono::duration_cast<std::chrono::duration<double>>(
                      t_ta1 - t_ta0).count();
    utility::print_par(world, "Time for TA contraction = ", time, "\n");

    auto func = [=](TA::Tensor<double> const &t) {
        tcc::tensor::DecomposedTensor<double> temp(1e-7, t);
        auto test_me = tensor::algebra::two_way_decomposition(temp);
        test_me = (test_me.empty()) ? temp : test_me;
        return tcc::tensor::Tile<decltype(test_me)>(t.range(),
                                                    std::move(test_me));
    };
    auto func2 = [](TA::Tensor<double> const &t) {
        return tcc::tensor::Tile<tcc::tensor::DecomposedTensor<double>>{
              t.range(), tcc::tensor::DecomposedTensor<double>{1e-7, t}};
    };
    auto Xab_lr = TA::to_new_tile_type(Xab, func);
    auto D_test = TA::to_new_tile_type(D_TA, func2);

    auto t_me0 = std::chrono::high_resolution_clock::now();
    decltype(Xab_lr) Xak_lr;
    Xak_lr("X,a,k") = Xab_lr("X,a,b") * D_test("b,k");
    auto t_me1 = std::chrono::high_resolution_clock::now();
    time = std::chrono::duration_cast<std::chrono::duration<double>>(
                 t_me1 - t_me0).count();

    utility::print_par(world, "Time for My contraction = ", time, "\n");

    /* Xak.truncate(); */
    /* world.gop.fence(); */
    /* auto Xak_lr = TA::to_new_tile_type( */
    /*       Xak, integrals::compute_functors::TaToLowRankTensor<3>( */
    /*                  low_rank_threshold)); */
    /* world.gop.fence(); */
    /* utility::print_size_info(Xak_lr, "Xab * D"); */


    madness::finalize();
    return 0;
}
