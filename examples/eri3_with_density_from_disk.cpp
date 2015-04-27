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

#include "../utility/parallel_break_point.h"
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
    utility::parallal_break_point(world, debug);

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
        for (auto i = 0ul; i < tile_norms.size(); ++i) {
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
        for (auto i = 0ul; i < pmap->size(); ++i) {
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

    world.gop.fence();

    libint2::init();
    auto eri_pool = ints::make_pool(ints::make_2body(basis, df_basis));

    utility::print_par(world, "\n");
    auto basis_array = utility::make_array(df_basis, basis, basis);
    auto Xab
          = BlockSparseIntegrals(world, eri_pool, basis_array,
                                 integrals::compute_functors::BtasToTaTensor{});

    auto func = [=](TA::Tensor<double> const &t) {

        tcc::tensor::DecomposedTensor<double> temp(low_rank_threshold,
                                                   t.clone());
        auto test_me = tensor::algebra::two_way_decomposition(temp);

        if (test_me.empty()) {
            auto const &extent = t.range().size();
            TA::Range new_range{extent[0], extent[1], extent[2]};
            test_me = tcc::tensor::DecomposedTensor<double>(
                  low_rank_threshold, TA::Tensor<double>{new_range, t.data()});
        }

        return tcc::tensor::Tile<decltype(test_me)>(t.range(),
                                                    std::move(test_me));
    };

    auto func2 = [=](TA::Tensor<double> const &t) {
        auto const &extent = t.range().size();
        TA::Range new_range{extent[0], extent[1]};
        TA::Tensor<double> new_tensor(new_range, t.data());

        return tcc::tensor::Tile<tcc::tensor::DecomposedTensor<double>>{
              t.range(), tcc::tensor::DecomposedTensor<double>(
                               low_rank_threshold, std::move(new_tensor))};
    };
    auto Xab_lr = TA::to_new_tile_type(Xab, func);
    auto D_test = TA::to_new_tile_type(D_TA, func2);
    world.gop.fence();
    utility::print_size_info(Xab_lr, "Xab_lr");
    utility::print_array_difference(Xab_lr, Xab, "Xab low rank",
                                    "Xab full rank");

    auto t_ta0 = std::chrono::high_resolution_clock::now();
    decltype(Xab) Xak;
    Xak("X, k, a") = Xab("X,a,b") * D_TA("b,k");
    auto t_ta1 = std::chrono::high_resolution_clock::now();
    auto time_ta = std::chrono::duration_cast<std::chrono::duration<double>>(
                         t_ta1 - t_ta0).count();
    utility::print_par(world, "\nTime for TA contraction = ", time_ta, "\n");

    auto t_me0 = std::chrono::high_resolution_clock::now();
    decltype(Xab_lr) Xak_lr;
    Xak_lr("X,k,a") = Xab_lr("X,a,b") * D_test("b,k");
    auto t_me1 = std::chrono::high_resolution_clock::now();
    auto time_me = std::chrono::duration_cast<std::chrono::duration<double>>(
                         t_me1 - t_me0).count();

    utility::print_par(world, "\nTime for My contraction = ", time_me, "\n");
    utility::print_par(world, "Speed up over fully dense = ", time_ta / time_me,
                       "\n\n");
    utility::print_size_info(Xak_lr, "Xak_lr no recompression");
    utility::print_array_difference(Xak_lr, Xak, "Xak low rank",
                                    "Xak full rank");

    // Look at going all the way to K
    {
        auto dfbasis_array = utility::make_array(df_basis, df_basis);
        auto eri2 = BlockSparseIntegrals(
              world, eri_pool, dfbasis_array,
              integrals::compute_functors::BtasToTaTensor{});

        auto inv_timer = tcc_time::make_timer(
              [&]() { return pure::inverse_sqrt(eri2); });
        auto eri2_inv = inv_timer.apply();
        utility::print_par(world, "\nEri2 inverse computation time = ",
                           inv_timer.time(), "\n");
        eri2_inv("i,j") = eri2_inv("i,k") * eri2_inv("k,j");
        eri2_inv.truncate();

        auto V_inv_to_lr = [=](TA::Tensor<double> const &t) {
            tcc::tensor::DecomposedTensor<double> temp(low_rank_threshold,
                                                       t.clone());
            auto test_me = tensor::algebra::two_way_decomposition(temp);
            /* if(true){ */
            if (test_me.empty()) { // was not low rank.
                auto const &extent = t.range().size();
                TA::Range new_range{extent[0], extent[1]};
                test_me = tcc::tensor::DecomposedTensor<double>(
                      low_rank_threshold,
                      TA::Tensor<double>{new_range, t.data()});
            }

            return tcc::tensor::Tile<decltype(test_me)>(t.range(),
                                                        std::move(test_me));
        };

        auto V_inv_lr = TA::to_new_tile_type(eri2_inv, V_inv_to_lr);

        utility::print_par(world, "\n");
        utility::print_array_difference(V_inv_lr, eri2_inv, "V^{-1} lr",
                                        "V^{-1}");
        utility::print_size_info(V_inv_lr, "V^{-1}");

        auto lr_0 = tcc_time::now();
        decltype(Xab_lr) W_lr;
        W_lr("X,a,b") = V_inv_lr("X,P") * Xab_lr("P,a,b");
        auto lr_1 = tcc_time::now();
        auto lr_time = tcc_time::duration_in_s(lr_0, lr_1);

        auto f_0 = tcc_time::now();
        decltype(Xab) W;
        W("X, a, b") = eri2_inv("X,P") * Xab("P, a, b");
        auto f_1 = tcc_time::now();
        auto f_time = tcc_time::duration_in_s(f_0, f_1);

        utility::print_par(world, "\nLow rank gemm time = ", lr_time,
                           " full rank time = ", f_time, " speed up = ",
                           f_time / lr_time, "\n");
        utility::print_array_difference(W_lr, W, "W_lr", "W full");
        utility::print_size_info(W_lr, "W");

        auto j_0 = tcc_time::now();
        TA::Array<double, 2, TA::Tensor<double>, TA::SparsePolicy> J;
        J("i,j") = W("X,i,j") * (Xak("X,a,b") * D_TA("a,b"));
        auto j_1 = tcc_time::now();
        auto j_time = tcc_time::duration_in_s(j_0, j_1);

        auto jlr_0 = tcc_time::now();
        decltype(D_test) J_lr;
        J_lr("i,j") = W_lr("X,i,j") * (Xab_lr("X,a,b") * D_test("a,b"));
        auto jlr_1 = tcc_time::now();
        auto jlr_time = tcc_time::duration_in_s(jlr_0, jlr_1);
        utility::print_par(world, "\nLow rank J time = ", jlr_time,
                           " full rank time = ", j_time, " speed up = ",
                           j_time / jlr_time, "\n");
        utility::print_array_difference(J_lr, J, "J_lr", "J full");
        utility::print_size_info(J_lr, "J");

        auto k_0 = tcc_time::now();
        TA::Array<double, 2, TA::Tensor<double>, TA::SparsePolicy> K;
        K("i,j") = Xak("X,k,i") * W("X,k,j");
        auto k_1 = tcc_time::now();
        auto k_time = tcc_time::duration_in_s(k_0, k_1);

        auto klr_0 = tcc_time::now();
        TA::Array<double, 2, tensor::Tile<tensor::DecomposedTensor<double>>,
                  TA::SparsePolicy> K_lr;
        K_lr("i,j") = Xak_lr("X,k,i") * W_lr("X,k,j");
        auto klr_1 = tcc_time::now();
        auto klr_time = tcc_time::duration_in_s(klr_0, klr_1);

        utility::print_par(world, "\nLow rank K time = ", klr_time,
                           " full rank time = ", k_time, " speed up = ",
                           k_time / klr_time, "\n");
        utility::print_array_difference(K_lr, K, "K_lr", "K full");
        utility::print_size_info(K_lr, "K");

        auto k_iter_time = k_time + time_ta;
        auto klr_iter_time = klr_time + time_me;
        utility::print_par(world, "\nLow rank total K time = ", klr_iter_time,
                           " full rank time = ", k_iter_time, " speed up = ",
                           k_iter_time / klr_iter_time, "\n");
    }

    madness::finalize();
    return 0;
}
