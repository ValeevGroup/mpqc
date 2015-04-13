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

int main(int argc, char **argv) {
    auto &world = madness::initialize(argc, argv);
    std::string mol_file = "";
    std::string basis_name = "";
    std::string df_basis_name = "";
    std::string density_file = "";
    int bs_nclusters = 0;
    int dfbs_nclusters = 0;
    if (argc >= 6) {
        mol_file = argv[1];
        basis_name = argv[2];
        df_basis_name = argv[3];
        bs_nclusters = std::stoi(argv[4]);
        dfbs_nclusters = std::stoi(argv[5]);
    } else if (argc >= 9 || argc < 6) {
        std::cout << "input is $./program mol_file "
                     "basis_file df_basis_file "
                     "bs_clusters dfbs_clusters low_rank_threshhold(optional) "
                     "debug(optional)\n";
        return 0;
    }

    double threshold = (argc == 7) ? std::stod(argv[6]) : 1e-11;
    auto low_rank_threshold = (argc == 8) ? std::stod(argv[7]) : 1e-7;
    volatile auto debug = (argc == 9) ? std::stoi(argv[8]) : 0;

    if (world.rank() == 0) {
        std::cout << "Mol file is " << mol_file << std::endl;
        std::cout << "basis is " << basis_name << std::endl;
        std::cout << "df basis is " << df_basis_name << std::endl;
        std::cout << "Using " << bs_nclusters << " bs clusters" << std::endl;
        std::cout << "Using " << dfbs_nclusters << " dfbs clusters"
                  << std::endl;
        std::cout << "low rank threshhold is " << low_rank_threshold
                  << std::endl;
        if (debug) {
            std::cout << "Debuging requested!" << std::endl;
        }
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

    libint2::init();
    auto eri_pool = ints::make_pool(ints::make_2body(basis, df_basis));

    utility::print_par(world, "\n");
    auto basis_array = utility::make_array(df_basis, basis, basis);
    auto Xab
          = BlockSparseIntegrals(world, eri_pool, basis_array,
                                 integrals::compute_functors::BtasToTaTensor{});

    {
        auto dfbasis_array = utility::make_array(df_basis, df_basis);
        auto eri2 = BlockSparseIntegrals(
              world, eri_pool, dfbasis_array,
              integrals::compute_functors::BtasToTaTensor{});

        auto eri2_lr = TA::to_new_tile_type(
              eri2, integrals::compute_functors::TaToLowRankTensor<2>(
                          low_rank_threshold));

        utility::print_par(world, "\n");
        utility::print_size_info(eri2_lr, "Eri2");

        auto inv_timer = tcc_time::make_timer(
              [&]() { return pure::inverse_sqrt(eri2); });
        auto eri2_inv = inv_timer.apply();
        utility::print_par(world, "\nEri2 inverse computation time = ",
                           inv_timer.time(), "\n");
        eri2_inv("i,j") = eri2_inv("i,k") * eri2_inv("k,j");
        eri2_inv.truncate();

        auto eri2_inv_lr = TA::to_new_tile_type(
              eri2_inv, integrals::compute_functors::TaToLowRankTensor<2>(
                              low_rank_threshold));

        utility::print_par(world, "\n");
        utility::print_size_info(eri2_inv_lr, "Eri2 inverse");

        Xab("X, a, b") = eri2_inv("X,P") * Xab("P, a, b");
        Xab.truncate();

        auto Xab_lr = TA::to_new_tile_type(
              Xab, integrals::compute_functors::TaToLowRankTensor<3>(
                         low_rank_threshold));

        utility::print_par(world, "\n");
        utility::print_size_info(Xab_lr, "\\tilde{X}ab");
    }

    madness::finalize();
    return 0;
}
