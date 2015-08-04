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
#include "../utility/parallel_break_point.h"
#include "../utility/array_storage.h"
#include "../utility/time.h"
#include "../utility/json_input.h"

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
#include "../integrals/make_engine.h"
#include "../integrals/integral_engine_pool.h"
#include "../integrals/sparse_task_integrals.h"

#include "../scf/soad.h"
#include "../scf/diagonalize_for_coffs.hpp"
#include "../scf/clusterd_coeffs.h"

#include "../ta_routines/array_to_eigen.h"

using namespace tcc;
namespace ints = integrals;

int try_main(int argc, char *argv[]) {
    auto &world = madness::initialize(argc, argv);

    rapidjson::Document in;
    parse_input(argc, argv, in);

    if (!in.HasMember("xyz file") || !in.HasMember("number of bs clusters")
        || !in.HasMember("number of dfbs clusters")) {
        if (world.rank() == 0) {
            std::cout << "At a minimum your input file must provide\n";
            std::cout << "\"xyz file\", which is path to an xyz input\n";
            std::cout << "\"number of bs clusters\", which is the number of "
                         "clusters in the obs\n";
            std::cout << "\"number of dfbs clusters\", which is the number of "
                         "clusters in the dfbs\n";
        }
    }
    // Get necessary info
    std::string mol_file = in["xyz file"].GetString();
    int bs_nclusters = in["number of bs clusters"].GetInt();
    int dfbs_nclusters = in["number of dfbs clusters"].GetInt();
    int occ_nclusters = in.HasMember("number of occupied clusters")
                              ? in["number of occupied clusters"].GetInt()
                              : 1;

    // Get basis info
    std::string basis_name = in.HasMember("basis") ? in["basis"].GetString()
                                                   : "cc-pvdz";
    std::string df_basis_name = in.HasMember("df basis")
                                      ? in["df basis"].GetString()
                                      : "cc-pvdz-ri";

    // Get thresh info
    auto threshold = in.HasMember("block sparse threshold")
                           ? in["block sparse threshold"].GetDouble()
                           : 1e-13;
    auto low_rank_threshold = in.HasMember("low rank threshold")
                                    ? in["low rank threshold"].GetDouble()
                                    : 1e-8;

    volatile int debug
          = in.HasMember("debug break") ? in["debug break"].GetInt() : 0;
    utility::parallal_break_point(world, debug);

    if (world.rank() == 0) {
        std::cout << "Mol file is " << mol_file << std::endl;
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
    auto charge = 0;
    auto occupation = mol.occupation(charge);
    auto repulsion_energy = mol.nuclear_repulsion();

    utility::print_par(world, "Nuclear repulsion_energy = ", repulsion_energy,
                       "\n");

    /* auto bs_clusters = molecule::attach_hydrogens_kmeans(mol, bs_nclusters); */
    /* auto dfbs_clusters = molecule::attach_hydrogens_kmeans(mol, dfbs_nclusters); */
    auto bs_clusters = molecule::kmeans(mol, bs_nclusters);
    auto dfbs_clusters = molecule::kmeans(mol, dfbs_nclusters);

    basis::BasisSet bs{basis_name};
    basis::BasisSet df_bs{df_basis_name};

    std::streambuf *cout_sbuf = std::cout.rdbuf(); // Silence libint printing.
    std::ofstream fout("/dev/null");
    std::cout.rdbuf(fout.rdbuf());
    basis::Basis basis{bs.create_basis(bs_clusters)};
    basis::Basis df_basis{df_bs.create_basis(dfbs_clusters)};
    std::cout.rdbuf(cout_sbuf);

    world.gop.fence();
    std::cout << "Libint Init" << std::endl;

    libint2::init();

    // Make a btas to decomposed tensor function
    struct convert_2d {
        double cut_;
        convert_2d(double thresh) : cut_{thresh} {}
        using TileType = tensor::Tile<tensor::DecomposedTensor<double>>;
        TileType operator()(tensor::ShallowTensor<2> const &bt) {
            auto range = bt.range();

            auto const extent = range.extent();
            const auto i = extent[0];
            const auto j = extent[1];
            auto local_range = TA::Range{i, j};

            auto tensor = TA::Tensor<double>(local_range);
            const auto b_data = bt.tensor().data();
            const auto size = bt.tensor().size();
            std::copy(b_data, b_data + size, tensor.data());

            auto dense
                  = tensor::DecomposedTensor<double>(cut_, std::move(tensor));

            return tensor::Tile<tensor::DecomposedTensor<double>>(
                  range, std::move(dense));
        }
    };

    auto bs_basis_array = utility::make_array(basis, basis);

    /* // Begin Two electron integrals section. */
    auto eri_pool = ints::make_pool(ints::make_2body(basis, df_basis));

    auto to_decomp_with_decompose = [=](TA::Tensor<double> const &t) {
        auto range = t.range();

        auto const extent = range.extent();
        const auto i = extent[0];
        const auto j = extent[1];
        auto local_range = TA::Range{i, j};

        auto tensor = TA::Tensor<double>(local_range, t.data());
        auto dense = tensor::DecomposedTensor<double>(low_rank_threshold,
                                                      std::move(tensor));

        auto test = tensor::algebra::two_way_decomposition(dense);
        if (!test.empty()) {
            dense = std::move(test);
        }

        return tensor::Tile<tensor::DecomposedTensor<double>>(range,
                                                              std::move(dense));
    };

    /* // Computing Eri2 */
    utility::print_par(world, "Starting 2 Center Integrals\n");
    TA::Array<double, 2, tensor::Tile<tensor::DecomposedTensor<double>>,
              TA::SparsePolicy> V_inv_oh;
    {
        utility::print_par(world, "\tStarting V\n");
        auto eri2 = ints::BlockSparseIntegrals(
              world, eri_pool, utility::make_array(df_basis, df_basis),
              integrals::compute_functors::BtasToTaTensor{});

        decltype(eri2) L_inv_TA;
        {
            utility::print_par(world, "\tReplicating to Eigen\n");
            auto eig_E2 = array_ops::array_to_eigen(eri2);
            Eig::LLT<decltype(eig_E2)> llt(eig_E2);
            eig_E2 = llt.matrixL();
            decltype(eig_E2) eig_L_inv = eig_E2.inverse();
            utility::print_par(world, "\tConverting back to TA\n");
            L_inv_TA = array_ops::eigen_to_array<TA::Tensor<double>>(
                  world, eig_L_inv, eri2.trange().data()[0],
                  eri2.trange().data()[1]);
        }

        utility::print_par(world, "\tDecomposing Tiles\n");
        V_inv_oh = TA::to_new_tile_type(L_inv_TA, to_decomp_with_decompose);
        utility::print_size_info(V_inv_oh, "V^{-1/2}");
    }
    decltype(V_inv_oh)::wait_for_lazy_cleanup(world, 60);

    struct convert_3d {
        double cut_;
        convert_3d(double thresh) : cut_{thresh} {}
        using TileType = tensor::Tile<tensor::DecomposedTensor<double>>;
        TileType operator()(tensor::ShallowTensor<3> const &bt) {
            auto range = bt.range();

            auto const extent = range.extent();
            const auto X = extent[0];
            const auto i = extent[1];
            const auto j = extent[2];
            auto local_range = TA::Range{X, i, j};

            auto tensor = TA::Tensor<double>(local_range);
            const auto b_data = bt.tensor().data();
            const auto size = bt.tensor().size();
            std::copy(b_data, b_data + size, tensor.data());

            auto dense
                  = tensor::DecomposedTensor<double>(cut_, std::move(tensor));

            auto test = tensor::algebra::two_way_decomposition(dense);
            if (!test.empty()) {
                dense = std::move(test);
            }

            return tensor::Tile<tensor::DecomposedTensor<double>>(
                  range, std::move(dense));
        }
    };

    // Compute center integrals
    world.gop.fence();
    utility::print_par(world, "\nStarting 3 Center Integrals\n");
    auto E0 = tcc_time::now();
    auto Xab = ints::BlockSparseIntegrals(
          world, eri_pool, utility::make_array(df_basis, basis, basis),
          convert_3d(low_rank_threshold));
    auto E1 = tcc_time::now();
    auto etime = tcc_time::duration_in_s(E0, E1);
    utility::print_par(world, "Time to compute 3 center integrals ", etime,
                       " s\n");
    utility::print_size_info(Xab, "E");
    decltype(Xab)::wait_for_lazy_cleanup(world, 60);
    world.gop.fence();
    libint2::cleanup();
    std::cout << "Starting E sleep" << std::endl;
    sleep(10);
    std::cout << "Ending E sleep" << std::endl;
    world.gop.fence();
    std::cout << "Cleaned up libint and fenced" << std::endl;
    world.gop.fence();

    // Make B tensor
    auto B0 = tcc_time::now();
    Xab("X,i,j") = V_inv_oh("X,P") * Xab("P,i,j");
    Xab.truncate();
    auto B1 = tcc_time::now();
    auto btime = tcc_time::duration_in_s(B0, B1);
    utility::print_par(world, "\nTime to compute B ", btime, " s\n");
    utility::print_size_info(Xab, "B Tensor");
    decltype(Xab)::wait_for_lazy_cleanup(world, 60);
    world.gop.fence();
    std::cout << "Finished B and E should be gone" << std::endl;
    world.gop.fence();
    std::cout << "Starting sleep" << std::endl;
    sleep(10);
    std::cout << "Ending sleep" << std::endl;
    madness::finalize();
    return 0;
}

int main(int argc, char **argv) {
    try {
        try_main(argc, argv);
    } catch (const madness::MadnessException &e) {
        std::cout << "MADNESS exception: " << e.what() << std::endl;
    } catch (const TiledArray::Exception &e) {
        std::cout << "TA exception: " << e.what() << std::endl;
    } catch (const std::invalid_argument &e) {
        std::cout << "std::invalid_argument: " << e.what() << std::endl;
        return 1;
    } catch (const std::exception &e) {
        std::cout << "std::exception: " << e.what() << std::endl;
    } catch (...) {
        std::cout << "unknown exception" << std::endl;
    }
    return 0;
}

