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
#include "../integrals/scf/soad.h"

#include "../purification/sqrt_inv.h"
#include "../purification/purification_devel.h"

using namespace tcc;
namespace ints = integrals;

void main_print_clusters(
      std::vector<std::shared_ptr<molecule::Cluster>> const &bs,
      std::ostream &os) {
    std::vector<std::vector<molecule::Atom>> clusters;
    auto total_atoms = 0ul;
    for (auto const &cluster : bs) {
        std::vector<molecule::Atom> temp;
        for (auto atom : molecule::collapse_to_atoms(*cluster)) {
            temp.push_back(std::move(atom));
            ++total_atoms;
        }
        clusters.push_back(std::move(temp));
    }

    os << total_atoms << std::endl;
    os << "Whole molecule" << std::endl;

    for (auto const &cluster : clusters) {
        for (auto const &atom : cluster) {
            auto center = 0.52917721092 * atom.center();
            os << atom.charge() << " " << center[0] << " " << center[1] << " "
               << center[2] << std::endl;
        }
    }
    os << std::endl;
    auto counter = 0ul;
    for (auto const &cluster : clusters) {
        os << cluster.size() << std::endl;
        os << "Cluster " << counter++ << std::endl;
        for (auto const &atom : cluster) {
            auto center = 0.52917721092 * atom.center();
            os << atom.charge() << " " << center[0] << " " << center[1] << " "
               << center[2] << std::endl;
        }
        os << std::endl;
    }
}

int main(int argc, char *argv[]) {
    auto &world = madness::initialize(argc, argv);

    std::string mol_file = (argc >= 2) ? argv[1] : "";
    int bs_nclusters = (argc >= 3) ? std::stoi(argv[2]) : 0;
    int dfbs_nclusters = (argc >= 4) ? std::stoi(argv[3]) : 0;

    if (mol_file.empty() || 0 == bs_nclusters || 0 == dfbs_nclusters) {
        std::cout << "input is $./program mol_file "
                     "bs_cluster dfbs_clusters "
                     "basis_set(cc-pvdz) df_basis_set(cc-pvdz-ri) "
                     "sparse_threshold(1e-11) "
                     "low_rank_threshhold(1e-7) print_cluster_xyz(true)\n";
        return 0;
    }

    std::string basis_name = (argc >= 5) ? argv[4] : "cc-pvdz";
    std::string df_basis_name = (argc >= 6) ? argv[5] : "cc-pvdz-ri";

    auto threshold = (argc >= 7) ? std::stod(argv[6]) : 1e-11;
    auto low_rank_threshold = (argc >= 8) ? std::stod(argv[7]) : 1e-7;
    bool print_clusters = (argc >= 9) ? std::stod(argv[8]) : true;

    if (world.rank() == 0) {
        std::cout << "Mol file is " << mol_file << std::endl;
        std::cout << "basis is " << basis_name << std::endl;
        std::cout << "df basis is " << df_basis_name << std::endl;
        std::cout << "Using " << bs_nclusters << " bs clusters" << std::endl;
        std::cout << "Using " << dfbs_nclusters << " dfbs clusters"
                  << std::endl;
        std::cout << "low rank threshhold is " << low_rank_threshold
                  << std::endl;
        if (print_clusters) {
            std::cout << "Printing clusters to clusters_bs.xyz and "
                         "cluster_dfbs.xyz." << std::endl;
        }
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

    auto bs_clusters = molecule::attach_hydrogens_kmeans(mol, bs_nclusters);
    auto dfbs_clusters = molecule::attach_hydrogens_kmeans(mol, dfbs_nclusters);

    if (world.rank() == 0) {
        if (print_clusters) {
            std::cout << "Printing clusters\n";
            std::ofstream bs_cluster_file("clusters_bs.xyz");
            main_print_clusters(bs_clusters, bs_cluster_file);
            bs_cluster_file.close();
            std::ofstream dfbs_cluster_file("clusters_dfbs.xyz");
            main_print_clusters(dfbs_clusters, dfbs_cluster_file);
            dfbs_cluster_file.close();
        }
    }
    world.gop.fence();

    basis::BasisSet bs{basis_name};
    basis::BasisSet df_bs{df_basis_name};

    std::streambuf *cout_sbuf = std::cout.rdbuf(); // Silence libint printing.
    std::ofstream fout("/dev/null");
    std::cout.rdbuf(fout.rdbuf());
    basis::Basis basis{bs.create_basis(bs_clusters)};
    basis::Basis df_basis{df_bs.create_basis(dfbs_clusters)};
    std::cout.rdbuf(cout_sbuf);


    libint2::init();

    // Make a btas to decomposed tensor function
    struct convert_2d {
        double cut_;
        convert_2d(double thresh) : cut_{thresh} {}
        using TileType = tensor::Tile<tensor::DecomposedTensor<double>>;
        TileType operator()(tensor::ShallowTensor<2> const &bt) {
            auto range = bt.range();

            auto const &extent = range.size();
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

    // Compute overlap.
    auto overlap_pool = ints::make_pool(ints::make_1body("overlap", basis));
    auto S = BlockSparseIntegrals(world, overlap_pool, bs_basis_array,
                                  convert_2d(low_rank_threshold));

    auto to_ta = [](tensor::Tile<tensor::DecomposedTensor<double>> t) {
        auto tensor = tensor::algebra::combine(t.tile());
        auto range = t.range();
        return TA::Tensor<double>(range, tensor.data());
    };
    auto S_TA = TA::to_new_tile_type(S, to_ta);

    // Invert overlap
    utility::print_par(world, "\nComputing overlap inverse\n");
    auto S_inv_sqrt = pure::inverse_sqrt(S_TA);

    // Compute T
    auto kinetic_pool = ints::make_pool(ints::make_1body("kinetic", basis));
    auto T = BlockSparseIntegrals(world, kinetic_pool, bs_basis_array,
                                  convert_2d(low_rank_threshold));

    /* // Compute V */
    auto nuclear_pool = ints::make_pool(ints::make_1body("nuclear", basis));
    auto V = BlockSparseIntegrals(world, nuclear_pool, bs_basis_array,
                                  convert_2d(low_rank_threshold));

    /* // Compute Hcore */
    utility::print_par(world, "Computing Hcore\n");
    decltype(V) H;
    H("i,j") = T("i,j") + V("i,j");
    world.gop.fence();

    auto H_TA = TA::to_new_tile_type(H, to_ta);


    auto to_decomp = [=](TA::Tensor<double> const &t) {
        auto range = t.range();

        auto const &extent = range.size();
        const auto i = extent[0];
        const auto j = extent[1];
        auto local_range = TA::Range{i, j};

        auto tensor = TA::Tensor<double>(local_range, t.data());
        auto dense = tensor::DecomposedTensor<double>(low_rank_threshold,
                                                      std::move(tensor));

        return tensor::Tile<tensor::DecomposedTensor<double>>(range,
                                                              std::move(dense));
    };

    /* // Begin Two electron integrals section. */
    auto eri_pool = ints::make_pool(ints::make_2body(basis, df_basis));

    /* // Computing Eri2 */
    utility::print_par(world, "\n");
    auto eri2 = ints::BlockSparseIntegrals(
          world, eri_pool, utility::make_array(df_basis, df_basis),
          integrals::compute_functors::BtasToTaTensor{});


    /* // Computing the sqrt inverse of Eri2 */
    utility::print_par(world, "\nComputing eri2 sqrt Inverse\n");
    auto inv_timer
          = tcc_time::make_timer([&]() { return pure::inverse_sqrt(eri2); });

    auto eri2_sqrt_inv = inv_timer.apply();
    utility::print_par(world, "Eri2 inverse computation time = ",
                       inv_timer.time(), "\n");
    utility::print_size_info(eri2_sqrt_inv, "Eri2 sqrt inverse");
    decltype(eri2_sqrt_inv) V_inv_TA;
    V_inv_TA("i,j") = eri2_sqrt_inv("i,k") * eri2_sqrt_inv("k,j");

    auto to_decomp_with_decompose = [=](TA::Tensor<double> const &t) {
        auto range = t.range();

        auto const &extent = range.size();
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
    auto V_inv = TA::to_new_tile_type(V_inv_TA, to_decomp_with_decompose);
    utility::print_size_info(V_inv, "V_inv");


    struct convert_3d {
        double cut_;
        convert_3d(double thresh) : cut_{thresh} {}
        using TileType = tensor::Tile<tensor::DecomposedTensor<double>>;
        TileType operator()(tensor::ShallowTensor<3> const &bt) {
            auto range = bt.range();

            auto const &extent = range.size();
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
    utility::print_par(world, "\n");
    auto Xab = ints::BlockSparseIntegrals(
          world, eri_pool, utility::make_array(df_basis, basis, basis),
          convert_3d(low_rank_threshold));

    utility::print_size_info(Xab, "E");


    decltype(Xab) W;
    W("X,i,j") = V_inv("X,P") * Xab("P,i,j");
    W.truncate();

    utility::print_size_info(W, "W");

    decltype(H) F;
    F = ints::scf::fock_from_minimal(world, basis, df_basis, eri_pool, V_inv, H,
                                     W, bs_clusters, low_rank_threshold,
                                     convert_3d(low_rank_threshold));
    auto F_TA = TA::to_new_tile_type(F, to_ta);

    auto purifier = pure::make_orthogonal_tr_reset_pure(S_inv_sqrt);
    auto D_TA = purifier(F_TA, occupation);
    auto energy = D_TA("i,j").dot(F_TA("i,j") + H_TA("i,j"), world).get();
    std::cout << "Initial energy = " << energy + repulsion_energy << std::endl;
    auto D = to_new_tile_type(D_TA, to_decomp);

    decltype(H) J, K;
    utility::print_par(world, "\nStarting SCF iterations\n");
    auto diis = TiledArray::DIIS<decltype(D_TA)>{3, 7};
    auto iter = 1;
    decltype(F_TA) Ferror;
    auto error = 1.0;
    const auto volume = double(F.trange().elements().volume());
    decltype(D) D_old = D;
    decltype(D) D_diff;
    while (error >= 1e-12 && iter <= 35) {
        auto t0 = tcc_time::now();
        D = to_new_tile_type(D_TA, to_decomp);
        auto j0 = tcc_time::now();
        J("i,j") = W("X,i,j") * (Xab("X,a,b") * D("a,b"));
        auto j1 = tcc_time::now();

        auto k0 = tcc_time::now();
        if (iter <= 5) {
            K("i,j") = W("X,a,i") * (Xab("X,j,b") * D("b,a"));
        } else {
            D_diff("i,j") = D("i,j") - D_old("i,j");
            D_diff.truncate();
            K("i,j") = K("i,j") + W("X,a,i") * (Xab("X,j,b") * D_diff("b,a"));
        }
        auto k1 = tcc_time::now();
        D_old = D;

        F("i,j") = H("i,j") + 2 * J("i,j") - K("i,j");

        F_TA = TA::to_new_tile_type(F, to_ta);

        Ferror("i,j") = F_TA("i,k") * D_TA("k,l") * S_TA("l,j")
                        - S_TA("i,k") * D_TA("k,l") * F_TA("l,j");

        error = Ferror("i,j").norm().get() / volume;
        diis.extrapolate(F_TA, Ferror);

        D_TA = purifier(F_TA, occupation);
        energy = D_TA("i,j").dot(F_TA("i,j") + H_TA("i,j"), world).get();

        auto t1 = tcc_time::now();
        if (iter % 5 == 0) {
            decltype(Xab) Xtemp;
            Xtemp("X,a,i") = Xab("X,i,b") * D("b,a");
            utility::print_par(world, "\nPrinting size info for Xtemp, iter ",
                               iter, "\n");
            utility::print_size_info(Xtemp, "Xtemp");
            utility::print_par(world, "\n");
            world.gop.fence();
        }

        auto time = tcc_time::duration_in_s(t0, t1);
        auto jtime = tcc_time::duration_in_s(j0, j1);
        auto ktime = tcc_time::duration_in_s(k0, k1);
        utility::print_par(world, "Iteration: ", iter++, " has energy ",
                           std::setprecision(11), energy + repulsion_energy,
                           " with error ", error, " in ", time, " s \n");
        utility::print_par(world, "\tJ time ", jtime, " s\n\tK time ", ktime,
                           " s\n");
    }

    utility::print_par(world, "Final energy = ", std::setprecision(11),
                       energy + repulsion_energy, "\n");

    world.gop.fence();
    libint2::cleanup();
    madness::finalize();
    return 0;
}
