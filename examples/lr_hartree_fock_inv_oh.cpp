#include <memory>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <chrono>

#include "../include/libint.h"
#include "../include/tiledarray.h"
#include "../include/btas.h"

#include "../utility/make_array.h"
#include "../utility/parallel_print.h"
#include "../utility/parallel_break_point.h"
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

#include "../density/sqrt_inv.h"
#include "../density/purification_devel.h"

#include "../array_ops/array_to_eigen.h"

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
    bool print_clusters = (argc >= 9) ? std::stoi(argv[8]) : true;
    volatile int debug = (argc >= 10) ? std::stod(argv[9]) : 0;
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

    if (world.rank() == 0) {
        std::cout << "Basis trange " << std::endl;
        TA::TiledRange1 bs_range = basis.create_flattend_trange1();
        std::cout << bs_range << std::endl;
        for (auto range : bs_range) {
            std::cout << range.first << " " << range.second << std::endl;
        }
        TA::TiledRange1 dfbs_range = df_basis.create_flattend_trange1();
        std::cout << "DF Basis trange " << std::endl;
        std::cout << dfbs_range << std::endl;
        for (auto range : dfbs_range) {
            std::cout << range.first << " " << range.second << std::endl;
        }
    }

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

    // Compute overlap.
    auto overlap_pool = ints::make_pool(ints::make_1body("overlap", basis));
    auto S = BlockSparseIntegrals(world, overlap_pool, bs_basis_array,
                                  convert_2d(low_rank_threshold));

    auto to_ta = [](tensor::Tile<tensor::DecomposedTensor<double>> const &t) {
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

        auto const extent = range.extent();
        const auto i = extent[0];
        const auto j = extent[1];
        auto local_range = TA::Range{i, j};

        auto tensor = TA::Tensor<double>(local_range, t.data());
        auto dense = tensor::DecomposedTensor<double>(low_rank_threshold,
                                                      std::move(tensor));

        return tensor::Tile<tensor::DecomposedTensor<double>>(range,
                                                              std::move(dense));
    };
    decltype(S_TA) sqrt_S_TA;
    sqrt_S_TA("i,j") = S_inv_sqrt("i,k") * S_TA("k,j");
    sqrt_S_TA.truncate();
    auto sqrt_S = TA::to_new_tile_type(sqrt_S_TA, to_decomp);
    auto S_inv_sqrt_lr = TA::to_new_tile_type(S_inv_sqrt, to_decomp);

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
    auto V_inv_oh
          = TA::to_new_tile_type(eri2_sqrt_inv, to_decomp_with_decompose);
    utility::print_size_info(V_inv_oh, "V^{-1/2}");


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
    utility::print_par(world, "\n");
    auto Xab = ints::BlockSparseIntegrals(
          world, eri_pool, utility::make_array(df_basis, basis, basis),
          convert_3d(low_rank_threshold));
    utility::print_size_info(Xab, "E");

    Xab("X,i,j") = V_inv_oh("X,P") * Xab("P,i,j");
    Xab.truncate();

    utility::print_size_info(Xab, "Xab with V^{-1/2}");
    utility::print_par(world, "\n");

    decltype(H) F;
    F = ints::scf::fock_from_minimal_v_oh(world, basis, df_basis, eri_pool, H,
                                          V_inv_oh, Xab, bs_clusters,
                                          low_rank_threshold,
                                          convert_3d(low_rank_threshold));
    auto F_TA = TA::to_new_tile_type(F, to_ta);

    auto purifier = pure::make_orthogonal_tr_reset_pure(S_inv_sqrt);
    auto D_TA = purifier(F_TA, occupation);
    auto energy = D_TA("i,j").dot(F_TA("i,j") + H_TA("i,j"), world).get();
    utility::print_par(world, "Initial energy = ", energy + repulsion_energy,
                       "\n");
    auto D = to_new_tile_type(D_TA, to_decomp);

    struct compress {
        double cut_;
        compress(double thresh) : cut_{thresh} {}
        using TileType = tensor::Tile<tensor::DecomposedTensor<double>>;
        double operator()(TileType &result) {
            if (result.tile().ndecomp() == 1) {
                auto test
                      = tensor::algebra::two_way_decomposition(result.tile());
                if (!test.empty()) {
                    result.tile() = std::move(test);
                }
            } else {
                tensor::algebra::recompress(result.tile());
            }

            return result.norm();
        }
    };

    auto nblocks = (dfbs_nclusters < (occupation / 2)) ? dfbs_nclusters
                                                       : occupation / 2;
    auto block_size = std::max(std::size_t((occupation / 2) / nblocks), 1ul);
    std::vector<std::size_t> blocks;
    blocks.reserve(nblocks + 1);
    blocks.push_back(0);
    for (auto i = block_size; i < occupation / 2; i += block_size) {
        blocks.push_back(i);
    }
    blocks.push_back(occupation / 2);
    auto tr1 = TA::TiledRange1(blocks.begin(), blocks.end());
    auto counter = 0;
    if (world.rank() == 0) {
        std::cout << "Range for occupied orbital\n";
        for (auto range : tr1) {
            std::cout << counter++ << ": " << range.first << " " << range.second
                      << std::endl;
        }
        std::cout << std::endl;
    }

    decltype(D) L;
    decltype(H) J, K;
    utility::print_par(world, "\nStarting SCF iterations\n");
    auto diis = TiledArray::DIIS<decltype(D_TA)>(1);
    auto iter = 1;
    decltype(F_TA) Ferror;
    auto error = 1.0;
    const auto volume = double(F.trange().elements().volume());
    decltype(D) D_old = D;
    decltype(D) D_diff;
    double d_diff_norm = 0.0;
    double time;
    double puretime, ktime, jtime;
    while (error >= 1e-12 && iter <= 35) {
        utility::print_par(world, "Iteration: ", iter, "\n");
        auto t0 = tcc_time::now();
        D = to_new_tile_type(D_TA, to_decomp);
        utility::print_par(world, "\tStarting Coulomb...\n");
        auto j0 = tcc_time::now();
        J("i,j") = Xab("X,i,j") * (Xab("X,a,b") * D("a,b"));
        auto j1 = tcc_time::now();

        utility::print_par(world, "\tStarting Exchange...\n");
        auto k0 = tcc_time::now();
        {
            decltype(D) D_orth;
            D_orth("i,j") = sqrt_S("i,k") * D("k,l") * sqrt_S("l,j");
            world.gop.fence();
            auto Eig_D = array_ops::array_to_eigen(D_orth);
            Eig::LDLT<decltype(Eig_D)> ldl(Eig_D);
            Eig::VectorXd d_vec = ldl.vectorD();
            bool use_cholesky = true;
            for (auto i = 0ul; i < occupation / 2; ++i) {
                if (d_vec[i] <= 0) {
                    use_cholesky = false;
                    break;
                }
            }

            array_ops::Matrix<double> L_eig;
            if (use_cholesky) {
                utility::print_par(world, "\tUsing LDLT for density\n");
                array_ops::Matrix<double> P
                      = ldl.transpositionsP()
                        * decltype(Eig_D)::Identity(Eig_D.rows(), Eig_D.cols());
                P.transposeInPlace();

                L_eig = P
                        * array_ops::Matrix<double>(ldl.matrixL())
                                .leftCols(occupation / 2);

                for (auto i = 0ul; i < occupation / 2; ++i) {
                    d_vec[i] = std::sqrt(d_vec[i]);
                }

                array_ops::Matrix<double> diag
                      = array_ops::Matrix<double>(d_vec.asDiagonal())
                              .block(0, 0, occupation / 2, occupation / 2);
                L_eig *= diag;
            } else {
                utility::print_par(world, "\tUsing EVD for density\n");
                Eig::SelfAdjointEigenSolver<decltype(Eig_D)> es(Eig_D);
                L_eig = es.eigenvectors().rightCols(occupation / 2);

                auto diag
                      = array_ops::Matrix<double>(es.eigenvalues().asDiagonal())
                              .bottomRightCorner(occupation / 2,
                                                 occupation / 2);

                for (auto i = 0; i < diag.cols(); ++i) {
                    diag(i, i) = std::sqrt(diag(i, i));
                }

                L_eig = L_eig * diag;
            }

            TA::TiledRange1 tr0 = D.trange().data().front();
            L = array_ops::
                  eigen_to_array<tensor::
                                       Tile<tensor::DecomposedTensor<double>>>(
                        world, L_eig, tr0, tr1);
            L("i,j") = S_inv_sqrt_lr("i,k") * L("k,j");
        }

        decltype(Xab) Eai;
        Eai("X,i,a") = Xab("X,a,b") * L("b,i");
        K("a,b") = Eai("X,i,a") * Eai("X,i,b");
        auto k1 = tcc_time::now();

        D_old = D;

        F("i,j") = H("i,j") + 2 * J("i,j") - K("i,j");


        F_TA = TA::to_new_tile_type(F, to_ta);
        F = to_new_tile_type(F_TA, to_decomp);

        Ferror("i,j") = F_TA("i,k") * D_TA("k,l") * S_TA("l,j")
                        - S_TA("i,k") * D_TA("k,l") * F_TA("l,j");

        error = Ferror("i,j").norm().get() / volume;
        diis.extrapolate(F_TA, Ferror);


        utility::print_par(world, "\tStarting Purification...\n");
        auto td0 = tcc_time::now();
        D_TA = purifier(F_TA, occupation);
        energy = D_TA("i,j").dot(F_TA("i,j") + H_TA("i,j"), world).get();
        auto td1 = tcc_time::now();

        auto t1 = tcc_time::now();
        time = tcc_time::duration_in_s(t0, t1);
        jtime = tcc_time::duration_in_s(j0, j1);
        ktime = tcc_time::duration_in_s(k0, k1);
        puretime = tcc_time::duration_in_s(td0, td1);
        utility::print_par(world, "\tHas energy ", std::setprecision(14),
                           energy + repulsion_energy, " with error ", error,
                           " in ", time, " s \n");
        utility::print_par(world, "\tJ time ", jtime, " s\n\tK time ", ktime,
                           " s\n\tPure time ", puretime, "\n");
        ++iter;
    }

    utility::print_par(world, "Final energy = ", std::setprecision(11),
                       energy + repulsion_energy, "\n");

    utility::print_par(world, "\n\n\nTEST REGION FOLLOWS\n\n\n");

    { // TEST REGION
        decltype(Xab) Eai;
        auto te0 = tcc_time::now();
        Eai("X,i,a") = Xab("X,a,b") * L("b,i");
        Eai.truncate();
        auto te1 = tcc_time::now();
        auto timee = tcc_time::duration_in_s(te0, te1);
        utility::print_par(world, "Made Eia in ", timee, " s\n");
        utility::print_size_info(Eai, "Eia");

        auto tec0 = tcc_time::now();
        TA::foreach_inplace(Eai, compress(low_rank_threshold));
        auto tec1 = tcc_time::now();
        auto timeec = tcc_time::duration_in_s(tec0, tec1);
        utility::print_par(world, "Compressed Eai in ", timeec, " s\n");
        utility::print_size_info(Eai, "Eai");

        utility::print_par(world, "\nForming K\n");
        decltype(K) Ktest;
        auto t0 = tcc_time::now();
        Ktest("a,b") = Eai("X,i,a") * Eai("X,i,b");
        auto t1 = tcc_time::now();
        auto time = tcc_time::duration_in_s(t0, t1);
        Ktest("i,j") = Ktest("i,j") - K("i,j");
        auto diff = Ktest("i,j").norm().get();

        utility::print_par(world, "\nK formation time (Eia^T * Eia) ", time,
                           " s with norm diff ", diff, "\n");
    } // End Test Region

    { // Begin MP2
        utility::print_par(world, "\nBegining MP2\n");
        auto F_eig = array_ops::array_to_eigen(F);
        auto S_eig = array_ops::array_to_eigen(S);
        Eig::GeneralizedSelfAdjointEigenSolver<decltype(S_eig)> es(F_eig,
                                                                   S_eig);
        Eig::VectorXd evals = es.eigenvalues();
        decltype(S_eig) C_occ = es.eigenvectors().leftCols(occupation / 2);
        decltype(S_eig) C_vir
              = es.eigenvectors().rightCols(S_eig.rows() - occupation / 2);

        auto nblocks = (dfbs_nclusters < (S_eig.rows() - occupation / 2))
                             ? dfbs_nclusters
                             : S_eig.rows() - occupation / 2;
        auto block_size
              = std::max(std::size_t((S_eig.rows() - occupation / 2) / nblocks),
                         1ul);
        std::vector<std::size_t> blocks;
        blocks.reserve(nblocks + 1);
        blocks.push_back(0);
        for (auto i = block_size; i < S_eig.rows() - occupation / 2;
             i += block_size) {
            blocks.push_back(i);
        }
        blocks.push_back(S_eig.rows() - occupation / 2);
        auto tr_vir = TA::TiledRange1(blocks.begin(), blocks.end());

        TA::TiledRange1 tr0 = D.trange().data().front();
        auto Ci = array_ops::
              eigen_to_array<tensor::Tile<tensor::DecomposedTensor<double>>>(
                    world, C_occ, tr0, tr1);
        auto Cv = array_ops::
              eigen_to_array<tensor::Tile<tensor::DecomposedTensor<double>>>(
                    world, C_vir, tr0, tr_vir);

        decltype(Xab) Xia;
        Xia("X,i,a") = Xab("X,mu,nu") * Ci("nu,i") * Cv("mu,a");

        utility::print_size_info(Xia, "Xia");

        auto Xia_TA = TA::to_new_tile_type(Xia, to_ta);
        TA::Array<double, 4, TA::Tensor<double>, TA::SparsePolicy> IAJB;
        IAJB("i,a,j,b") = Xia_TA("X,i,a") * Xia_TA("X,j,b");
        utility::print_size_info(IAJB, "IAJB");
        auto vec_ptr = std::make_shared<Eig::VectorXd>(std::move(evals));
        struct Mp2Red {
            using result_type = double;
            using argument_type = TA::Tensor<double>;

            std::shared_ptr<Eig::VectorXd> vec_;
            unsigned int n_occ_;

            Mp2Red(std::shared_ptr<Eig::VectorXd> vec, int n_occ)
                    : vec_(std::move(vec)), n_occ_(n_occ) {}
            Mp2Red(Mp2Red const &) = default;

            result_type operator()() const { return 0.0; }
            result_type operator()(result_type const &t) const { return t; }
            void operator()(result_type &me, result_type const &other) const {
                me += other;
            }

            void operator()(result_type &me, argument_type const &tile) const {
                auto const &range = tile.range();
                auto const &vec = *vec_;
                auto const st = range.start();
                auto const fn = range.finish();
                auto tile_idx = 0;
                for (auto i = st[0]; i < fn[0]; ++i) {
                    const auto e_i = vec[i];
                    for (auto a = st[1]; a < fn[1]; ++a) {
                        const auto e_ia = e_i - vec[a + n_occ_];
                        for (auto j = st[2]; j < fn[2]; ++j) {
                            const auto e_iaj = e_ia + vec[j];
                            for (auto b = st[3]; b < fn[3]; ++b, ++tile_idx)
                            {
                                const auto e_iajb = e_iaj - vec[b + n_occ_];
                                me += 1 / (e_iajb)*tile.data()[tile_idx];
                            }
                        }
                    }
                }
                // for (auto i = 0ul; i < range.volume(); ++i) {
                //     auto idx = range.idx(i);
                //     // vals are ordered i, a, j, b
                //     // need e_i + e_j - e_a - e_b
                //     auto ival = vec[idx[0]];
                //     auto aval = vec[n_occ_ + idx[1]];
                //     auto jval = vec[idx[2]];
                //     auto bval = vec[n_occ_ + idx[3]];
                //     auto total = 1.0 / (ival + jval - aval - bval);
                //     me += total * tile.data()[i];
                // }
            }
        };

        double energy_mp2
              = (IAJB("i,a,j,b") * (2 * IAJB("i,a,j,b") - IAJB("i,b,j,a")))
                      .reduce(Mp2Red(vec_ptr, occupation / 2));

        utility::print_par(world, "MP2 energy = ", energy_mp2,
                           " total energy = ",
                           energy + energy_mp2 + repulsion_energy, "\n");
    }

    world.gop.fence();
    libint2::cleanup();
    madness::finalize();
    return 0;
}
