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
#include "../cc/ccsd.h"
#include "../cc/two_electron_int_mo.h"
#include "../mp2/trange1_engine.h"
#include "../mp2/mp2.h"
#include "../ta_routines/array_to_eigen.h"

using namespace tcc;
namespace ints = integrals;

static std::map<int, std::string> atom_names = { {1 , "H"}
                                             , {2 , "He"}
                                             , {3 , "Li"}
                                             , {4 , "Be"}
                                             , {5 , "B"}
                                             , {6 , "C"}
                                             , {7 , "N"}
                                             , {8 , "O"}
                                             , {9 , "F"}
                                             , {10 , "Ne"}
                                             , {11 , "Na"}
                                             , {12 , "Mg"}
};

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
            os << atom_names[atom.charge()] << " " << center[0] << " " << center[1] << " "
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
            os << atom_names[atom.charge()] << " " << center[0] << " " << center[1] << " "
               << center[2] << std::endl;
        }
        os << std::endl;
    }
}

int main(int argc, char *argv[]) {
    auto &world = madness::initialize(argc, argv);


    // declare variables needed for ccsd
    std::shared_ptr<tcc::cc::TwoElectronIntMO<TA::Tensor<double>, TA::SparsePolicy>> g;
    std::shared_ptr<tcc::TRange1Engine> tre;
    Eigen::MatrixXd ens;
    TA::Array<double, 2, TA::Tensor<double>, TA::SparsePolicy> fock_mo;
    {
        std::string mol_file = (argc >= 2) ? argv[1] : "";
        std::size_t blocksize = (argc >= 3) ? std::stoi(argv[2]) : 16;
        int bs_nclusters = (argc >= 4) ? std::stoi(argv[3]) : 0;
        int dfbs_nclusters = (argc >= 5) ? std::stoi(argv[4]) : 0;

        if (mol_file.empty() || 0 == bs_nclusters || 0 == dfbs_nclusters) {
            std::cout << "input is $./program mol_file blocksize"
                    "bs_cluster dfbs_clusters "
                    "basis_set(cc-pvdz) df_basis_set(cc-pvdz-ri) "
                    "sparse_threshold(1e-11) "
                    "low_rank_threshhold(1e-7) print_cluster_xyz(true)\n";
            return 0;
        }

        std::string basis_name = (argc >= 6) ? argv[5] : "cc-pvdz";
        std::string df_basis_name = (argc >= 7) ? argv[6] : "cc-pvdz-ri";

        auto threshold = (argc >= 8) ? std::stod(argv[7]) : 1e-11;
        auto low_rank_threshold = (argc >= 9) ? std::stod(argv[8]) : 1e-7;
        bool print_clusters = (argc >= 10) ? std::stoi(argv[9]) : true;
        volatile int debug = (argc >= 11) ? std::stod(argv[10]) : 0;
        utility::parallal_break_point(world, debug);

        if (world.rank() == 0) {
            std::cout << "Mol file is " << mol_file << std::endl;
            std::cout << "basis is " << basis_name << std::endl;
            std::cout << "df basis is " << df_basis_name << std::endl;
            std::cout << "Using " << bs_nclusters << " bs clusters" <<
            std::endl;
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

        utility::print_par(world, "Nuclear repulsion_energy = ",
                           repulsion_energy,
                           "\n");

        auto bs_clusters = molecule::attach_hydrogens_kmeans(mol, bs_nclusters);
        auto dfbs_clusters = molecule::attach_hydrogens_kmeans(mol,
                                                               dfbs_nclusters);

        world.gop.fence();
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
            TA::TiledRange1 dfbs_range = df_basis.create_flattend_trange1();
            std::cout << "DF Basis trange " << std::endl;
            std::cout << dfbs_range << std::endl;
        }

        libint2::init();

        // Make a btas to decomposed tensor function
        struct convert_2d {
            double cut_;

            convert_2d(double thresh) : cut_{thresh} { }

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
                        = tensor::DecomposedTensor<double>(cut_,
                                                           std::move(tensor));

                return tensor::Tile<tensor::DecomposedTensor<double>>(
                        range, std::move(dense));
            }
        };

        auto bs_basis_array = utility::make_array(basis, basis);

        // Compute overlap.
        auto overlap_pool = ints::make_pool(ints::make_1body("overlap", basis));
        auto S = BlockSparseIntegrals(world, overlap_pool, bs_basis_array,
                                      convert_2d(low_rank_threshold));

        auto to_ta = [](
                tensor::Tile<tensor::DecomposedTensor<double>> const &t) {
          auto tensor = tensor::algebra::combine(t.tile());
          auto range = t.range();
          return TA::Tensor<double>(range, tensor.data());
        };
        auto S_TA = TA::to_new_tile_type(S, to_ta);

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

        auto to_decomp = [=](TA::Tensor <double> const &t) {
          auto range = t.range();

          auto const extent = range.extent();
          const auto i = extent[0];
          const auto j = extent[1];
          auto local_range = TA::Range{i, j};

          auto tensor = TA::Tensor<double>(local_range, t.data());
          auto dense = tensor::DecomposedTensor<double>(low_rank_threshold,
                                                        std::move(tensor));

          return tensor::Tile<tensor::DecomposedTensor<double>>(range,
                                                                std::move(
                                                                        dense));
        };

        /* // Begin Two electron integrals section. */
        auto eri_pool = ints::make_pool(ints::make_2body(basis, df_basis));

        auto to_decomp_with_decompose = [=](TA::Tensor <double> const &t) {
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
                                                                std::move(
                                                                        dense));
        };
        /* // Computing Eri2 */
        utility::print_par(world, "Starting 2 Center Integrals\n");
        TA::Array <double, 2, tensor::Tile<tensor::DecomposedTensor<double>>,
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
                L_inv_TA = array_ops::eigen_to_array < TA::Tensor < double >> (
                        world, eig_L_inv, eri2.trange().data()[0],
                                eri2.trange().data()[1]);
            }

            utility::print_par(world, "\tDecomposing Tiles\n");
            V_inv_oh = TA::to_new_tile_type(L_inv_TA, to_decomp_with_decompose);
            utility::print_size_info(V_inv_oh, "V^{-1/2}");
        }
        decltype(V_inv_oh)::wait_for_lazy_cleanup(world, 2);


        struct convert_3d {
            double cut_;

            convert_3d(double thresh) : cut_{thresh} { }

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
                        = tensor::DecomposedTensor<double>(cut_,
                                                           std::move(tensor));

                auto test = tensor::algebra::two_way_decomposition(dense);
                if (!test.empty()) {
                    dense = std::move(test);
                }

                return tensor::Tile<tensor::DecomposedTensor<double>>(
                        range, std::move(dense));
            }
        };

        // Compute center integrals
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

        // Make B tensor
        auto B0 = tcc_time::now();
        Xab("X,i,j") = V_inv_oh("X,P") * Xab("P,i,j");
        Xab.truncate();
        auto B1 = tcc_time::now();
        auto btime = tcc_time::duration_in_s(B0, B1);
        utility::print_par(world, "\nTime to compute B ", btime,
                           " s\n");
        utility::print_size_info(Xab, "B Tensor");
        decltype(Xab)::wait_for_lazy_cleanup(world, 60);

        decltype(H) F;
        utility::print_par(world, "\nStarting SOAD guess");
        auto soad0 = tcc_time::now();
        F = ints::scf::fock_from_minimal_v_oh(world, basis, df_basis, eri_pool,
                                              H,
                                              V_inv_oh, Xab, bs_clusters,
                                              low_rank_threshold * 100,
                                              convert_3d(low_rank_threshold));
        auto soad1 = tcc_time::now();
        auto soad_time = tcc_time::duration_in_s(soad0, soad1);
        utility::print_par(world, "\nSoad time ", soad_time, " s\n");
        decltype(F)::wait_for_lazy_cleanup(world);

        utility::print_par(world, "\nConverting Fock to TA::Tensor...\n");
        auto F_TA = TA::to_new_tile_type(F, to_ta);
        world.gop.fence();

        auto n_occ = occupation / 2;
        auto tr_i = scf::tr_occupied(dfbs_nclusters, n_occ);
        utility::print_par(world, "Computing MO coeffs...\n");
        auto Coeffs_TA = scf::Coeffs_from_fock(F_TA, S_TA, tr_i, n_occ);
        utility::print_par(world, "Converting Coeffs to Decomp Form...\n");
        auto Coeffs = TA::to_new_tile_type(Coeffs_TA, to_decomp);

        decltype(Coeffs_TA) D_TA;
        utility::print_par(world, "Forming Density...\n");
        D_TA("i,j") = Coeffs_TA("i,a") * Coeffs_TA("j,a");

        utility::print_par(world, "Computing Initial energy...\n");
        auto energy = D_TA("i,j").dot(F_TA("i,j") + H_TA("i,j"), world).get();
        utility::print_par(world, "Initial energy = ",
                           energy + repulsion_energy,
                           "\n");

        auto D = to_new_tile_type(D_TA, to_decomp);

        decltype(D) J, K;
        utility::print_par(world, "\nStarting SCF iterations\n");
        auto diis = TiledArray::DIIS<decltype(D_TA)>(1);
        auto iter = 1;
        decltype(F_TA) Ferror;
        auto error = 1.0;
        auto old_e = energy;
        auto delta_e = energy;
        const auto volume = double(F.trange().elements().volume());
        decltype(Xab) Eai;
        double time;
        double ktime, jtime;
        while ((error >= 1e-13 || std::abs(delta_e) >= 1e-12) && iter <= 35) {
            utility::print_par(world, "Iteration: ", iter, "\n");
            auto t0 = tcc_time::now();
            D = to_new_tile_type(D_TA, to_decomp);

            utility::print_par(world, "\tStarting Coulomb...  ");
            auto j0 = tcc_time::now();
            J("i,j") = Xab("X,i,j") * (Xab("X,a,b") * D("a,b"));
            auto j1 = tcc_time::now();
            jtime = tcc_time::duration_in_s(j0, j1);
            utility::print_par(world, jtime, " s\n");

            utility::print_par(world, "\tStarting Exchange... ");
            auto k0 = tcc_time::now();
            Eai("X,i,a") = Xab("X,a,b") * Coeffs("b,i");
            auto w1 = tcc_time::now();
            K("a,b") = Eai("X,i,a") * Eai("X,i,b");
            auto k1 = tcc_time::now();
            ktime = tcc_time::duration_in_s(k0, k1);
            auto wtime = tcc_time::duration_in_s(k0, w1);
            utility::print_par(world, ktime, " s\n");
            utility::print_par(world, "\t\tW time ", wtime, " s\n");


            F("i,j") = H("i,j") + 2 * J("i,j") - K("i,j");
            F_TA = TA::to_new_tile_type(F, to_ta);
            // to make energy variational need to use same density for F and energy!
            energy = D_TA("i,j").dot(F_TA("i,j") + H_TA("i,j"), world).get();

            Ferror("i,j") = F_TA("i,k") * D_TA("k,l") * S_TA("l,j")
                            - S_TA("i,k") * D_TA("k,l") * F_TA("l,j");

            error = Ferror("i,j").norm().get() / volume;
            diis.extrapolate(F_TA, Ferror);

            Coeffs_TA = scf::Coeffs_from_fock(F_TA, S_TA, tr_i, n_occ);
            Coeffs = TA::to_new_tile_type(Coeffs_TA, to_decomp);
            D_TA("i,j") = Coeffs_TA("i,a") * Coeffs_TA("j,a");

            delta_e = energy - old_e;
            old_e = energy;

            auto t1 = tcc_time::now();
            time = tcc_time::duration_in_s(t0, t1);

            utility::print_par(world, "\tHas energy ", std::setprecision(17),
                               energy + repulsion_energy, " with error ", error,
                               " in ", time, " s \n");
            utility::print_par(world, "\t\tDelta e = ", delta_e, "\n");
            if (iter == 5) {
                utility::print_par(world, "\n\nTemporary Info:\n");
                utility::print_size_info(Eai, "Exch Temp");
                utility::print_par(world, "\n\n");
            }
            ++iter;
        }

        utility::print_par(world, "\nFinal energy = ", std::setprecision(17),
                           energy + repulsion_energy, "\n");

    // end of SCF
    // prepare CC
        utility::print_par(world, "\nCC Test\n");

      S_TA = TA::to_new_tile_type(S, to_ta);
      F_TA = TA::to_new_tile_type(F, to_ta);
      auto X_ab_TA = TA::to_new_tile_type(Xab, to_ta);

      auto F_eig = array_ops::array_to_eigen(F_TA);
      auto S_eig = array_ops::array_to_eigen(S_TA);
      Eig::GeneralizedSelfAdjointEigenSolver<decltype(S_eig)> es(F_eig,
                                                                 S_eig);
      ens = es.eigenvalues();
      auto C_all = es.eigenvectors();
      decltype(S_eig) C_occ = C_all.leftCols(occupation / 2);
      decltype(S_eig) C_vir = C_all.rightCols(S_eig.rows() - occupation / 2);

      std::size_t all = S.trange().elements().extent()[0];
      tre = std::make_shared<TRange1Engine>(occupation / 2, all, blocksize);

      // start mp2
//            MP2<TA::Tensor<double>, TA::SparsePolicy> mp2(F_TA, S_TA, X_ab_TA, *tre);
//
//            auto two_e = mp2.get_g();
//            mp2.compute();

    auto tr_0 = Xab.trange().data().back();
    auto tr_all = tre->get_all_tr1();
    auto tr_i0 = tre->get_occ_tr1();
    auto tr_vir = tre->get_vir_tr1();

    if (world.rank() ==0){
      std::cout << "TiledRange1 All   ";
      std::cout << tr_all << std::endl;
      std::cout << "TiledRange1 Vir   ";
      std::cout << tr_vir << std::endl;
    }

    auto Ci = array_ops::eigen_to_array<TA::Tensor<double>>(world,C_occ,tr_0, tr_i0);
    auto Cv = array_ops::eigen_to_array<TA::Tensor<double>>(world,C_vir,tr_0, tr_vir);
    auto Call = array_ops::eigen_to_array<TA::Tensor<double>>(world,C_all,tr_0, tr_all);
    g = std::make_shared<tcc::cc::TwoElectronIntMO<TA::Tensor<double>, TA::SparsePolicy>>(X_ab_TA,Ci, Cv);

    fock_mo("p,q") = F_TA("mu,nu")*Call("mu,p")*Call("nu,q");
    }

    world.gop.fence();

    utility::print_par(world, "\nBegining CC\n");


    tcc::cc::CCSD<TA::TensorD, TA::SparsePolicy> ccsd(fock_mo, ens, tre, g);

//            ccsd.compute_cc2();
    ccsd.compute_ccsd();


    world.gop.fence();
    libint2::cleanup();
    madness::finalize();
    return 0;
}


