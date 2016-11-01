#include <memory>
#include <fstream>
#include <iomanip>

#include <rapidjson/document.h>
#include <rapidjson/stringbuffer.h>
#include <rapidjson/writer.h>
#include <tiledarray.h>

#include "../utility/make_array.h"
#include "../utility/parallel_print.h"
#include "../utility/parallel_break_point.h"

#include "../utility/array_info.h"
#include "../utility/print_size_info.h"

#include "mpqc/util/misc/time.h"
#include "../utility/json_handling.h"

#include "../molecule/atom.h"
#include "../molecule/cluster.h"
#include "../molecule/molecule.h"
#include "../molecule/clustering_functions.h"
#include "../molecule/make_clusters.h"

#include "../basis/atom_basisset.h"
#include "../basis/basis_set.h"
#include "../basis/cluster_shells.h"
#include "../basis/basis.h"

#include "../scf/diagonalize_for_coffs.hpp"
#include "../scf/scf.h"
#include "../scf/traditional_df_fock_builder.h"
#include "../scf/purification_density_build.h"
#include "../scf/eigen_solve_density_builder.h"

#include "../scf/traditional_four_center_fock_builder.h"

#include "../cc/ccsd_t.h"
#include "../cc/lazy_tile.h"
#include "../cc/ccsd_intermediates.h"
#include "../utility/trange1_engine.h"
#include "mpqc/math/external/eigen/eigen.h"
#include "../scf/soad.h"
#include "../eom_cc/eom_ccsd.h"

using namespace mpqc;
namespace ints = integrals;


// TODO test case that verify the result automatic
int try_main(int argc, char *argv[], madness::World &world) {


    // parse the input
    rapidjson::Document in;
    json::parse_input(argc, argv, in);

    std::cout << std::setprecision(15);
    rapidjson::Document cc_in;
    if (in.HasMember("CCSD")) {
        cc_in = json::get_nested(in, "CCSD");
    }

    if (!in.HasMember("xyz file") || !in.HasMember("number of clusters")) {
        if (world.rank() == 0) {
            std::cout << "At a minimum your input file must provide\n";
            std::cout << "\"xyz file\", which is path to an xyz input\n";
            std::cout << "\"number of clusters\", which is the number of "
                         "clusters in the obs\n";
        }
    }

    // declare variables needed for ccsd
    //    std::shared_ptr<mpqc::cc::CCSDIntermediate<TA::TensorD,TA::SparsePolicy>>
    //    intermidiate;
    std::shared_ptr<mpqc::cc::
                          CCSDIntermediate<TA::TensorD, TA::SparsePolicy,
                                           cc::DirectTwoElectronSparseArray>>
          intermidiate;

    std::shared_ptr<mpqc::TRange1Engine> tre;

    Eigen::MatrixXd ens;

    TA::Array<double, 2, TA::TensorD, TA::SparsePolicy> fock_mo;

    cc::DirectTwoElectronSparseArray lazy_two_electron_int;

    //    mpqc::integrals::DirArray<4,
    //    integrals::IntegralBuilder<4,libint2::TwoBodyEngine<libint2::Coulomb>,integrals::TensorPassThrough>>
    //    lazy_two_electron_int;
    {

        // Get necessary info
        std::string mol_file = in["xyz file"].GetString();

        std::string ghost_atoms
              = in.HasMember("GhostAtoms") ? in["GhostAtoms"].GetString() : "";


        int nclusters = in["number of clusters"].GetInt();
        std::size_t mo_blocksize = cc_in["BlockSize"].GetInt();
        std::size_t ao_blocksize = in.HasMember("AOBlockSize")
                                         ? in["AOBlockSize"].GetInt()
                                         : mo_blocksize;

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

        // get SCF converge
        double scf_converge = in.HasMember("SCFConverge")
                                    ? in["SCFConverge"].GetDouble()
                                    : 1.0e-7;


        // get SCF max iteration
        int scf_max_iter
              = in.HasMember("SCFMaxIter") ? in["SCFMaxIter"].GetInt() : 30;

        // get other info
        bool frozen_core = cc_in.HasMember("FrozenCore")
                                 ? cc_in["FrozenCore"].GetBool()
                                 : false;

        TiledArray::SparseShape<float>::threshold(threshold);

        auto mol = mpqc::molecule::read_xyz(mol_file);
        auto charge = 0;
        auto occ = mol.occupation(charge);
        auto repulsion_energy = mol.nuclear_repulsion();


        if (world.rank() == 0) {
            std::cout << "Mol file is " << mol_file << std::endl;
            utility::print_file(world, mol_file);
            std::cout << "basis is " << basis_name << std::endl;
            std::cout << "df basis is " << df_basis_name << std::endl;
            std::cout << "Using " << nclusters << " clusters" << std::endl;
        }

        utility::print_par(world, "Sparse threshold is ",
                           TiledArray::SparseShape<float>::threshold(), "\n");

        utility::print_par(world, "Nuclear repulsion_energy = ",
                           repulsion_energy, "\n");


        world.gop.fence();

        // use clustered_mol to generate basis
        molecule::Molecule clustered_mol{};
        if (!ghost_atoms.empty()) {
            auto ghost_molecue = mpqc::molecule::read_xyz(ghost_atoms);
            auto ghost_elements = ghost_molecue.clusterables();

            if (world.rank() == 0) {
                std::cout << "Ghost Atom file: " << ghost_atoms << std::endl;
                utility::print_file(world, ghost_atoms);
            }

            auto mol_elements = mol.clusterables();

            mol_elements.insert(mol_elements.end(), ghost_elements.begin(),
                                ghost_elements.end());

            clustered_mol = mpqc::molecule::kmeans(mol_elements, nclusters);

            clustered_mol.set_charge(mol.charge());
            clustered_mol.set_mass(mol.mass());

        } else {

            if (world.rank() == 0) {
                std::cout << "Ghost Atom file: None" << std::endl;
            }
            clustered_mol
                  = mpqc::molecule::kmeans(mol.clusterables(), nclusters);
        }

        mpqc::basis::BasisSet bs{basis_name};
        mpqc::basis::BasisSet df_bs{df_basis_name};

        std::streambuf *cout_sbuf
              = std::cout.rdbuf(); // Silence libint printing.
        std::ofstream fout("/dev/null");
        std::cout.rdbuf(fout.rdbuf());
        mpqc::basis::Basis basis{bs.get_cluster_shells(clustered_mol)};
        mpqc::basis::Basis df_basis{df_bs.get_cluster_shells(clustered_mol)};
        std::cout.rdbuf(cout_sbuf);

        // reblock basis
        bool if_reblock = in.HasMember("Reblock") ? in["Reblock"].GetBool()
                                                  : false;
        if (if_reblock) {

            utility::print_par(world, "AOBlockSize:  ", ao_blocksize, "\n");
            basis = reblock(basis, cc::reblock_basis, ao_blocksize);
            df_basis = reblock(df_basis, cc::reblock_basis, ao_blocksize);
        }
        if (world.rank() == 0) {
            TA::TiledRange1 bs_range = basis.create_trange1();
            auto minmax_block = cc::minmax_blocksize(bs_range);
            auto average_block = cc::average_blocksize(bs_range);
            std::cout << "Basis trange " << std::endl;
            std::cout << bs_range << std::endl;
            std::cout << "Min and Max block size: " << minmax_block.first << " "
                      << minmax_block.second << std::endl;
            std::cout << "Average: " << average_block << std::endl;
            TA::TiledRange1 dfbs_range = df_basis.create_trange1();
            minmax_block = cc::minmax_blocksize(dfbs_range);
            average_block = cc::average_blocksize(dfbs_range);
            std::cout << "DF Basis trange " << std::endl;
            std::cout << dfbs_range << std::endl;
            std::cout << "Min and Max block size: " << minmax_block.first << " "
                      << minmax_block.second << std::endl;
            std::cout << "Average: " << average_block << std::endl;
        }

        // start SCF
        libint2::init();

        const auto bs_array = utility::make_array(basis, basis);

        // Overlap ints
        auto time0 = mpqc_time::fenced_now(world);
        auto overlap_e = ints::make_1body_shr_pool("overlap", basis, mol);
        auto S = ints::sparse_integrals(world, overlap_e, bs_array);
        auto time1 = mpqc_time::fenced_now(world);
        auto time = mpqc_time::duration_in_s(time0, time1);
        mpqc::utility::print_par(world, "Overlap Time:  ", time, "\n");

        // Kinetic ints
        time0 = mpqc_time::fenced_now(world);
        auto kinetic_e = ints::make_1body_shr_pool("kinetic", basis, mol);
        auto T = ints::sparse_integrals(world, kinetic_e, bs_array);
        time1 = mpqc_time::fenced_now(world);
        time = mpqc_time::duration_in_s(time0, time1);
        mpqc::utility::print_par(world, "Kinetic Time:  ", time, "\n");

        time0 = mpqc_time::fenced_now(world);
        auto nuclear_e = ints::make_1body_shr_pool("nuclear", basis, mol);
        auto V = ints::sparse_integrals(world, nuclear_e, bs_array);
        time1 = mpqc_time::fenced_now(world);
        time = mpqc_time::duration_in_s(time0, time1);
        mpqc::utility::print_par(world, "Nuclear Time:  ", time, "\n");

        time0 = mpqc_time::fenced_now(world);
        decltype(T) H;
        H("i,j") = T("i,j") + V("i,j");
        time1 = mpqc_time::fenced_now(world);
        time = mpqc_time::duration_in_s(time0, time1);
        mpqc::utility::print_par(world, "Core Time:  ", time, "\n");

        time0 = mpqc_time::fenced_now(world);
        auto eri_e = ints::make_2body_shr_pool(df_basis, basis);
        auto F_soad
              = scf::fock_from_soad(world, clustered_mol, basis, eri_e, H);
        time1 = mpqc_time::fenced_now(world);
        time = mpqc_time::duration_in_s(time0, time1);
        mpqc::utility::print_par(world, "Soad Time:  ", time, "\n");

        time0 = mpqc_time::fenced_now(world);
        auto three_c_array = utility::make_array(df_basis, basis, basis);
        auto eri3 = ints::sparse_integrals(world, eri_e, three_c_array);
        time1 = mpqc_time::fenced_now(world);
        time = mpqc_time::duration_in_s(time0, time1);
        mpqc::utility::print_par(world, "Three Center Time:  ", time, "\n");

        time0 = mpqc_time::fenced_now(world);
        const auto dfbs_array = utility::make_array(df_basis, df_basis);
        auto Metric = ints::sparse_integrals(world, eri_e, dfbs_array);
        scf::DFFockBuilder<decltype(eri3)> builder(Metric, eri3);
        time1 = mpqc_time::fenced_now(world);
        time = mpqc_time::duration_in_s(time0, time1);
        mpqc::utility::print_par(world, "Two Center Time:  ", time, "\n");

        std::unique_ptr<scf::FockBuilder> f_builder;
        if (in.HasMember("Fock Builder")
            && in["Fock Builder"].GetString() == std::string("four center")) {
            auto four_c_array = utility::make_array(basis, basis, basis, basis);

            time0 = mpqc_time::fenced_now(world);
            auto eri4 = ints::sparse_integrals(world, eri_e, four_c_array);
            time1 = mpqc_time::fenced_now(world);
            time = mpqc_time::duration_in_s(time0, time1);
            mpqc::utility::print_par(world, "Four Center Time: ", time, "\n");

            auto builder
                  = scf::FourCenterBuilder<decltype(eri4)>(std::move(eri4));

            f_builder = make_unique<decltype(builder)>(std::move(builder));
        } else {


            f_builder = make_unique<decltype(builder)>(std::move(builder));
        }

        time0 = mpqc_time::fenced_now(world);
        auto multi_pool
              = ints::make_1body_shr_pool("emultipole2", basis, clustered_mol);
        auto r_xyz = ints::sparse_xyz_integrals(world, multi_pool, bs_array);
        time1 = mpqc_time::fenced_now(world);
        time = mpqc_time::duration_in_s(time0, time1);
        mpqc::utility::print_par(world, "Multipole Integral Time:  ", time,
                                 "\n");

        std::unique_ptr<scf::DensityBuilder> d_builder;
        if (in.HasMember("Density Builder")
            && in["DensityBuilder"].GetString()
               == std::string("purification")) {
            auto db = scf::PurificationDensityBuilder(S, r_xyz, occ / 2,
                                                      nclusters, 0.0, false);

            d_builder
                  = make_unique<scf::PurificationDensityBuilder>(std::move(db));
        } else {
            auto db = scf::ESolveDensityBuilder(S, r_xyz, occ / 2, nclusters,
                                                0.0, "cholesky inverse", false);
            d_builder = make_unique<scf::ESolveDensityBuilder>(std::move(db));
        }

        scf::ClosedShellSCF scf(H, S, repulsion_energy, std::move(f_builder),
                                std::move(d_builder), F_soad);

        scf.solve(scf_max_iter, scf_converge);
        // end SCF

        // start ccsd prepration

        utility::print_par(world, "\nCC Calculation\n");

        int n_frozen_core = 0;
        if (frozen_core) {
            n_frozen_core = mol.core_electrons();
            utility::print_par(world, "Frozen Core: ", n_frozen_core,
                               " electrons", "\n");
            n_frozen_core = n_frozen_core / 2;
        }

        TA::Array<double, 2, TA::TensorD, TA::SparsePolicy> F;
        F = scf.fock();

        decltype(Metric) L_inv;
        {
            auto M_eig = array_ops::array_to_eigen(Metric);
            MatrixD L_inv_eig
                  = MatrixD(Eig::LLT<MatrixD>(M_eig).matrixL()).inverse();
            auto tr_M = Metric.trange().data()[0];
            L_inv = array_ops::eigen_to_array<TA::TensorD>(world, L_inv_eig,
                                                           tr_M, tr_M);
        }

        TA::Array<double, 3, TA::TensorD, TA::SparsePolicy> Xab;
        Xab("X,a,b") = L_inv("X,Y") * eri3("Y,a,b");

        auto F_eig = array_ops::array_to_eigen(F);
        auto S_eig = array_ops::array_to_eigen(S);

        // check the condition number in Overlap
        Eig::SelfAdjointEigenSolver<decltype(S_eig)> S_es(S_eig);
        // eigen value in increasing order
        auto cond = S_es.eigenvalues()(S_es.eigenvalues().size() - 1)
                    / S_es.eigenvalues()(0);
        utility::print_par(world, "Condition Number in Overlap: ", cond, "\n");

        // solve C
        Eig::GeneralizedSelfAdjointEigenSolver<decltype(S_eig)> es(F_eig,
                                                                   S_eig);
        ens = es.eigenvalues().bottomRows(S_eig.rows() - n_frozen_core);

        auto C_all = es.eigenvectors();
        decltype(S_eig) C_occ = C_all.block(0, n_frozen_core, S_eig.rows(),
                                            occ / 2 - n_frozen_core);
        decltype(S_eig) C_vir = C_all.rightCols(S_eig.rows() - occ / 2);
        C_all = C_all.rightCols(S_eig.rows() - n_frozen_core);

        std::size_t all = S.trange().elements_range().extent()[0];


        // check block size
        std::size_t occ_blocksize = cc_in.HasMember("OccBlockSize")
                                          ? cc_in["OccBlockSize"].GetInt()
                                          : mo_blocksize;
        std::size_t vir_blocksize = cc_in.HasMember("VirBlockSize")
                                          ? cc_in["VirBlockSize"].GetInt()
                                          : mo_blocksize;

        tre = std::make_shared<TRange1Engine>(occ / 2, all, occ_blocksize,
                                              vir_blocksize, n_frozen_core);

        auto tr_0 = Xab.trange().data().back();
        auto tr_all = tre->get_all_tr1();
        auto tr_i0 = tre->get_occ_tr1();
        auto tr_vir = tre->get_vir_tr1();

        utility::print_par(world, "Block Size in Occupied     ", occ_blocksize,
                           "\n");
        utility::print_par(world, "TiledRange1 Occupied ", tr_i0, "\n");
        utility::print_par(world, "Average: ", cc::average_blocksize(tr_i0),
                           "\n");
        auto min_max = cc::minmax_blocksize(tr_i0);
        utility::print_par(world, "Min and Max block size: ", min_max.first,
                           " ", min_max.second, "\n");


        utility::print_par(world, "Block Size in Virtual     ", vir_blocksize,
                           "\n");
        utility::print_par(world, "TiledRange1 Virtual  ", tr_vir, "\n");
        utility::print_par(world, "Average: ", cc::average_blocksize(tr_vir),
                           "\n");
        min_max = cc::minmax_blocksize(tr_vir);
        utility::print_par(world, "Min and Max block size: ", min_max.first,
                           " ", min_max.second, "\n");

        auto Ci = array_ops::eigen_to_array<TA::Tensor<double>>(world, C_occ,
                                                                tr_0, tr_i0);

        auto Cv = array_ops::eigen_to_array<TA::Tensor<double>>(world, C_vir,
                                                                tr_0, tr_vir);

        auto Call = array_ops::eigen_to_array<TA::Tensor<double>>(world, C_all,
                                                                  tr_0, tr_all);

        std::vector<TA::TiledRange1> tr_04(4, basis.create_trange1());
        TA::TiledRange trange_4(tr_04.begin(), tr_04.end());


        world.gop.fence();

        fock_mo("p,q") = F("mu,nu") * Cv("mu,p") * Ci("nu,q");

        //        if (do_screen) {
        //            auto screen_builder =
        //            ints::init_schwarz_screen(1e-10);
        //            auto shr_screen =
        //            std::make_shared<ints::SchwarzScreen>(screen_builder(world,
        //            eri_e, basis));
        //
        //            const auto bs4_array = utility::make_array(basis,
        //            basis,
        //            basis, basis);
        //            auto lazy_two_electron_int =
        //            mpqc_ints::direct_sparse_integrals(world, eri_e,
        //            bs4_array, shr_screen);
        //            intermidiate =
        //            std::make_shared<mpqc::cc::CCSDIntermediate<TA::TensorD,
        //            TA::SparsePolicy>>
        //                    (Xab, Ci, Cv, lazy_two_electron_int);
        //        } else {

        //            const auto bs4_array = utility::make_array(basis,
        //            basis,
        //            basis, basis);
        //            auto lazy_two_electron_int =
        //            mpqc_ints::direct_sparse_integrals(world, eri_e,
        //            bs4_array);

        std::string screen
              = cc_in.HasMember("Screen") ? cc_in["Screen"].GetString() : "";
        int screen_option = 0;
        if (screen == "schwarz") {
            screen_option = 1;
        } else if (screen == "qqr") {
            screen_option = 2;
        }


        auto direct = cc_in.HasMember("Direct") ? cc_in["Direct"].GetBool()
                                                : true;

        if (direct) {

            auto time0 = mpqc_time::now();
            lazy_two_electron_int = cc::make_lazy_two_electron_sparse_array(
                  world, basis, trange_4, screen_option);
            auto time1 = mpqc_time::now();
            auto duration = mpqc_time::duration_in_s(time0, time1);
            if (world.rank() == 0) {
                std::cout << "Time to initialize direct two electron sparse "
                             "integral: " << duration << std::endl;
            }
        }

        intermidiate = std::
              make_shared<mpqc::cc::
                                CCSDIntermediate<TA::TensorD, TA::SparsePolicy,
                                                 cc::DirectTwoElectronSparseArray>>(
                    Ci, Cv, Xab, lazy_two_electron_int);
        //        }
    }

    // clean up all temporary from HF
    world.gop.fence();

    utility::print_par(world, "\nBegining CC Calculation\n");
    utility::parallal_break_point(world, 0);


    if (in.HasMember("EOM-CCSD")) {
      // call ccsd
      mpqc::cc::CCSD<TA::Tensor<double>, TA::SparsePolicy> ccsd(
            fock_mo, ens, tre, intermidiate, cc_in);
      ccsd.compute();

      // call eom-cc
      rapidjson::Document eomccsd_in;
      eomccsd_in = json::get_nested(in, "EOM-CCSD");

      mpqc::EOM_CCSD eomccsd(ccsd, intermidiate);
      eomccsd.read_guess_vectors(eomccsd_in);

      std:size_t MaxIter = eomccsd_in.HasMember("MaxIter")
                          ? eomccsd_in["MaxIter"].GetInt()
                          : 50;
      if (world.rank() == 0) {
        std::cout << "MaxIter: " << MaxIter << std::endl;
      }
      double Convergence = eomccsd_in.HasMember("Convergence")
                                  ? eomccsd_in["Convergence"].GetDouble()
                                  : 1.0e-6;

      if (world.rank() == 0) {
        std::cout << "Convergence: " << Convergence << std::endl;
      }
      //eomccsd.compute_energy(MaxIter, Convergence);



    }

    world.gop.fence();
    libint2::cleanup();
    return 0;
}


int main(int argc, char *argv[]) {

    int rc = 0;

    auto &world = madness::initialize(argc, argv);
    mpqc::utility::print_par(world, "MADNESS process total size: ",
                             world.size(), "\n");

    try {

        try_main(argc, argv, world);

    } catch (TiledArray::Exception &e) {
        std::cerr << "!! TiledArray exception: " << e.what() << "\n";
        rc = 1;
    } catch (madness::MadnessException &e) {
        std::cerr << "!! MADNESS exception: " << e.what() << "\n";
        rc = 1;
    } catch (SafeMPI::Exception &e) {
        std::cerr << "!! SafeMPI exception: " << e.what() << "\n";
        rc = 1;
    } catch (std::exception &e) {
        std::cerr << "!! std exception: " << e.what() << "\n";
        rc = 1;
    } catch (...) {
        std::cerr << "!! exception: unknown exception\n";
        rc = 1;
    }


    madness::finalize();
    return rc;
}
