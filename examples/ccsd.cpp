#include <memory>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <chrono>
#include <rapidjson/document.h>
#include <rapidjson/stringbuffer.h>
#include <rapidjson/writer.h>

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

#include "../integrals/integrals.h"

#include "../scf/diagonalize_for_coffs.hpp"
#include "../cc/ccsd_t.h"
#include "../cc/lazy_tile.h"
#include "../cc/ccsd_intermediates.h"
#include "../cc/trange1_engine.h"
#include "../ta_routines/array_to_eigen.h"
#include "../scf/soad.h"

using namespace mpqc;
namespace ints = integrals;

class ThreeCenterScf {
private:
    using array_type = TA::TSpArrayD;
    array_type H_;
    array_type S_;

    array_type F_;
    array_type D_;
    array_type C_;
    array_type L_invV_;
    TiledArray::DIIS<array_type> diis_;

    std::vector<double> k_times_;
    std::vector<double> j_times_;
    std::vector<double> w_times_;
    std::vector<double> scf_times_;

    int64_t occ_;
    double repulsion_;


    void compute_density(int64_t occ) {
        auto F_eig = tcc::array_ops::array_to_eigen(F_);
        auto S_eig = tcc::array_ops::array_to_eigen(S_);

        Eig::GeneralizedSelfAdjointEigenSolver<decltype(S_eig)> es(F_eig,
                                                                   S_eig);
        decltype(S_eig) C = es.eigenvectors().leftCols(occ);
        auto tr_ao = S_.trange().data()[0];

        auto occ_nclusters = (occ_ < 10) ? occ_ : 10;
        auto tr_occ = tcc::scf::tr_occupied(occ_nclusters, occ_);

        C_ = tcc::array_ops::eigen_to_array<TA::TensorD>(H_.get_world(), C, tr_ao, tr_occ);

        D_("i,j") = C_("i,k") * C_("j,k");
    }

    template <typename Integral>
    void form_fock(Integral const &eri3) {
        auto &world = F_.get_world();

        world.gop.fence();
        auto w0 = tcc::utility::time::now();
        TA::Array<double, 3, TA::TensorD, TA::SparsePolicy> W;
        W("X, mu, i") = L_invV_("X,Y") * (eri3("Y, mu, nu") * C_("nu, i"));
        world.gop.fence();
        auto w1 = tcc::utility::time::now();
        w_times_.push_back(tcc::utility::time::duration_in_s(w0,w1));


        array_type J;
        J("mu, nu") = eri3("X, mu, nu")
                      * (L_invV_("Y, X") * (W("Y, rho, i") * C_("rho, i")));
        world.gop.fence();
        auto j1 = tcc::utility::time::now();
        j_times_.push_back(tcc::utility::time::duration_in_s(w1,j1));


        // Permute W
        W("X,i,nu") = W("X,nu,i");
        array_type K;
        K("mu, nu") = W("X, i, mu") * W("X, i, nu");
        world.gop.fence();
        auto k1 = tcc::utility::time::now();
        k_times_.push_back(tcc::utility::time::duration_in_s(j1,k1));

        F_("i,j") = H_("i,j") + 2 * J("i,j") - K("i,j");
    }


public:

    ThreeCenterScf(array_type const &H, array_type const &F_guess, array_type const &S,
                   array_type const &L_invV, int64_t occ, double rep)
            : H_(H), F_(F_guess), S_(S), L_invV_(L_invV), occ_(occ), repulsion_(rep) {
        compute_density(occ_);
    }


    const array_type get_overlap() const {
        return S_;
    }

    const array_type get_fock() const {
        return F_;
    }

    template <typename Integral>
    void solve(int64_t max_iters, double thresh, Integral const &eri3) {

        if(F_.get_world().rank() == 0){
            std::cout << "Start SCF" << std::endl;
            std::cout << "Convergence : " << thresh << std::endl;
            std::cout << "Max Iteration : " << max_iters << std::endl;
        }

        auto iter = 0;
        auto error = std::numeric_limits<double>::max();
        auto old_energy = 0.0;

        while (iter < max_iters && thresh < error) {
            auto s0 = tcc_time::now();
            F_.get_world().gop.fence();
            form_fock(eri3);

            auto current_energy = energy();
            error = std::abs(old_energy - current_energy);
            old_energy = current_energy;

            array_type Grad;
            Grad("i,j") = F_("i,k") * D_("k,l") * S_("l,j")
                          - S_("i,k") * D_("k,l") * F_("l,j");

            diis_.extrapolate(F_, Grad);

            // Lastly update density
            compute_density(occ_);

            F_.get_world().gop.fence();
            auto s1 = tcc_time::now();
            scf_times_.push_back(tcc_time::duration_in_s(s0, s1));

            if(F_.get_world().rank() == 0){
                std::cout << "Iteration: " << (iter + 1)
                << " energy: " << old_energy << " error: " << error
                << std::endl;
                std::cout << "\tW time: " << w_times_.back() << std::endl;
                std::cout << "\tJ time: " << j_times_.back()
                << " s K time: " << k_times_.back()
                << " s iter time: " << scf_times_.back() << std::endl;
            }
            ++iter;
        }
    }

    double energy() {
        return repulsion_
               + D_("i,j").dot(F_("i,j") + H_("i,j"), D_.get_world()).get();
    }
};

// TODO test case that verify the result automatic
int try_main(int argc, char *argv[], madness::World &world) {


    // parse the input
    rapidjson::Document in;
    parse_input(argc, argv, in);

    std::cout << std::setprecision(15);
    Document cc_in;
    if (in.HasMember("CCSD")){
        cc_in = get_nested(in,"CCSD");
    }
    else if(in.HasMember("CCSD(T)")){
        cc_in = get_nested(in, "CCSD(T)");
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
//    std::shared_ptr<mpqc::cc::CCSDIntermediate<TA::TensorD,TA::SparsePolicy>> intermidiate;
    std::shared_ptr<mpqc::cc::CCSDIntermediate<TA::TensorD,TA::SparsePolicy,cc::DirectTwoElectronSparseArray>> intermidiate;

    std::shared_ptr<mpqc::TRange1Engine> tre;

    Eigen::MatrixXd ens;

    TA::Array<double, 2, TA::TensorD, TA::SparsePolicy> fock_mo;

    cc::DirectTwoElectronSparseArray lazy_two_electron_int;

//    mpqc::integrals::DirArray<4, integrals::IntegralBuilder<4,libint2::TwoBodyEngine<libint2::Coulomb>,integrals::TensorPassThrough>> lazy_two_electron_int;
    {

        // Get necessary info
        std::string mol_file = in["xyz file"].GetString();

        std::string ghost_atoms = in.HasMember("GhostAtoms") ? in["GhostAtoms"].GetString() : "";


        int nclusters = in["number of clusters"].GetInt();
        std::size_t mo_blocksize = cc_in["BlockSize"].GetInt();
        std::size_t ao_blocksize = in.HasMember("AOBlockSize") ? in["AOBlockSize"].GetInt(): mo_blocksize;

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
        double scf_converge = in.HasMember("SCFConverge") ? in["SCFConverge"].GetDouble() : 1.0e-7;


        // get SCF max iteration
        int scf_max_iter = in.HasMember("SCFMaxIter") ? in["SCFMaxIter"].GetInt() : 30;

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
            tcc::utility::print_file(world, mol_file);
            std::cout << "basis is " << basis_name << std::endl;
            std::cout << "df basis is " << df_basis_name << std::endl;
            std::cout << "Using " << nclusters << " clusters"
                      << std::endl;
        }

        tcc::utility::print_par(world, "Sparse threshold is ",
                           TiledArray::SparseShape<float>::threshold(), "\n");

        tcc::utility::print_par(world, "Nuclear repulsion_energy = ",
                           repulsion_energy, "\n");


        world.gop.fence();

        // use clustered_mol to generate basis
        molecule::Molecule clustered_mol{};
        if(!ghost_atoms.empty()){
            auto ghost_molecue = mpqc::molecule::read_xyz(ghost_atoms);
            auto ghost_elements = ghost_molecue.clusterables();

            if (world.rank() == 0){
                std::cout << "Ghost Atom file: " << ghost_atoms << std::endl;
                tcc::utility::print_file(world,ghost_atoms);
            }

            auto mol_elements = mol.clusterables();

            mol_elements.insert(mol_elements.end(),ghost_elements.begin(), ghost_elements.end());

            clustered_mol = mpqc::molecule::kmeans(mol_elements, nclusters);

            clustered_mol.set_charge(mol.charge());
            clustered_mol.set_mass(mol.mass());

        }
        else{

            if (world.rank() == 0){
                std::cout << "Ghost Atom file: None" << std::endl;
            }
            clustered_mol = mpqc::molecule::kmeans(mol.clusterables(), nclusters);
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
        bool if_reblock = in.HasMember("Reblock") ? in["Reblock"].GetBool() : false;
        if(if_reblock){

            tcc::utility::print_par(world,"AOBlockSize:  ",ao_blocksize, "\n");
            basis = reblock(basis,cc::reblock_basis,ao_blocksize);
            df_basis = reblock(df_basis,cc::reblock_basis,ao_blocksize);

        }
        if (world.rank() == 0) {
            TA::TiledRange1 bs_range = basis.create_trange1();
            auto minmax_block = cc::minmax_blocksize(bs_range);
            auto average_block = cc::average_blocksize(bs_range);
            std::cout << "Basis trange " << std::endl;
            std::cout << bs_range << std::endl;
            std::cout << "Min and Max block size: " << minmax_block.first << " " << minmax_block.second << std::endl;
            std::cout << "Average: " << average_block << std::endl;
            TA::TiledRange1 dfbs_range = df_basis.create_trange1();
            minmax_block = cc::minmax_blocksize(dfbs_range);
            average_block = cc::average_blocksize(dfbs_range);
            std::cout << "DF Basis trange " << std::endl;
            std::cout << dfbs_range << std::endl;
            std::cout << "Min and Max block size: " << minmax_block.first << " " << minmax_block.second << std::endl;
            std::cout << "Average: " << average_block << std::endl;
        }

        // start SCF
        libint2::init();

        const auto bs_array = tcc::utility::make_array(basis, basis);

        // Overlap ints
        auto overlap_e = ints::make_1body_shr_pool("overlap", basis, mol);
        auto S = ints::sparse_integrals(world, overlap_e, bs_array);

        // Overlap ints
        auto kinetic_e = ints::make_1body_shr_pool("kinetic", basis, mol);
        auto T = ints::sparse_integrals(world, kinetic_e, bs_array);

        auto nuclear_e = ints::make_1body_shr_pool("nuclear", basis, mol);
        auto V = ints::sparse_integrals(world, nuclear_e, bs_array);

        decltype(T) H;
        H("i,j") = T("i,j") + V("i,j");

        const auto dfbs_array = tcc::utility::make_array(df_basis, df_basis);
        auto eri_e = ints::make_2body_shr_pool(df_basis, basis);

        decltype(H) L_inv;
        {
            auto Vmetric = ints::sparse_integrals(world, eri_e, dfbs_array);
            auto V_eig = tcc::array_ops::array_to_eigen(Vmetric);
            MatrixD Leig = Eig::LLT<MatrixD>(V_eig).matrixL();
            MatrixD L_inv_eig = Leig.inverse();

            auto tr_V = Vmetric.trange().data()[0];
            L_inv = tcc::array_ops::eigen_to_array<TA::TensorD>(world, L_inv_eig,
                                                                tr_V, tr_V);
        }

        auto three_c_array = tcc::utility::make_array(df_basis, basis, basis);
        auto eri3 = ints::sparse_integrals(world, eri_e, three_c_array);


        auto soad0 = tcc::utility::time::now();
        auto F_soad
                = scf::fock_from_soad(world, mol, basis, eri_e, H);

        auto soad1 = tcc::utility::time::now();
        auto soad_time = tcc::utility::time::duration_in_s(soad0, soad1);
        tcc::utility::print_par(world, "Soad Time: " , soad_time, "\n");

        ThreeCenterScf scf(H, F_soad, S, L_inv, occ / 2, repulsion_energy);
        scf.solve(scf_max_iter, scf_converge, eri3);

        // end SCF


        // start ccsd prepration

        tcc::utility::print_par(world, "\nCC Calculation\n");

        int n_frozen_core = 0;
        if (frozen_core) {
            n_frozen_core = mol.core_electrons();
            tcc::utility::print_par(world, "Frozen Core: ", n_frozen_core,
                               " electrons", "\n");
            n_frozen_core = n_frozen_core / 2;
        }

        TA::Array<double,2,TA::TensorD,TA::SparsePolicy> F;
        F = scf.get_fock();

        TA::Array<double,3,TA::TensorD, TA::SparsePolicy> Xab;
        Xab("X,a,b") = L_inv("X,Y")*eri3("Y,a,b");


        auto F_eig = tcc::array_ops::array_to_eigen(F);
        auto S_eig = tcc::array_ops::array_to_eigen(S);

        // check the condition number in Overlap
        Eig::SelfAdjointEigenSolver<decltype(S_eig)> S_es(S_eig);
        // eigen value in increasing order
        auto cond = S_es.eigenvalues()(S_es.eigenvalues().size()-1)/S_es.eigenvalues()(0);
        tcc::utility::print_par(world,"Condition Number in Overlap: ", cond, "\n");

        // solve C
        Eig::GeneralizedSelfAdjointEigenSolver<decltype(S_eig)> es(F_eig, S_eig);
        ens = es.eigenvalues().bottomRows(S_eig.rows() - n_frozen_core);

        auto C_all = es.eigenvectors();
        decltype(S_eig) C_occ = C_all.block(0, n_frozen_core, S_eig.rows(), occ / 2 - n_frozen_core);
        decltype(S_eig) C_vir = C_all.rightCols(S_eig.rows() - occ / 2);
        C_all = C_all.rightCols(S_eig.rows() - n_frozen_core);

        std::size_t all = S.trange().elements().extent()[0];


        // check block size
        std::size_t occ_blocksize = cc_in.HasMember("OccBlockSize") ? cc_in["OccBlockSize"].GetInt() : mo_blocksize;
        std::size_t vir_blocksize = cc_in.HasMember("VirBlockSize") ? cc_in["VirBlockSize"].GetInt() : mo_blocksize;

        tre = std::make_shared<TRange1Engine>(occ / 2, all, occ_blocksize, vir_blocksize, n_frozen_core);

        auto tr_0 = Xab.trange().data().back();
        auto tr_all = tre->get_all_tr1();
        auto tr_i0 = tre->get_occ_tr1();
        auto tr_vir = tre->get_vir_tr1();

        tcc::utility::print_par(world, "Block Size in Occupied     ", occ_blocksize, "\n");
        tcc::utility::print_par(world, "TiledRange1 Occupied ", tr_i0, "\n");
        tcc::utility::print_par(world, "Average: ", cc::average_blocksize(tr_i0), "\n");
        auto min_max = cc::minmax_blocksize(tr_i0);
        tcc::utility::print_par(world, "Min and Max block size: ",min_max.first, " ", min_max.second, "\n");


        tcc::utility::print_par(world, "Block Size in Virtual     ", vir_blocksize, "\n");
        tcc::utility::print_par(world, "TiledRange1 Virtual  ", tr_vir, "\n");
        tcc::utility::print_par(world, "Average: ", cc::average_blocksize(tr_vir), "\n");
        min_max = cc::minmax_blocksize(tr_vir);
        tcc::utility::print_par(world, "Min and Max block size: ",min_max.first, " ", min_max.second, "\n");

        auto Ci = tcc::array_ops::eigen_to_array<TA::Tensor<double>>(world, C_occ, tr_0, tr_i0);

        auto Cv = tcc::array_ops::eigen_to_array<TA::Tensor<double>>(world, C_vir, tr_0, tr_vir);

        auto Call = tcc::array_ops::eigen_to_array<TA::Tensor<double>>(world, C_all, tr_0, tr_all);

        std::vector<TA::TiledRange1> tr_04(4, basis.create_trange1());
        TA::TiledRange trange_4(tr_04.begin(), tr_04.end());


        world.gop.fence();

        fock_mo("p,q") = F("mu,nu") * Cv("mu,p") * Ci("nu,q");

//        if (do_screen) {
//            auto screen_builder = ints::init_schwarz_screen(1e-10);
//            auto shr_screen = std::make_shared<ints::SchwarzScreen>(screen_builder(world, eri_e, basis));
//
//            const auto bs4_array = tcc::utility::make_array(basis, basis, basis, basis);
//            auto lazy_two_electron_int = mpqc_ints::direct_sparse_integrals(world, eri_e, bs4_array, shr_screen);
//            intermidiate = std::make_shared<mpqc::cc::CCSDIntermediate<TA::TensorD, TA::SparsePolicy>>
//                    (Xab, Ci, Cv, lazy_two_electron_int);
//        } else {

//            const auto bs4_array = tcc::utility::make_array(basis, basis, basis, basis);
//            auto lazy_two_electron_int = mpqc_ints::direct_sparse_integrals(world, eri_e, bs4_array);

        std::string screen = cc_in.HasMember("Screen") ? cc_in["Screen"].GetString() : "";
        int screen_option = 0;
        if(screen == "schwarz"){
            screen_option = 1;
        }
        else if(screen =="qqr"){
            screen_option = 2;
        }


        auto direct = cc_in.HasMember("Direct") ? cc_in["Direct"].GetBool(): true;

        if(direct){

            auto time0 = tcc_time::now();
            lazy_two_electron_int = cc::make_lazy_two_electron_sparse_array(world, basis, trange_4,screen_option);
            auto time1 = tcc_time::now();
            auto duration = tcc_time::duration_in_s(time0,time1);
            if(world.rank() == 0){
                std::cout << "Time to initialize direct two electron sparse integral: " << duration << std::endl;

            }
        }

        intermidiate = std::make_shared<mpqc::cc::CCSDIntermediate<TA::TensorD, TA::SparsePolicy, cc::DirectTwoElectronSparseArray>>
                (Ci, Cv, Xab, lazy_two_electron_int);
//        }
    }

    // clean up all temporary from HF
    world.gop.fence();

    tcc::utility::print_par(world, "\nBegining CC Calculation\n");
    tcc::utility::parallal_break_point(world, 0);


    if(in.HasMember("CCSD(T)")){
        mpqc::cc::CCSD_T<TA::Tensor < double>, TA::SparsePolicy > ccsd_t(fock_mo, ens, tre, intermidiate, cc_in);
        ccsd_t.compute();
    }
    else if(in.HasMember("CCSD")){
        mpqc::cc::CCSD<TA::Tensor < double>, TA::SparsePolicy > ccsd(fock_mo, ens, tre, intermidiate, cc_in);
        ccsd.compute();
    }

    world.gop.fence();
    libint2::cleanup();
    return 0;
}


int main(int argc, char *argv[]) {

    int rc = 0;

    auto &world = madness::initialize(argc, argv);

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
