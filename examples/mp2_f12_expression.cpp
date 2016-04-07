#include "../common/namespaces.h"
#include "../common/typedefs.h"

#include "../include/libint.h"
#include "../include/tiledarray.h"

#include "../utility/make_array.h"
#include "../clustering/kmeans.h"

#include "../molecule/atom.h"
#include "../molecule/atom_based_cluster.h"
#include "../molecule/molecule.h"
#include "../molecule/clustering_functions.h"
#include "../molecule/make_clusters.h"

#include "../scf/scf.h"

#include "../basis/atom_basisset.h"
#include "../basis/basis_set.h"
#include "../basis/basis.h"

#include "../f12/mp2f12.h"
#include "../utility/cc_utility.h"
#include "../integrals/integrals.h"
#include "../integrals/atomic_integral.h"
#include "../expression/orbital_registry.h"

#include "../utility/time.h"
#include "../utility/wcout_utf8.h"
#include "../utility/array_info.h"
#include "../ta_routines/array_to_eigen.h"
#include "../scf/traditional_df_fock_builder.h"
#include "../scf/eigen_solve_density_builder.h"

#include "../mp2/mp2.h"
#include "../scf/traditional_four_center_fock_builder.h"

#include <memory>
#include <locale>

using namespace mpqc;
namespace ints = mpqc::integrals;


TA::TensorD ta_pass_through(TA::TensorD &&ten) { return std::move(ten); }

int main(int argc, char *argv[]) {
    auto &world = madness::initialize(argc, argv);
    std::wcout.imbue(std::locale("en_US.UTF-8"));
    std::string mol_file = "";
    std::string basis_name = "";
    std::string aux_basis_name = "";
    std::string df_basis_name = "";
    int nclusters = 0;
    std::cout << std::setprecision(15);
    double threshold = 1e-30;
    double well_sep_threshold = 0.1;
    integrals::QQR::well_sep_threshold(well_sep_threshold);
    if (argc == 6) {
        mol_file = argv[1];
        basis_name = argv[2];
        aux_basis_name = argv[3];
        df_basis_name = argv[4];
        nclusters = std::stoi(argv[5]);
    } else {
        std::cout << "input is $./program mol_file basis_file aux_basis df_basis nclusters ";
        return 0;
    }
    TiledArray::SparseShape<float>::threshold(threshold);

    auto clustered_mol = molecule::kmeans(
            molecule::read_xyz(mol_file).clusterables(), nclusters);

    std::cout << "Molecule " << std::endl;

    auto atoms = clustered_mol.atoms();

    for (auto& atom: atoms){
        std::cout << atom << std::endl;
    }

    std::cout << std::endl;

    auto repulsion_energy = clustered_mol.nuclear_repulsion();
    std::cout << "Nuclear Repulsion Energy: " << repulsion_energy << std::endl;
    auto occ = clustered_mol.occupation(0);

    basis::BasisSet bs(basis_name);
    basis::Basis basis(bs.get_cluster_shells(clustered_mol));

    basis::BasisSet dfbs(df_basis_name);
    basis::Basis df_basis(dfbs.get_cluster_shells(clustered_mol));

    basis::BasisSet abs(aux_basis_name);
    basis::Basis abs_basis(abs.get_cluster_shells(clustered_mol));

    basis::Basis ri_basis = basis.join(abs_basis);

    auto bs_registry = std::make_shared<OrbitalBasisRegistry>();
    bs_registry->add(OrbitalIndex(L"κ"),std::make_shared<basis::Basis>(basis));
    bs_registry->add(OrbitalIndex(L"Κ"),std::make_shared<basis::Basis>(df_basis));
    bs_registry->add(OrbitalIndex(L"α"), std::make_shared<basis::Basis>(abs_basis));
    bs_registry->add(OrbitalIndex(L"ρ"), std::make_shared<basis::Basis>(ri_basis));


    f12::GTGParams gtg_params(1.0, 6);


    auto param = gtg_params.compute();

    for(auto& pair : param){
        std::cout << pair.first << " " << pair.second << std::endl;
    }

    libint2::init();

    integrals::AtomicIntegral<TA::TensorD, TA::SparsePolicy> ao_int
            (world,
             ta_pass_through,
             std::make_shared<molecule::Molecule>(clustered_mol),
             bs_registry,
             gtg_params.compute()
            );

    // Overlap ints
    auto S = ao_int.compute(L"(κ|λ)");
    auto H = ao_int.compute(L"(κ|H|λ)");

    auto eri4 = ao_int.compute(L"( κ1 λ1 | G|κ1 λ1)");
    scf::FourCenterBuilder<decltype(eri4)> builder(eri4);
    world.gop.fence();

    const auto bs_array = utility::make_array(basis, basis);
    auto multi_pool = ints::make_1body_shr_pool("emultipole2", basis, clustered_mol);
    auto r_xyz = ints::sparse_xyz_integrals(world, multi_pool, bs_array);

    auto db = scf::ESolveDensityBuilder(S, r_xyz, occ / 2, nclusters, 0.0, "cholesky inverse", false);

    scf::ClosedShellSCF scf(H, S, repulsion_energy, std::move(builder), std::move(db));
    scf.solve(50, 1e-12);

    // obs fock build
    std::size_t all = S.trange().elements().extent()[0];
    auto tre = TRange1Engine(occ / 2, all, 1, 4, 0);

    auto F = scf.fock();
    ao_int.registry().insert(Formula(L"(μ|F|ν)"),F);


    //mp2
    // solve Coefficient
    std::size_t n_frozen_core = 0;
    auto F_eig = array_ops::array_to_eigen(F);
    auto S_eig = array_ops::array_to_eigen(S);
    Eig::GeneralizedSelfAdjointEigenSolver<decltype(S_eig)> es(F_eig, S_eig);
    Eig::VectorXd ens = es.eigenvalues().bottomRows(S_eig.rows() - n_frozen_core);

    Eig::MatrixXd C_all = es.eigenvectors();
    decltype(S_eig) C_occ = C_all.block(0, 0, S_eig.rows(),occ / 2);
    decltype(S_eig) C_occ_corr = C_all.block(0, n_frozen_core, S_eig.rows(),occ / 2 - n_frozen_core);
    decltype(S_eig) C_vir = C_all.rightCols(S_eig.rows() - occ / 2);

    auto tr_0 = eri4.trange().data().back();
    auto tr_all = tre.get_all_tr1();
    auto tr_i0 = tre.get_occ_tr1();
    auto tr_vir = tre.get_vir_tr1();


    auto Ci = array_ops::eigen_to_array<TA::Tensor<double>>(world, C_occ_corr, tr_0, tr_i0);
    auto Cv = array_ops::eigen_to_array<TA::Tensor<double>>(world, C_vir, tr_0, tr_vir);
    auto Call = array_ops::eigen_to_array<TA::Tensor<double>>(world, C_all, tr_0, tr_all);

    auto orbital_registry = OrbitalSpaceRegistry<decltype(Ci)>();
    using OrbitalSpace = OrbitalSpace<decltype(Ci)>;
    auto occ_space = OrbitalSpace(OrbitalIndex(L"m"), OrbitalIndex(L"κ"), Ci);
    orbital_registry.add(occ_space);

    auto corr_occ_space = OrbitalSpace(OrbitalIndex(L"i"), OrbitalIndex(L"κ"), Ci);
    orbital_registry.add(corr_occ_space);

    auto vir_space = OrbitalSpace(OrbitalIndex(L"a"), OrbitalIndex(L"κ"), Cv);
    orbital_registry.add(vir_space);

    auto obs_space = OrbitalSpace(OrbitalIndex(L"p"),OrbitalIndex(L"κ"), Call);
    orbital_registry.add(obs_space);

    auto sp_orbital_registry = std::make_shared<decltype(orbital_registry)>(orbital_registry);

    // clean atomic integral
//    ao_int.registry().clear();

    auto mo_integral = integrals::MolecularIntegral<TA::TensorD,TA::SparsePolicy>(ao_int,sp_orbital_registry);
//    mo_integral.atomic_integral().registry().print_formula();


    // test mp2


    // mp2
    {
        auto g_iajb = mo_integral.compute(L"(i a|G|j b)");
        auto mp2 = MP2<TA::TensorD, TA::SparsePolicy>(g_iajb,ens,std::make_shared<TRange1Engine>(tre));
        mp2.compute();
    }

    // df-mp2

    {
        auto g_iajb = mo_integral.compute(L"(i a|G|j b)[df]");
        auto mp2 = MP2<TA::TensorD, TA::SparsePolicy>(g_iajb,ens,std::make_shared<TRange1Engine>(tre));
        mp2.compute();
    }

    // physics notation
    {
        auto g_ijab = mo_integral.compute(L"<i j|G|a b>[df]");
        g_ijab("i,a,j,b") = g_ijab("i,j,a,b");
        auto mp2 = MP2<TA::TensorD, TA::SparsePolicy>(g_ijab,ens,std::make_shared<TRange1Engine>(tre));
        mp2.compute();
    }


    {
        auto g_ijab = mo_integral.compute(L"<i j|G|a b>");
        g_ijab("i,a,j,b") = g_ijab("i,j,a,b");
        auto mp2 = MP2<TA::TensorD, TA::SparsePolicy>(g_ijab,ens,std::make_shared<TRange1Engine>(tre));
        mp2.compute();

    }

//    mo_integral.atomic_integral().registry().print_formula();
    // CABS fock build

    // integral

    auto S_cabs = ao_int.compute(L"(α|β)");
    auto S_ribs = ao_int.compute(L"(ρ|σ)");
    auto S_obs_ribs = ao_int.compute(L"(μ|σ)");
    auto S_obs = scf.overlap();

//    std::cout << S_obs << std::endl;
//    std::cout << S_ribs << std::endl;

    // construct cabs
    decltype(S_obs) C_cabs, C_ri;
    {
        auto S_obs_eigen = array_ops::array_to_eigen(S_obs);
//        MatrixD X_obs_eigen_inv = MatrixD(Eigen::LLT<MatrixD>(S_obs_eigen).matrixL()).inverse();
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(S_obs_eigen);
        MatrixD X_obs_eigen_inv = es.operatorInverseSqrt();

        auto S_ribs_eigen = array_ops::array_to_eigen(S_ribs);
//        MatrixD X_ribs_eigen_inv = MatrixD(Eigen::LLT<MatrixD>(S_ribs_eigen).matrixL()).inverse();
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es2(S_ribs_eigen);
        MatrixD X_ribs_eigen_inv = es2.operatorInverseSqrt();

        MatrixD S_obs_ribs_eigen = array_ops::array_to_eigen(S_obs_ribs);

        MatrixD S_obs_ribs_ortho_eigen = X_obs_eigen_inv.transpose() * S_obs_ribs_eigen * X_ribs_eigen_inv;

        Eigen::JacobiSVD<MatrixD> svd(S_obs_ribs_ortho_eigen, Eigen::ComputeFullV);
        MatrixD V_eigen = svd.matrixV();
        // but there could be more!! (homework!)
        size_t nbf_ribs = S_obs_ribs_ortho_eigen.cols();
        //size_t nbf_obs = S_obs_ribs_ortho.rows();
        auto nbf_cabs = nbf_ribs - svd.nonzeroSingularValues();
        MatrixD Vnull(nbf_ribs, nbf_cabs);

        //Populate Vnull with vectors of V that are orthogonal to AO space
        Vnull = V_eigen.block(0, svd.nonzeroSingularValues(), nbf_ribs, nbf_cabs);

        //Un-orthogonalize coefficients
        MatrixD C_cabs_eigen = X_ribs_eigen_inv * Vnull;

        auto tr_cabs = S_cabs.trange().data()[0];
        auto tr_ribs = S_ribs.trange().data()[0];

//        std::cout << C_cabs_eigen << std::endl;

        C_cabs = array_ops::eigen_to_array<TA::TensorD>(world, C_cabs_eigen, tr_ribs, tr_cabs);
        C_ri = array_ops::eigen_to_array<TA::TensorD>(world, X_ribs_eigen_inv, tr_ribs, tr_ribs);


        auto C_cabs_space = OrbitalSpace(OrbitalIndex(L"a'"), OrbitalIndex(L"ρ"), C_cabs);
        auto C_ribs_space = OrbitalSpace(OrbitalIndex(L"P'"), OrbitalIndex(L"ρ"), C_ri);
        mo_integral.orbital_space()->add(C_cabs_space);
        mo_integral.orbital_space()->add(C_ribs_space);

//        std::cout << C_cabs << std::endl;
//        std::cout << C_ri << std::endl;
    }


    f12::MP2F12 mp2f12(mo_integral, std::make_shared<TRange1Engine>(tre), ens);

//    mp2f12.compute_mp2_f12_c();
//    mo_integral.registry().clear();
//    ao_int.registry().clear();
    mp2f12.compute_mp2_f12_c_df();

//    ao_int.registry().print_formula(world);
//    mo_integral.registry().print_formula(world);

//    mo_integral.registry().clear();
//    ao_int.registry().clear();


    madness::finalize();
    libint2::cleanup();
    return 0;
}
