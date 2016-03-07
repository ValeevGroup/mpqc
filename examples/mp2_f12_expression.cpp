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

#include "../f12/utility.h"
#include "../integrals/integrals.h"
#include "../integrals/atomic_integral.h"
#include "../integrals/molecular_integral.h"
#include "../expression/orbital_space_registry.h"

#include "../utility/time.h"
#include "../utility/wcout_utf8.h"
#include "../utility/array_info.h"
#include "../ta_routines/array_to_eigen.h"
#include "../scf/traditional_df_fock_builder.h"

#include "../mp2/mp2.h"

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
    int nclusters = 0;
    std::cout << std::setprecision(15);
    double threshold = 1e-12;
    double well_sep_threshold = 0.1;
    integrals::QQR::well_sep_threshold(well_sep_threshold);
    if (argc == 4) {
        mol_file = argv[1];
        basis_name = argv[2];
        nclusters = std::stoi(argv[3]);
    } else {
        std::cout << "input is $./program mol_file basis_file nclusters ";
        return 0;
    }
    TiledArray::SparseShape<float>::threshold(threshold);

    auto clustered_mol = molecule::kmeans(
            molecule::read_xyz(mol_file).clusterables(), nclusters);

    auto repulsion_energy = clustered_mol.nuclear_repulsion();
    std::cout << "Nuclear Repulsion Energy: " << repulsion_energy << std::endl;
    auto occ = clustered_mol.occupation(0);

    basis::BasisSet bs(basis_name);
    basis::Basis basis(bs.get_cluster_shells(clustered_mol));

    basis::BasisSet dfbs("cc-pVTZ-RI");
    basis::Basis df_basis(dfbs.get_cluster_shells(clustered_mol));

    basis::BasisSet abs("aug-cc-pVDZ-CABS");
    basis::Basis abs_basis(abs.get_cluster_shells(clustered_mol));

    f12::GTGParams gtg_params(1.2, 6);

    libint2::init();

    integrals::AtomicIntegral<TA::TensorD, TA::SparsePolicy> ao_int
            (world,
             ta_pass_through,
             std::make_shared<molecule::Molecule>(clustered_mol),
             std::make_shared<basis::Basis>(basis),
             std::make_shared<basis::Basis>(df_basis),
             std::make_shared<basis::Basis>(abs_basis),
             gtg_params.compute()
            );

    // Overlap ints
    auto S = ao_int.compute(L"(κ|λ)");
    auto T = ao_int.compute(L"(κ| T|λ)");
    auto V = ao_int.compute(L"(κ| V|λ)");

    decltype(T) H;
    H("i,j") = T("i,j") + V("i,j");

    // Unscreened four center stored RHF.
    auto metric = ao_int.compute(L"(Κ|G|Λ)");
    scf::DFFockBuilder builder(metric);
    world.gop.fence();

    auto eri3 = ao_int.compute(L"( Κ| G|κ1 λ1)");
    scf::ClosedShellSCF<scf::DFFockBuilder> scf(H, S, occ / 2, repulsion_energy, std::move(builder));
    scf.solve(50, 1e-8, eri3);

    // obs fock build
    std::size_t all = S.trange().elements().extent()[0];
    auto tre = TRange1Engine(occ / 2, all, 4, 8, 0);

    auto F = scf.fock();
    auto L_inv = builder.inv();

    {
        TA::Array<double, 3, TA::TensorD, TA::SparsePolicy> Xab;
        Xab("X,a,b") = L_inv("X,Y") * eri3("Y,a,b");

        //ri-mp2
        auto mp2 = MP2<TA::TensorD, TA::SparsePolicy>(F, S, Xab, std::make_shared<TRange1Engine>(tre));
        mp2.compute();
    }

    //mp2
    // solve Coefficient
    std::size_t n_frozen_core = 0;
    auto F_eig = array_ops::array_to_eigen(F);
    auto S_eig = array_ops::array_to_eigen(S);
    Eig::GeneralizedSelfAdjointEigenSolver<decltype(S_eig)> es(F_eig, S_eig);
    auto ens = es.eigenvalues().bottomRows(S_eig.rows() - n_frozen_core);

    auto C_all = es.eigenvectors();
    decltype(S_eig) C_occ = C_all.block(0, 0, S_eig.rows(),occ / 2);
    decltype(S_eig) C_occ_corr = C_all.block(0, n_frozen_core, S_eig.rows(),occ / 2 - n_frozen_core);
    decltype(S_eig) C_vir = C_all.rightCols(S_eig.rows() - occ / 2);

    auto tr_0 = eri3.trange().data().back();
    auto tr_all = tre.get_all_tr1();
    auto tr_i0 = tre.get_occ_tr1();
    auto tr_vir = tre.get_vir_tr1();


    auto Ci = array_ops::eigen_to_array<TA::Tensor<double>>(world, C_occ_corr, tr_0, tr_i0);
    auto Cv = array_ops::eigen_to_array<TA::Tensor<double>>(world, C_vir, tr_0, tr_vir);
    auto Call = array_ops::eigen_to_array<TA::Tensor<double>>(world, C_all, tr_0, tr_all);

    auto orbital_registry = OrbitalSpaceRegistry<decltype(Ci)>();
    using OrbitalSpace = OrbitalSpace<decltype(Ci)>;
    auto occ_space = OrbitalSpace(OrbitalIndex(L"m"),Ci);
    orbital_registry.add(occ_space);

    auto corr_occ_space = OrbitalSpace(OrbitalIndex(L"i"),Ci);
    orbital_registry.add(corr_occ_space);

    auto vir_space = OrbitalSpace(OrbitalIndex(L"a"),Cv);
    orbital_registry.add(vir_space);

    auto obs_space = OrbitalSpace(OrbitalIndex(L"p"),Call);
    orbital_registry.add(obs_space);

    auto sp_orbital_registry = std::make_shared<decltype(orbital_registry)>(orbital_registry);

    // clean atomic integral
    ao_int.registry().clear();

    auto mo_integral = integrals::MolecularIntegral<TA::TensorD,TA::SparsePolicy>(ao_int,sp_orbital_registry);
//    mo_integral.atomic_integral().registry().print_formula();


    // test fock build
    {
        auto fock = ao_int.compute(L"(μ |F| ν)[df]");
        {
            TA::Array<double, 3, TA::TensorD, TA::SparsePolicy> Xab;
            Xab("X,a,b") = L_inv("X,Y") * eri3("Y,a,b");

            //ri-mp2
            auto mp2 = MP2<TA::TensorD, TA::SparsePolicy>(fock, S, Xab, std::make_shared<TRange1Engine>(tre));
            mp2.compute();
        }
    }


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

    {
        auto g_ijab = mo_integral.compute(L"<i j|G|a b>[df]");
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


    // construct cabs
    decltype(S_obs) C_cabs, C_ri;
    {
        auto S_obs_eigen = array_ops::array_to_eigen(S_obs);
        MatrixD X_obs_eigen = Eigen::LLT<MatrixD>(S_obs_eigen).matrixL();
        MatrixD X_obs_eigen_inv = X_obs_eigen.inverse();

        auto S_ribs_eigen = array_ops::array_to_eigen(S_ribs);
        MatrixD X_ribs_eigen = Eigen::LLT<MatrixD>(S_ribs_eigen).matrixL();
        MatrixD X_ribs_eigen_inv = X_ribs_eigen.inverse();


        MatrixD S_obs_ribs_eigen = array_ops::array_to_eigen(S_obs_ribs);

        auto S_obs_ribs_ortho_eigen = X_obs_eigen_inv.transpose() * S_obs_ribs_eigen * X_ribs_eigen_inv;

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

        C_cabs = array_ops::eigen_to_array<TA::TensorD>(world, C_cabs_eigen, tr_ribs, tr_cabs);
        C_ri = array_ops::eigen_to_array<TA::TensorD>(world,X_ribs_eigen_inv, tr_ribs, tr_ribs);


        auto C_cabs_space = OrbitalSpace(OrbitalIndex(L"a'"), C_cabs);
        auto C_ribs_space = OrbitalSpace(OrbitalIndex(L"P'"), C_ri);
        mo_integral.orbital_space()->add(C_cabs_space);
        mo_integral.orbital_space()->add(C_ribs_space);

//    std::cout << "C_cabs" << std::endl;
//    std::cout << C_cabs << std::endl;
    }

    // compute fock mo
    decltype(S_obs) F_ribs;
    {
        auto F_ribs_ribs = mo_integral.compute(L"(P'|F|Q')[df]");
        auto F_cabs_ribs = mo_integral.compute(L"(a'|F|P')[df]");
        auto F_pq = mo_integral.compute(L"(p|F|q)[df]");
        auto F_ij = mo_integral.compute(L"(i|F|j)[df]");
        auto F_mn = mo_integral.compute(L"(m|F|n)[df]");
        auto F_ri_m = mo_integral.compute(L"(P'|F|m)[df]");
        auto F_cabs_p = mo_integral.compute(L"(a'|F|p)[df]");
    }

    // V term
    decltype(S_obs) V_ijij;
    {
        V_ijij("i1,j1,i2,j2") = mo_integral(L"(i1 i2|R|j1 j2)[df]");
        V_ijij("i1,j1,i2,j2") -= mo_integral(L"(i1 p|G|j1 q)[df]")*mo_integral(L"(i2 p|R|j2 q)[df]");
        V_ijij("i1,j1,i2,j2") -= mo_integral(L"(i1 m|G|j1 a')[df]")*mo_integral(L"(j2 m|R|i2 a')[df]");
        V_ijij("i1,j1,i2,j2") -= mo_integral(L"(j1 m|G|i1 a')[df]")*mo_integral(L"(i2 m|R|j2 a')[df]");
    }



    // compute r12 integral
//    auto JK_ribs = ao_int.compute(L"(ρ1 σ1|G|ρ2 σ2)[df]");
//    auto F12_ribs = ao_int.compute(L"(ρ1 σ1|R|ρ2 σ2)");
//    auto F12_sq_ribs = ao_int.compute(L"(ρ1 σ1|R2|ρ2 σ2)");


//    auto JKF12_obs = ao_int.compute(L"(κ1 λ1 |GR|μ1 ν1)");
//    auto Comm_obs = ao_int.compute(L"( κ2 λ2 | dR2 |μ2 ν2)");


    // Coefficients

    // AO to MO transform


//    ao_int.registry().print_formula(world);
    mo_integral.registry().print_formula(world);

    madness::finalize();
    libint2::cleanup();
    return 0;
}
