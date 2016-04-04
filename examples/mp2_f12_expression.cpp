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
#include "../cc/utility.h"
#include "../integrals/integrals.h"
#include "../integrals/atomic_integral.h"
#include "../integrals/molecular_integral.h"
#include "../expression/orbital_space_registry.h"

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
             std::make_shared<basis::Basis>(basis),
             std::make_shared<basis::Basis>(df_basis),
             std::make_shared<basis::Basis>(abs_basis),
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


        auto C_cabs_space = OrbitalSpace(OrbitalIndex(L"a'"), C_cabs);
        auto C_ribs_space = OrbitalSpace(OrbitalIndex(L"P'"), C_ri);
        mo_integral.orbital_space()->add(C_cabs_space);
        mo_integral.orbital_space()->add(C_ribs_space);

//        std::cout << C_cabs << std::endl;
//        std::cout << C_ri << std::endl;
    }


    auto occ_tr1 = tre.get_occ_tr1();
    TiledArray::TiledRange occ4_trange({occ_tr1,occ_tr1,occ_tr1,occ_tr1});
    auto ijij_ijji_shape = f12::make_ijij_ijji_shape(occ4_trange);

//    {
//        auto F_mn = mo_integral.compute(L"(m|F|n)");
//        std::cout << "F_mn" << std::endl;
//        std::cout << F_mn << std::endl;
//
//        auto F_pm = mo_integral.compute(L"(P'|F|n)");
//        std::cout << "F_P'm" << std::endl;
//        std::cout << F_pm << std::endl;
//
//        auto F_Ap = mo_integral.compute(L"(a'|F|p)");
//        std::cout << "F_a'p" << std::endl;
//        std::cout << F_Ap << std::endl;
//    }

    decltype(S_obs) V_ijij_shape(world,occ4_trange,ijij_ijji_shape);
    {

        std::cout << "V_ijij_ijji_shape" << std::endl;
        V_ijij_shape("i1,j1,i2,j2") = (mo_integral(L"(Κ |GR|i2 i1)")*mo_integral(L"(Κ|GR|Λ)[inv]")*mo_integral(L"(Λ |GR|j1 j2)")).set_shape(ijij_ijji_shape);
//        std::cout << V_ijij_shape << std::endl;
        V_ijij_shape("i1,j1,i2,j2") -= (mo_integral(L"(i1 p|G|j1 q)[df]")*mo_integral(L"(i2 p|R|j2 q)[df]")).set_shape(ijij_ijji_shape);
//        std::cout << V_ijij_shape << std::endl;
        V_ijij_shape("i1,j1,i2,j2") -= (mo_integral(L"(i1 m|G|j1 a')[df]")*mo_integral(L"(m i2|R|a' j2)[df]")).set_shape(ijij_ijji_shape);
//        std::cout << V_ijij_shape << std::endl;
        V_ijij_shape("i1,j1,i2,j2") -= (mo_integral(L"(j1 m|G|i1 a')[df]")*mo_integral(L"(m j2|R|a' i2)[df]")).set_shape(ijij_ijji_shape);
        std::cout << V_ijij_shape << std::endl;

    }



    decltype(S_obs) C_iajb;
    {
        C_iajb("i,a,j,b") = mo_integral(L"(i a|R|j a')[df]")*mo_integral(L"(b|F|a')[df]");
        C_iajb("i,a,j,b") += mo_integral(L"(i a'|R|j b)[df]")*mo_integral(L"(a|F|a')[df]");
    }

    decltype(S_obs) t2;
    {
        decltype(S_obs) g_iajb;
        g_iajb = mo_integral.compute(L"(i a|G|j b)[df]");
        g_iajb("a,b,i,j") = g_iajb("i,a,j,b");
        t2 = mpqc::cc::d_abij(g_iajb,ens,occ/2);
    }


    decltype(V_ijij_shape) V_bar_ijij_shape;
    {

        V_bar_ijij_shape("i1,j1,i2,j2") = (V_ijij_shape("i1,j1,i2,j2") + C_iajb("i1,a,j1,b")*t2("a,b,i2,j2")).set_shape(ijij_ijji_shape);
        std::cout << "V_bar_ijij_ijji" << std::endl;
        std::cout << V_bar_ijij_shape << std::endl;
    }

    double E21 = 0.0;
    {
        // diagonal sum
        E21 = V_bar_ijij_shape("i1,j1,i2,j2").reduce(f12::DiagonalSum<TA::TensorD>());
        std::cout << E21 << std::endl;

//        V_bar_ijji_shape("i1,j1,i2,j2") = V_bar_ijji_shape("i1,j1,j2,i2");
        // off diagonal sum
        E21 += 0.5*(5*V_bar_ijij_shape("i1,j1,i2,j2")-V_bar_ijij_shape("i1,j1,j2,i2")).reduce(f12::OffDiagonalSum<TA::TensorD>());
        std::cout << "E21: " << E21 << std::endl;

    }

    decltype(V_ijij_shape) X_ijij_shape;
    {
        std::cout << "X_ijij_shape" << std::endl;
        X_ijij_shape("i1,j1,i2,j2") = (mo_integral(L"(Κ |R2|i1 i2)")*mo_integral(L"(Κ|R2|Λ)[inv]")*mo_integral(L"(Λ |R2|j1 j2)")).set_shape(ijij_ijji_shape);
//        std::cout << X_ijij_shape << std::endl;
        X_ijij_shape("i1,j1,i2,j2") -= (mo_integral(L"(i1 p|R|j1 q)[df]")*mo_integral(L"(i2 p|R|j2 q)[df]")).set_shape(ijij_ijji_shape);
//        std::cout << X_ijij_shape << std::endl;
        X_ijij_shape("i1,j1,i2,j2") -= (mo_integral(L"(i1 m|R|j1 a')[df]")*mo_integral(L"(m i2|R|a' j2)[df]")).set_shape(ijij_ijji_shape);
//        std::cout << X_ijij_shape << std::endl;
        X_ijij_shape("i1,j1,i2,j2") -= (mo_integral(L"(j1 m|R|i1 a')[df]")*mo_integral(L"(m j2|R|a' i2)[df]")).set_shape(ijij_ijji_shape);
        std::cout << X_ijij_shape << std::endl;
    }

    decltype(S_obs) B_ijij_shape;
    {

//        hJ = mo_integral(L"(i|V|P')") + mo_integral(L"(i|T|P')") + mo_integral(L"(i|J|P')[df]");
        auto hJ = mo_integral.compute(L"(P' | hJ | i)[df]");

        std::cout << "B_ijij_shape" << std::endl;
        B_ijij_shape("i1,j1,i2,j2") = (mo_integral(L"(Κ |dR2|i1 i2)")*mo_integral(L"(Κ|dR2|Λ)[inv]")*mo_integral(L"(Λ |dR2|j1 j2)")).set_shape(ijij_ijji_shape);

//        std::cout << B_ijij_shape << std::endl;
        B_ijij_shape("i1,j1,i2,j2") += (mo_integral(L"(i1 P'|R2|j1 j2)[df]")*hJ("P',i2")).set_shape(ijij_ijji_shape);
        B_ijij_shape("i1,j1,i2,j2") += (mo_integral(L"(j1 P'|R2|i1 i2)[df]")*hJ("P',j2")).set_shape(ijij_ijji_shape);

//        std::cout << B_ijij_shape << std::endl;

        B_ijij_shape("i1,j1,i2,j2") -= (mo_integral(L"(i1 P'|R|j1 Q')[df]")*mo_integral(L"(P'|K|R')[df]")*mo_integral(L"(R' i2|R|Q' j2)[df]")).set_shape(ijij_ijji_shape);
        B_ijij_shape("i1,j1,i2,j2") -= (mo_integral(L"(j1 P'|R|i1 Q')[df]")*mo_integral(L"(P'|K|R')[df]")*mo_integral(L"(R' j2|R|Q' i2)[df]")).set_shape(ijij_ijji_shape);

//        std::cout << B_ijij_shape << std::endl;
        B_ijij_shape("i1,j1,i2,j2") -= (mo_integral(L"(i1 P'|R|j1 m)[df]")*mo_integral(L"(P'|F|R')[df]")*mo_integral(L"(R' i2|R|m j2)[df]")).set_shape(ijij_ijji_shape);
        B_ijij_shape("i1,j1,i2,j2") -= (mo_integral(L"(j1 P'|R|i1 m)[df]")*mo_integral(L"(P'|F|R')[df]")*mo_integral(L"(R' j2|R|m i2)[df]")).set_shape(ijij_ijji_shape);


//        std::cout << B_ijij_shape << std::endl;
        B_ijij_shape("i1,j1,i2,j2") -= (mo_integral(L"(i1 p|R|j1 a)[df]")*mo_integral(L"(p|F|r)[df]")*mo_integral(L"(i2 r|R|j2 a)[df]")).set_shape(ijij_ijji_shape);
        B_ijij_shape("i1,j1,i2,j2") -= (mo_integral(L"(j1 p|R|i1 a)[df]")*mo_integral(L"(p|F|r)[df]")*mo_integral(L"(j2 r|R|i2 a)[df]")).set_shape(ijij_ijji_shape);


//        std::cout << B_ijij_shape << std::endl;
        B_ijij_shape("i1,j1,i2,j2") += (mo_integral(L"(i1 m|R|j1 b')[df]")*mo_integral(L"(m|F|n)[df]")*mo_integral(L"(n i2|R|b' j2)[df]")).set_shape(ijij_ijji_shape);
        B_ijij_shape("i1,j1,i2,j2") += (mo_integral(L"(j1 m|R|i1 b')[df]")*mo_integral(L"(m|F|n)[df]")*mo_integral(L"(n j2|R|b' i2)[df]")).set_shape(ijij_ijji_shape);

//        std::cout << B_ijij_shape << std::endl;

        B_ijij_shape("i1,j1,i2,j2") -= (2.0*mo_integral(L"(i1 m|R|j1 b')[df]")*mo_integral(L"(m|F|P')[df]")*mo_integral(L"(P' i2|R|b' j2)[df]")).set_shape(ijij_ijji_shape);
        B_ijij_shape("i1,j1,i2,j2") -= (2.0*mo_integral(L"(j1 m|R|i1 b')[df]")*mo_integral(L"(m|F|P')[df]")*mo_integral(L"(P' j2|R|b' i2)[df]")).set_shape(ijij_ijji_shape);

//        std::cout << B_ijij_shape << std::endl;

        B_ijij_shape("i1,j1,i2,j2") -= (2.0*mo_integral(L"(i1 p|R|j1 a)[df]")*mo_integral(L"(p|F|a')[df]")*mo_integral(L"(i2 a'|R|j2 a)[df]")).set_shape(ijij_ijji_shape);
        B_ijij_shape("i1,j1,i2,j2") -= (2.0*mo_integral(L"(j1 p|R|i1 a)[df]")*mo_integral(L"(p|F|a')[df]")*mo_integral(L"(j2 a'|R|i2 a)[df]")).set_shape(ijij_ijji_shape);

        std::cout << B_ijij_shape << std::endl;
    }


    {
        auto Fij = mo_integral.compute(L"(i|F|j)[df]");
        auto Fij_eigen = array_ops::array_to_eigen(Fij);
        f12::convert_X_ijkl(X_ijij_shape, Fij_eigen);

        auto C_bar_iajb = f12::convert_C_iajb(C_iajb, occ/2, ens);

        B_ijij_shape("i1,j1,i2,j2") = B_ijij_shape("i1,j1,i2,j2") - X_ijij_shape("i1,j1,i2,j2") - (C_iajb("i1,a,j1,b")*C_bar_iajb("i2,a,j2,b")).set_shape(ijij_ijji_shape);
        std::cout << "B bar ijij_ijji Shape" << std::endl;
        std::cout << B_ijij_shape << std::endl;
    }

    double E22 = 0.0;
    {
        // diagonal sum
        E22 = 0.25*B_ijij_shape("i1,j1,i2,j2").reduce(f12::DiagonalSum<TA::TensorD>());
        // off diagonal sum
        E22 += 0.0625 * (7 * B_ijij_shape("i1,j1,i2,j2") - B_ijij_shape("i1,j1,j2,i2")).reduce(f12::OffDiagonalSum<TA::TensorD>());
        std::cout << "E22: " << E22 << std::endl;
    }
    std::cout << "E_F12: " << E22+E21 << std::endl;

//    ao_int.registry().print_formula(world);
//    mo_integral.registry().print_formula(world);

    mo_integral.registry().clear();
    ao_int.registry().clear();

    // without df

    decltype(S_obs) V_ijij_nodf_shape(world,occ4_trange,ijij_ijji_shape);
    {

        std::cout << "V_ijij_nodf_shape" << std::endl;
        V_ijij_nodf_shape("i1,j1,i2,j2") = mo_integral(L"(i1 i2|GR|j1 j2)").set_shape(ijij_ijji_shape);
//        std::cout << V_ijij_nodf_shape << std::endl;
        V_ijij_nodf_shape("i1,j1,i2,j2") -= (mo_integral(L"(i1 p|G|j1 q)")*mo_integral(L"(i2 p|R|j2 q)")).set_shape(ijij_ijji_shape);
//        std::cout << V_ijij_nodf_shape << std::endl;
        V_ijij_nodf_shape("i1,j1,i2,j2") -= (mo_integral(L"(i1 m|G|j1 a')")*mo_integral(L"(m i2|R|a' j2)")).set_shape(ijij_ijji_shape);
//        std::cout << V_ijij_nodf_shape << std::endl;
        V_ijij_nodf_shape("i1,j1,i2,j2") -= (mo_integral(L"(j1 m|G|i1 a')")*mo_integral(L"(m j2|R|a' i2)")).set_shape(ijij_ijji_shape);
//        std::cout << V_ijij_nodf_shape << std::endl;

    }

    decltype(S_obs) C_iajb_nodf;
    {
        C_iajb_nodf("i,a,j,b") = mo_integral(L"(i a|R|j a')")*mo_integral(L"(b|F|a')");
        C_iajb_nodf("i,a,j,b") += mo_integral(L"(i a'|R|j b)")*mo_integral(L"(a|F|a')");
    }

    decltype(S_obs) t2_nodf;
    {
        decltype(S_obs) g_iajb;
        g_iajb = mo_integral.compute(L"(i a|G|j b)");
        g_iajb("a,b,i,j") = g_iajb("i,a,j,b");
        t2_nodf = mpqc::cc::d_abij(g_iajb,ens,occ/2);
    }


    decltype(V_ijij_shape) V_bar_ijij_shape_nodf;
    {
        V_bar_ijij_shape_nodf("i1,j1,i2,j2") = (V_ijij_nodf_shape("i1,j1,i2,j2") + C_iajb_nodf("i1,a,j1,b")*t2_nodf("a,b,i2,j2")).set_shape(ijij_ijji_shape);
        std::cout << "V_bar_ijij_ijji_nodf" << std::endl;
        std::cout << V_bar_ijij_shape_nodf << std::endl;
    }


    E21 = 0.0;
    {
        // diagonal sum
        E21 = V_bar_ijij_shape_nodf("i1,j1,i2,j2").reduce(f12::DiagonalSum<TA::TensorD>());
        std::cout << E21 << std::endl;

        // off diagonal sum
        E21 += 0.5*(5*V_bar_ijij_shape_nodf("i1,j1,i2,j2")-V_bar_ijij_shape_nodf("i1,j1,j2,i2")).reduce(f12::OffDiagonalSum<TA::TensorD>());
        std::cout << "E21: " << E21 << std::endl;

    }

    decltype(V_ijij_shape) X_ijij_nodf_shape;
    {
        std::cout << "X_ijij_nodf_shape" << std::endl;
        X_ijij_nodf_shape("i1,j1,i2,j2") = mo_integral(L"(i1 i2 |R2|j1 j2)").set_shape(ijij_ijji_shape);
//        std::cout << X_ijij_nodf_shape << std::endl;
        X_ijij_nodf_shape("i1,j1,i2,j2") -= (mo_integral(L"(i1 p|R|j1 q)")*mo_integral(L"(i2 p|R|j2 q)")).set_shape(ijij_ijji_shape);
//        std::cout << X_ijij_nodf_shape << std::endl;
        X_ijij_nodf_shape("i1,j1,i2,j2") -= (mo_integral(L"(i1 m|R|j1 a')")*mo_integral(L"(m i2|R|a' j2)")).set_shape(ijij_ijji_shape);
//        std::cout << X_ijij_shape << std::endl;
        X_ijij_nodf_shape("i1,j1,i2,j2") -= (mo_integral(L"(j1 m|R|i1 a')")*mo_integral(L"(m j2|R|a' i2)")).set_shape(ijij_ijji_shape);
//        std::cout << X_ijij_nodf_shape << std::endl;
    }

    decltype(S_obs) B_ijij_shape_nodf;
    {

//        hJ = mo_integral(L"(i|V|P')") + mo_integral(L"(i|T|P')") + mo_integral(L"(i|J|P')[df]");
        auto hJ = mo_integral.compute(L"(P' | hJ | i)");

        std::cout << "B_ijij_shape_nodf" << std::endl;
        B_ijij_shape_nodf("i1,j1,i2,j2") = (mo_integral(L"(i1 i2 |dR2|j1 j2)")).set_shape(ijij_ijji_shape);

//        std::cout << B_ijij_shape << std::endl;
        B_ijij_shape_nodf("i1,j1,i2,j2") += (mo_integral(L"(i1 P'|R2|j1 j2)")*hJ("P',i2")).set_shape(ijij_ijji_shape);
        B_ijij_shape_nodf("i1,j1,i2,j2") += (mo_integral(L"(j1 P'|R2|i1 i2)")*hJ("P',j2")).set_shape(ijij_ijji_shape);

//        std::cout << B_ijij_shape << std::endl;

        B_ijij_shape_nodf("i1,j1,i2,j2") -= (mo_integral(L"(i1 P'|R|j1 Q')")*mo_integral(L"(P'|K|R')")*mo_integral(L"(R' i2|R|Q' j2)")).set_shape(ijij_ijji_shape);
        B_ijij_shape_nodf("i1,j1,i2,j2") -= (mo_integral(L"(j1 P'|R|i1 Q')")*mo_integral(L"(P'|K|R')")*mo_integral(L"(R' j2|R|Q' i2)")).set_shape(ijij_ijji_shape);

//        std::cout << B_ijij_shape << std::endl;
        B_ijij_shape_nodf("i1,j1,i2,j2") -= (mo_integral(L"(i1 P'|R|j1 m)")*mo_integral(L"(P'|F|R')")*mo_integral(L"(R' i2|R|m j2)")).set_shape(ijij_ijji_shape);
        B_ijij_shape_nodf("i1,j1,i2,j2") -= (mo_integral(L"(j1 P'|R|i1 m)")*mo_integral(L"(P'|F|R')")*mo_integral(L"(R' j2|R|m i2)")).set_shape(ijij_ijji_shape);


//        std::cout << B_ijij_shape << std::endl;
        B_ijij_shape_nodf("i1,j1,i2,j2") -= (mo_integral(L"(i1 p|R|j1 a)")*mo_integral(L"(p|F|r)")*mo_integral(L"(i2 r|R|j2 a)")).set_shape(ijij_ijji_shape);
        B_ijij_shape_nodf("i1,j1,i2,j2") -= (mo_integral(L"(j1 p|R|i1 a)")*mo_integral(L"(p|F|r)")*mo_integral(L"(j2 r|R|i2 a)")).set_shape(ijij_ijji_shape);


//        std::cout << B_ijij_shape << std::endl;
        B_ijij_shape_nodf("i1,j1,i2,j2") += (mo_integral(L"(i1 m|R|j1 b')")*mo_integral(L"(m|F|n)")*mo_integral(L"(n i2|R|b' j2)")).set_shape(ijij_ijji_shape);
        B_ijij_shape_nodf("i1,j1,i2,j2") += (mo_integral(L"(j1 m|R|i1 b')")*mo_integral(L"(m|F|n)")*mo_integral(L"(n j2|R|b' i2)")).set_shape(ijij_ijji_shape);

//        std::cout << B_ijij_shape << std::endl;

        B_ijij_shape_nodf("i1,j1,i2,j2") -= (2.0*mo_integral(L"(i1 m|R|j1 b')")*mo_integral(L"(m|F|P')")*mo_integral(L"(P' i2|R|b' j2)")).set_shape(ijij_ijji_shape);
        B_ijij_shape_nodf("i1,j1,i2,j2") -= (2.0*mo_integral(L"(j1 m|R|i1 b')")*mo_integral(L"(m|F|P')")*mo_integral(L"(P' j2|R|b' i2)")).set_shape(ijij_ijji_shape);

//        std::cout << B_ijij_shape << std::endl;

        B_ijij_shape_nodf("i1,j1,i2,j2") -= (2.0*mo_integral(L"(i1 p|R|j1 a)")*mo_integral(L"(p|F|a')")*mo_integral(L"(i2 a'|R|j2 a)")).set_shape(ijij_ijji_shape);
        B_ijij_shape_nodf("i1,j1,i2,j2") -= (2.0*mo_integral(L"(j1 p|R|i1 a)")*mo_integral(L"(p|F|a')")*mo_integral(L"(j2 a'|R|i2 a)")).set_shape(ijij_ijji_shape);

        std::cout << B_ijij_shape_nodf << std::endl;
    }

    {
        auto Fij = mo_integral.compute(L"(i|F|j)");
        auto Fij_eigen = array_ops::array_to_eigen(Fij);
        f12::convert_X_ijkl(X_ijij_nodf_shape, Fij_eigen);

        auto C_bar_iajb = f12::convert_C_iajb(C_iajb_nodf, occ/2, ens);

        B_ijij_shape_nodf("i1,j1,i2,j2") = B_ijij_shape_nodf("i1,j1,i2,j2") - X_ijij_nodf_shape("i1,j1,i2,j2") - (C_iajb_nodf("i1,a,j1,b")*C_bar_iajb("i2,a,j2,b")).set_shape(ijij_ijji_shape);
        std::cout << "B bar ijij_ijji Shape" << std::endl;
        std::cout << B_ijij_shape_nodf << std::endl;
    }

    E22 = 0.0;
    {
        // diagonal sum
        E22 = 0.25*B_ijij_shape_nodf("i1,j1,i2,j2").reduce(f12::DiagonalSum<TA::TensorD>());
        // off diagonal sum
        E22 += 0.0625 * (7 * B_ijij_shape_nodf("i1,j1,i2,j2") - B_ijij_shape_nodf("i1,j1,j2,i2")).reduce(f12::OffDiagonalSum<TA::TensorD>());
        std::cout << "E22: " << E22 << std::endl;
    }
    std::cout << "E_F12: " << E22+E21 << std::endl;

    madness::finalize();
    libint2::cleanup();
    return 0;
}
