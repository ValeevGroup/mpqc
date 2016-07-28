//
// Created by Chong Peng on 4/14/16.
//

#ifndef MPQC_MO_BUILD_H
#define MPQC_MO_BUILD_H


#include <rapidjson/document.h>

#include "../../../../../include/tiledarray.h"
#include "../../../../../common/namespaces.h"
#include "../../../../../utility/trange1_engine.h"
#include <mpqc/chemistry/qc/integrals/molecular_integral.h>
#include <mpqc/chemistry/qc/expression/orbital_registry.h>

namespace mpqc {


template <typename Tile, typename Policy>
std::shared_ptr<TRange1Engine> closed_shell_obs_mo_build_eigen_solve(
        integrals::MolecularIntegral<Tile, Policy> &mo_int,
        Eigen::VectorXd& ens,
        const rapidjson::Document& in,
        const molecule::Molecule& mols,
        int occ)
{
    auto& ao_int = mo_int.atomic_integral();
    auto orbital_registry = mo_int.orbital_space();
    auto& world = ao_int.world();
    using TArray = TA::DistArray<Tile, Policy>;


    auto mo_time0 = mpqc_time::fenced_now(world);
    utility::print_par(world, "\nBuilding ClosedShell OBS MO Orbital\n");

    // find fock matrix
    TArray F;
    if(ao_int.registry().have(Formula(L"<μ|F|ν>"))){
        F = ao_int.registry().retrieve(Formula(L"<μ|F|ν>"));
    }
    else{
        F = ao_int.registry().retrieve(Formula(L"<μ|F|ν>[df]"));
    }

    auto S = ao_int.compute(L"<κ|λ>");

    MatrixD F_eig = array_ops::array_to_eigen(F);
    MatrixD S_eig = array_ops::array_to_eigen(S);

    // check the condition number in Overlap
    Eig::SelfAdjointEigenSolver<MatrixD> S_es(S_eig);
    // eigen value in increasing order
    auto cond = S_es.eigenvalues()(S_es.eigenvalues().size() - 1) / S_es.eigenvalues()(0);
    utility::print_par(world, "Condition Number in Overlap: ", cond, "\n");


    // solve mo coefficients
    Eig::GeneralizedSelfAdjointEigenSolver<MatrixD> es(F_eig, S_eig);

    // start to solve coefficient
    bool frozen_core = in.HasMember("FrozenCore") ? in["FrozenCore"].GetBool() : false;
    std::size_t n_frozen_core = 0;
    if (frozen_core) {
        n_frozen_core = mols.core_electrons();
        utility::print_par(world, "Frozen Core: ", n_frozen_core," electrons", "\n");
        n_frozen_core = n_frozen_core / 2;
    }

    ens = es.eigenvalues();
    MatrixD C_all = es.eigenvectors();
    MatrixD C_occ = C_all.block(0, 0, S_eig.rows(),occ);
    MatrixD C_corr_occ = C_all.block(0, n_frozen_core, S_eig.rows(), occ - n_frozen_core);
    MatrixD C_vir = C_all.rightCols(S_eig.rows() - occ);

    // get all the sizes
    std::size_t mo_blocksize = in.HasMember("MoBlockSize") ? in["MoBlockSize"].GetInt() : 24;
    std::size_t occ_blocksize = in.HasMember("OccBlockSize") ? in["OccBlockSize"].GetInt() : mo_blocksize;
    std::size_t vir_blocksize = in.HasMember("VirBlockSize") ? in["VirBlockSize"].GetInt() : mo_blocksize;

    utility::print_par(world,"OccBlockSize: ", occ_blocksize, "\n");
    utility::print_par(world,"VirBlockSize: ", vir_blocksize, "\n");


    std::size_t all = S.trange().elements().extent()[0];
    auto tre = std::make_shared<TRange1Engine>(occ, all, occ_blocksize, vir_blocksize, n_frozen_core);

    // get all the trange1s
    auto tr_obs = S.trange().data().back();
    auto tr_corr_occ = tre->get_occ_tr1();
    auto tr_occ = tre->compute_range(occ, occ_blocksize);
    auto tr_vir = tre->get_vir_tr1();
    auto tr_all = tre->get_all_tr1();

    utility::parallel_print_range_info(world, tr_occ, "Occ");
    utility::parallel_print_range_info(world, tr_corr_occ, "CorrOcc");
    utility::parallel_print_range_info(world, tr_vir, "Vir");
    utility::parallel_print_range_info(world, tr_all, "Obs");

    // convert to TA
    auto C_occ_ta = array_ops::eigen_to_array<Tile>(world, C_occ, tr_obs, tr_occ);
    auto C_corr_occ_ta = array_ops::eigen_to_array<Tile>(world, C_corr_occ, tr_obs, tr_corr_occ);
    auto C_vir_ta = array_ops::eigen_to_array<Tile>(world, C_vir, tr_obs, tr_vir);
    auto C_all_ta = array_ops::eigen_to_array<Tile>(world, C_all, tr_obs, tr_all);

    // insert to registry
    using OrbitalSpaceTArray = OrbitalSpace<TA::DistArray<Tile,Policy>>;
    auto occ_space = OrbitalSpaceTArray(OrbitalIndex(L"m"), OrbitalIndex(L"κ"), C_occ_ta);
    orbital_registry->add(occ_space);

    auto corr_occ_space = OrbitalSpaceTArray(OrbitalIndex(L"i"), OrbitalIndex(L"κ"), C_corr_occ_ta);
    orbital_registry->add(corr_occ_space);

    auto vir_space = OrbitalSpaceTArray(OrbitalIndex(L"a"), OrbitalIndex(L"κ"), C_vir_ta);
    orbital_registry->add(vir_space);

    auto obs_space = OrbitalSpaceTArray(OrbitalIndex(L"p"),OrbitalIndex(L"κ"), C_all_ta);
    orbital_registry->add(obs_space);

    auto mo_time1 = mpqc_time::fenced_now(world);
    auto mo_time = mpqc_time::duration_in_s(mo_time0,mo_time1);
    utility::print_par(world,"ClosedShell OBS MO Build Time: ", mo_time, " S\n");

    return tre;
};

template <typename Tile, typename Policy>
void closed_shell_cabs_mo_build_svd(
        integrals::MolecularIntegral <Tile, Policy> &mo_int,
        const rapidjson::Document &in,
        const std::shared_ptr<TRange1Engine> tre) {
    auto& ao_int = mo_int.atomic_integral();
    auto orbital_registry = mo_int.orbital_space();
    auto &world = ao_int.world();
    // CABS fock build
    auto mo_time0 = mpqc_time::fenced_now(world);
    utility::print_par(world, "\nBuilding ClosedShell CABS MO Orbital\n");

    // integral
    auto S_cabs = ao_int.compute(L"<α|β>");
    auto S_ribs = ao_int.compute(L"<ρ|σ>");
    auto S_obs_ribs = ao_int.compute(L"<μ|σ>");
    auto S_obs = ao_int.compute(L"<κ|λ>");

    // construct cabs
    TA::DistArray <Tile, Policy> C_cabs, C_ri, C_allvir;
    {
        MatrixD S_obs_eigen = array_ops::array_to_eigen(S_obs);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(S_obs_eigen);
        MatrixD X_obs_eigen_inv = es.operatorInverseSqrt();

        MatrixD S_ribs_eigen = array_ops::array_to_eigen(S_ribs);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es2(S_ribs_eigen);
        MatrixD X_ribs_eigen_inv = es2.operatorInverseSqrt();

        // orthogonalize
        MatrixD S_obs_ribs_eigen = array_ops::array_to_eigen(S_obs_ribs);
        MatrixD S_obs_ribs_ortho_eigen = X_obs_eigen_inv.transpose() * S_obs_ribs_eigen * X_ribs_eigen_inv;

        // SVD solve
        Eigen::JacobiSVD<MatrixD> svd(S_obs_ribs_ortho_eigen, Eigen::ComputeFullV);
        MatrixD V_eigen = svd.matrixV();
        size_t nbf_ribs = S_obs_ribs_ortho_eigen.cols();
        auto nbf_cabs = nbf_ribs - svd.nonzeroSingularValues();
        MatrixD Vnull(nbf_ribs, nbf_cabs);
        Vnull = V_eigen.block(0, svd.nonzeroSingularValues(), nbf_ribs, nbf_cabs);
        MatrixD C_cabs_eigen = X_ribs_eigen_inv * Vnull;



        // solve orbitals for all virtual

        auto n_vir = tre->get_vir();
        MatrixD C_allvirtual_eigen = MatrixD::Zero(nbf_ribs, n_vir+nbf_cabs);

        {

          auto C_vir = orbital_registry->retrieve(OrbitalIndex(L"a")).array();
          MatrixD C_vir_eigen = array_ops::array_to_eigen(C_vir);

          auto n_obs = C_vir_eigen.rows();

          C_allvirtual_eigen.block(0,0,n_obs,n_vir) << C_vir_eigen;
          C_allvirtual_eigen.block(0,n_vir,nbf_ribs, nbf_cabs) << C_cabs_eigen;
        }

        // get mo block size
        std::size_t mo_blocksize = in.HasMember("MoBlockSize") ? in["MoBlockSize"].GetInt() : 24;
        std::size_t vir_blocksize = in.HasMember("VirBlockSize") ? in["VirBlockSize"].GetInt() : mo_blocksize;
        utility::print_par(world, "VirBlockSize: ", vir_blocksize, "\n");

        auto tr_ribs = S_ribs.trange().data()[0];
        auto tr_cabs = S_cabs.trange().data()[0];
        auto tr_cabs_mo = tre->compute_range(nbf_cabs, vir_blocksize);
        auto tr_ribs_mo = tre->compute_range(nbf_ribs, vir_blocksize);
        auto tr_allvir_mo = tre->compute_range(nbf_cabs+n_vir, vir_blocksize);


        utility::parallel_print_range_info(world, tr_cabs_mo, "CABS MO");
        utility::parallel_print_range_info(world, tr_allvir_mo, "All Virtual MO");
        utility::parallel_print_range_info(world, tr_ribs_mo, "RIBS MO");

        C_cabs = array_ops::eigen_to_array<TA::TensorD>(world, C_cabs_eigen, tr_ribs, tr_cabs_mo);
        C_ri = array_ops::eigen_to_array<TA::TensorD>(world, X_ribs_eigen_inv, tr_ribs, tr_ribs_mo);
        C_allvir = array_ops::eigen_to_array<TA::TensorD>(world, C_allvirtual_eigen, tr_ribs, tr_allvir_mo);

        // insert to orbital space
        using OrbitalSpaceTArray = OrbitalSpace<TA::DistArray < Tile, Policy>>;
        auto C_cabs_space = OrbitalSpaceTArray(OrbitalIndex(L"a'"), OrbitalIndex(L"ρ"), C_cabs);
        auto C_ribs_space = OrbitalSpaceTArray(OrbitalIndex(L"P'"), OrbitalIndex(L"ρ"), C_ri);
        auto C_allvir_space = OrbitalSpaceTArray(OrbitalIndex(L"A'"), OrbitalIndex(L"ρ"), C_allvir);

        orbital_registry->add(C_cabs_space);
        orbital_registry->add(C_ribs_space);
        orbital_registry->add(C_allvir_space);

        auto mo_time1 = mpqc_time::fenced_now(world);
        auto mo_time = mpqc_time::duration_in_s(mo_time0, mo_time1);
        utility::print_par(world, "ClosedShell CABS MO Build Time: ", mo_time, " S\n");
    }

};

template <typename Tile, typename Policy>
std::shared_ptr<TRange1Engine> closed_shell_dualbasis_mo_build_eigen_solve_svd(
        integrals::MolecularIntegral <Tile, Policy> &mo_int,
        Eigen::VectorXd &ens, const rapidjson::Document &in,
        const molecule::Molecule &mols, int occ) {
        auto &ao_int = mo_int.atomic_integral();
        auto &world = ao_int.world();
        using TArray = TA::DistArray<Tile, Policy>;

        utility::print_par(world, "\nBuilding ClosedShell Dual Basis MO Orbital\n");
        auto mo_time0 = mpqc_time::fenced_now(world);

        // solving occupied orbitals
        TArray F;
        if (ao_int.registry().have(Formula(L"<μ|F|ν>"))) {
            F = ao_int.registry().retrieve(Formula(L"<μ|F|ν>"));
        } else {
            F = ao_int.registry().retrieve(Formula(L"<μ|F|ν>[df]"));
        }

        auto S = ao_int.compute(L"<κ|λ>");

        MatrixD F_eig = array_ops::array_to_eigen(F);
        MatrixD S_eig = array_ops::array_to_eigen(S);

        // solve mo coefficients
        Eig::GeneralizedSelfAdjointEigenSolver<MatrixD> es(F_eig, S_eig);

        bool frozen_core =
                in.HasMember("FrozenCore") ? in["FrozenCore"].GetBool() : false;
        std::size_t n_frozen_core = 0;
        if (frozen_core) {
            n_frozen_core = mols.core_electrons();
            utility::print_par(world, "Frozen Core: ", n_frozen_core, " electrons",
                               "\n");
            n_frozen_core = n_frozen_core / 2;
        }

        Eigen::VectorXd ens_occ = es.eigenvalues().segment(0, occ);
//        std::cout << ens_occ << std::endl;
        MatrixD C_all = es.eigenvectors();
        MatrixD C_occ = C_all.block(0, 0, S_eig.rows(), occ);
        MatrixD C_corr_occ = C_all.block(0, n_frozen_core, S_eig.rows(), occ - n_frozen_core);

        // finished solving occupied orbitals


        // start to solve virtual orbitals
        auto S_vbs = ao_int.compute(L"<Α|Β>");
        auto S_obs_vbs = ao_int.compute(L"<μ|Α>");

        // construct C_vbs
        MatrixD C_vbs;
        std::size_t nbf_vbs;
        std::size_t nbf_v;
        {

            // S_A^B -(1/2)
            MatrixD S_vbs_eigen = array_ops::array_to_eigen(S_vbs);
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es2(S_vbs_eigen);
            MatrixD X_vbs_eigen_inv = es2.operatorInverseSqrt();

            // orthogonalize
            MatrixD S_obs_vbs_eigen = array_ops::array_to_eigen(S_obs_vbs);
            MatrixD S_obs_vbs_ortho_eigen = S_obs_vbs_eigen * X_vbs_eigen_inv;

            MatrixD X_i_mu = C_occ.transpose()*S_obs_vbs_ortho_eigen;

            // SVD solve
            Eigen::JacobiSVD<MatrixD> svd(X_i_mu, Eigen::ComputeFullV);
            MatrixD V_eigen = svd.matrixV();
            nbf_vbs = S_obs_vbs_ortho_eigen.cols();
            nbf_v = nbf_vbs - svd.nonzeroSingularValues();
            std::cout << "Virtual Orbital Size:  " << nbf_v << std::endl;
            MatrixD Vnull(nbf_vbs, nbf_v);
            Vnull = V_eigen.block(0, svd.nonzeroSingularValues(), nbf_vbs, nbf_v);
            C_vbs = X_vbs_eigen_inv.transpose() * Vnull;

        }


        // get all the sizes
        std::size_t mo_blocksize = in.HasMember("MoBlockSize") ? in["MoBlockSize"].GetInt() : 24;
        std::size_t vir_blocksize = in.HasMember("VirBlockSize")
                                    ? in["VirBlockSize"].GetInt()
                                    : mo_blocksize;
        std::size_t occ_blocksize = in.HasMember("OccBlockSize")
                                    ? in["OccBlockSize"].GetInt()
                                    : mo_blocksize;
        utility::print_par(world, "OccBlockSize: ", occ_blocksize, "\n");
        utility::print_par(world, "VirBlockSize: ", vir_blocksize, "\n");

        auto tre = std::make_shared<TRange1Engine>(occ, nbf_vbs, occ_blocksize,
                                                   vir_blocksize, n_frozen_core);
        auto tr_obs = S.trange().data().back();
        auto tr_vbs = S_vbs.trange().data().back();
        auto tr_occ = tre->compute_range(occ, occ_blocksize);
        auto tr_corr_occ = tre->get_occ_tr1();
        auto tr_vir = tre->get_vir_tr1();

        utility::parallel_print_range_info(world, tr_occ, "Occ");
        utility::parallel_print_range_info(world, tr_corr_occ, "CorrOcc");
        utility::parallel_print_range_info(world, tr_vir, "Vir");

        // convert to TA
        auto C_occ_ta = array_ops::eigen_to_array<Tile>(world, C_occ, tr_obs, tr_occ);
        auto C_corr_occ_ta = array_ops::eigen_to_array<Tile>(world, C_corr_occ, tr_obs, tr_corr_occ);
        auto C_vir_ta = array_ops::eigen_to_array<Tile>(world, C_vbs, tr_vbs, tr_vir);

        // insert to registry
        using OrbitalSpaceTArray = OrbitalSpace<TA::DistArray<Tile, Policy>>;
        auto occ_space = OrbitalSpaceTArray(OrbitalIndex(L"m"), OrbitalIndex(L"κ"), C_occ_ta);
        mo_int.orbital_space()->add(occ_space);

        auto corr_occ_space = OrbitalSpaceTArray(OrbitalIndex(L"i"), OrbitalIndex(L"κ"), C_corr_occ_ta);
        mo_int.orbital_space()->add(corr_occ_space);

        auto vir_space = OrbitalSpaceTArray(OrbitalIndex(L"a"), OrbitalIndex(L"Α"), C_vir_ta);
        mo_int.orbital_space()->add(vir_space);


        // solve energy in virtual orbital
        TArray F_vbs;
        if (ao_int.registry().have(Formula(L"<μ|F|ν>"))) {
            F_vbs = mo_int.compute(L"<a|F|b>");
        } else {
            F_vbs = mo_int.compute(L"<a|F|b>[df]");
        }
        MatrixD F_vbs_mo_eigen = array_ops::array_to_eigen(F_vbs);
//    std::cout << "F_vbs MO" << std::endl;
//    std::cout << F_vbs_mo_eigen << std::endl;

        Eigen::SelfAdjointEigenSolver<MatrixD> es3(F_vbs_mo_eigen);
        auto ens_vir = es3.eigenvalues();

        ens = Eigen::VectorXd(nbf_vbs);
        ens << ens_occ, ens_vir;

//        std::cout << "Energy of Orbitals " << std::endl;
//        std::cout << ens << std::endl;

        // resolve the virtual orbitals
        MatrixD C_vir_rotate = es3.eigenvectors();
        C_vbs = C_vbs*C_vir_rotate;
        TArray C_vir_ta_new = array_ops::eigen_to_array<Tile>(world, C_vbs, tr_vbs, tr_vir);

        // remove old virtual orbitals
        mo_int.orbital_space()->remove(OrbitalIndex(L"a"));
        mo_int.registry().purge_index(world,OrbitalIndex(L"a"));

        // add new virtual orbial
        vir_space = OrbitalSpaceTArray(OrbitalIndex(L"a"), OrbitalIndex(L"Α"), C_vir_ta_new);
        mo_int.orbital_space()->add(vir_space);

      auto mo_time1 = mpqc_time::fenced_now(world);
      auto mo_time = mpqc_time::duration_in_s(mo_time0,mo_time1);
      utility::print_par(world,"ClosedShell Dual Basis MO Build Time: ", mo_time, " S\n");

        return tre;
    }


template <typename Tile, typename Policy>
void closed_shell_dualbasis_cabs_mo_build_svd(
        integrals::MolecularIntegral <Tile, Policy> &mo_int,
        const rapidjson::Document &in,
        const std::shared_ptr<TRange1Engine> tre)
{
  auto& ao_int = mo_int.atomic_integral();
  auto orbital_registry = mo_int.orbital_space();
  auto &world = ao_int.world();
  // CABS fock build
  auto mo_time0 = mpqc_time::fenced_now(world);
  utility::print_par(world, "\nBuilding ClosedShell Dual Basis CABS MO Orbital\n");

  // integral
  auto S_ribs = ao_int.compute(L"<ρ|σ>");
  auto S_obs_ribs = ao_int.compute(L"<μ|σ>");
  auto S_vbs_ribs = ao_int.compute(L"<Α|σ>");

  // construct cabs
  MatrixD C_cabs_eigen;
  MatrixD C_allvir_eigen;
  std::size_t nbf_ribs_minus_occ;
  {
    MatrixD S_ribs_eigen = array_ops::array_to_eigen(S_ribs);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(S_ribs_eigen);
    MatrixD X_ribs_eigen_inv = es.operatorInverseSqrt();

    MatrixD S_obs_ribs_eigen = array_ops::array_to_eigen(S_obs_ribs);

    // C_mu^i
    TA::DistArray<Tile,Policy> Ci = orbital_registry->retrieve(OrbitalIndex(L"m")).array();
    MatrixD Ci_eigen = array_ops::array_to_eigen(Ci);

    MatrixD X1 = Ci_eigen.transpose() * S_obs_ribs_eigen * X_ribs_eigen_inv;

    // SVD solve to project out occupied
    Eigen::JacobiSVD<MatrixD> svd1(X1, Eigen::ComputeFullV);
    MatrixD V_eigen = svd1.matrixV();
    size_t nbf_ribs = S_ribs_eigen.cols();
    nbf_ribs_minus_occ = nbf_ribs - svd1.nonzeroSingularValues();
    MatrixD Vnull1(nbf_ribs, nbf_ribs_minus_occ);
    Vnull1 = V_eigen.block(0, svd1.nonzeroSingularValues(), nbf_ribs, nbf_ribs_minus_occ);
    C_allvir_eigen = X_ribs_eigen_inv * Vnull1;


    MatrixD S_vbs_ribs_eigen = array_ops::array_to_eigen(S_vbs_ribs);
    // C_a
    TA::DistArray<Tile,Policy> Ca = orbital_registry->retrieve(OrbitalIndex(L"a")).array();
    MatrixD Ca_eigen = array_ops::array_to_eigen(Ca);
    MatrixD X2 = Ca_eigen.transpose() * S_vbs_ribs_eigen * C_allvir_eigen;

    Eigen::JacobiSVD<MatrixD> svd2(X2, Eigen::ComputeFullV);
    MatrixD V_eigen2 = svd2.matrixV();
    auto nbf_cabs = nbf_ribs_minus_occ - svd2.nonzeroSingularValues();
    MatrixD Vnull2(nbf_ribs_minus_occ, nbf_cabs);
    Vnull2 = V_eigen2.block(0, svd2.nonzeroSingularValues(), nbf_ribs_minus_occ, nbf_cabs);
    C_cabs_eigen = C_allvir_eigen*Vnull2;
  }

 // get mo block size
 std::size_t mo_blocksize = in.HasMember("MoBlockSize") ? in["MoBlockSize"].GetInt() : 24;
 utility::print_par(world, "MoBlockSize: ", mo_blocksize, "\n");

  // get cabs trange
  auto tr_cabs = mo_int.atomic_integral().orbital_basis_registry()->retrieve(OrbitalIndex(L"α")).create_trange1();
  auto tr_ribs = S_ribs.trange().data().back();
  auto tr_cabs_mo = tre->compute_range(tr_cabs.elements().second, mo_blocksize);
  auto tr_allvir_mo = tre->compute_range(nbf_ribs_minus_occ, mo_blocksize);


 utility::parallel_print_range_info(world, tr_cabs_mo, "CABS MO");
  TA::DistArray<Tile,Policy>  C_cabs = array_ops::eigen_to_array<TA::TensorD>(world, C_cabs_eigen, tr_ribs, tr_cabs_mo);

    // insert to orbital space
 using OrbitalSpaceTArray = OrbitalSpace<TA::DistArray < Tile, Policy>>;
 auto C_cabs_space = OrbitalSpaceTArray(OrbitalIndex(L"a'"), OrbitalIndex(L"ρ"), C_cabs);
 orbital_registry->add(C_cabs_space);

  utility::parallel_print_range_info(world, tr_allvir_mo, "All Virtual MO");
  TA::DistArray<Tile,Policy>  C_allvir = array_ops::eigen_to_array<TA::TensorD>(world, C_allvir_eigen, tr_ribs, tr_allvir_mo);

  // insert to orbital space
  auto C_allvir_space = OrbitalSpaceTArray(OrbitalIndex(L"A'"), OrbitalIndex(L"ρ"), C_allvir);
  orbital_registry->add(C_allvir_space);

 auto mo_time1 = mpqc_time::fenced_now(world);
 auto mo_time = mpqc_time::duration_in_s(mo_time0, mo_time1);
 utility::print_par(world, "ClosedShell Dual Basis CABS MO Build Time: ", mo_time, " S\n");


};

} // end of namespace mpqc


#endif //MPQC_MO_BUILD_H
