//
// Created by Chong Peng on 4/14/16.
//

#ifndef MPQC_MO_BUILD_H
#define MPQC_MO_BUILD_H


#include <rapidjson/document.h>

#include "../include/tiledarray.h"
#include "../common/namespaces.h"
#include "../integrals/atomic_integral.h"
#include "../expression/orbital_registry.h"
#include "../utility/trange1_engine.h"

namespace mpqc {


template <typename Tile, typename Policy>
std::shared_ptr<TRange1Engine> closed_shell_obs_mo_build_eigen_solve(
        integrals::AtomicIntegral<Tile, Policy> &ao_int,
        OrbitalSpaceRegistry<TA::DistArray < Tile, Policy>>& orbital_registry,
        Eigen::VectorXd& ens,
        const rapidjson::Document& in,
        const molecule::Molecule& mols,
        int occ)
{
    auto& world = ao_int.get_world();
    using TArray = TA::DistArray<Tile, Policy>;


    auto mo_time0 = mpqc_time::fenced_now(world);
    utility::print_par(world, "\nBuilding ClosedShell OBS MO Orbital\n");

    // find fock matrix
    TArray F;
    if(ao_int.registry().have(Formula(L"(μ|F|ν)"))){
        F = ao_int.registry().retrieve(Formula(L"(μ|F|ν)"));
    }
    else{
        F = ao_int.registry().retrieve(Formula(L"(μ|F|ν)[df]"));
    }

    auto S = ao_int.compute(L"(κ|λ)");

    auto F_eig = array_ops::array_to_eigen(F);
    auto S_eig = array_ops::array_to_eigen(S);

    // check the condition number in Overlap
    Eig::SelfAdjointEigenSolver<decltype(S_eig)> S_es(S_eig);
    // eigen value in increasing order
    auto cond = S_es.eigenvalues()(S_es.eigenvalues().size() - 1) / S_es.eigenvalues()(0);
    utility::print_par(world, "Condition Number in Overlap: ", cond, "\n");


    // solve mo coefficients
    Eig::GeneralizedSelfAdjointEigenSolver<decltype(S_eig)> es(F_eig, S_eig);

    // start to solve coefficient
    bool frozen_core = in.HasMember("FrozenCore") ? in["FrozenCore"].GetBool() : false;
    std::size_t n_frozen_core = 0;
    if (frozen_core) {
        n_frozen_core = mols.core_electrons();
        utility::print_par(world, "Frozen Core: ", n_frozen_core," electrons", "\n");
        n_frozen_core = n_frozen_core / 2;
    }

    ens = es.eigenvalues().bottomRows(S_eig.rows() - n_frozen_core);
    Eig::MatrixXd C_all = es.eigenvectors();
    Eig::MatrixXd C_occ = C_all.block(0, 0, S_eig.rows(),occ);
    Eig::MatrixXd C_corr_occ = C_all.block(0, n_frozen_core, S_eig.rows(), occ - n_frozen_core);
    Eig::MatrixXd C_vir = C_all.rightCols(S_eig.rows() - occ);

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
    using OrbitalSpace = OrbitalSpace<decltype(C_occ_ta)>;
    auto occ_space = OrbitalSpace(OrbitalIndex(L"m"), OrbitalIndex(L"κ"), C_occ_ta);
    orbital_registry.add(occ_space);

    auto corr_occ_space = OrbitalSpace(OrbitalIndex(L"i"), OrbitalIndex(L"κ"), C_corr_occ_ta);
    orbital_registry.add(corr_occ_space);

    auto vir_space = OrbitalSpace(OrbitalIndex(L"a"), OrbitalIndex(L"κ"), C_vir_ta);
    orbital_registry.add(vir_space);

    auto obs_space = OrbitalSpace(OrbitalIndex(L"p"),OrbitalIndex(L"κ"), C_all_ta);
    orbital_registry.add(obs_space);

    auto mo_time1 = mpqc_time::fenced_now(world);
    auto mo_time = mpqc_time::duration_in_s(mo_time0,mo_time1);
    utility::print_par(world,"ClosedShell OBS MO Build Time: ", mo_time, " S\n");

    return tre;
};


template <typename Tile, typename Policy>
void closed_shell_cabs_mo_build_eigen_solve(
        integrals::AtomicIntegral<Tile, Policy> &ao_int,
        OrbitalSpaceRegistry<TA::DistArray < Tile, Policy>>& orbital_registry,
        const rapidjson::Document& in,
        const std::shared_ptr<TRange1Engine> tre)
{
    auto& world = ao_int.get_world();
    // CABS fock build
    auto mo_time0 = mpqc_time::fenced_now(world);
    utility::print_par(world, "\nBuilding ClosedShell CABS MO Orbital\n");

    // integral
    auto S_cabs = ao_int.compute(L"(α|β)");
    auto S_ribs = ao_int.compute(L"(ρ|σ)");
    auto S_obs_ribs = ao_int.compute(L"(μ|σ)");
    auto S_obs = ao_int.compute(L"(κ|λ)");

    // construct cabs
    decltype(S_obs) C_cabs, C_ri;
    {
        auto S_obs_eigen = array_ops::array_to_eigen(S_obs);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(S_obs_eigen);
        MatrixD X_obs_eigen_inv = es.operatorInverseSqrt();

        auto S_ribs_eigen = array_ops::array_to_eigen(S_ribs);
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

        // get mo block size
        std::size_t mo_blocksize = in.HasMember("MoBlockSize") ? in["MoBlockSize"].GetInt() : 24;
        utility::print_par(world,"MoBlockSize: ", mo_blocksize, "\n");

        auto tr_ribs = S_ribs.trange().data()[0];
        auto tr_cabs = S_cabs.trange().data()[0];
        auto tr_cabs_mo = tre->compute_range(tr_cabs.elements().second, mo_blocksize);
        auto tr_ribs_mo = tre->compute_range(tr_ribs.elements().second, mo_blocksize);


        utility::parallel_print_range_info(world, tr_cabs_mo, "CABS MO");
        utility::parallel_print_range_info(world, tr_ribs_mo, "RIBS MO");

        C_cabs = array_ops::eigen_to_array<TA::TensorD>(world, C_cabs_eigen, tr_ribs, tr_cabs_mo);
        C_ri = array_ops::eigen_to_array<TA::TensorD>(world, X_ribs_eigen_inv, tr_ribs, tr_ribs_mo);

        // insert to orbital space
        using OrbitalSpace = OrbitalSpace<decltype(C_cabs)>;
        auto C_cabs_space = OrbitalSpace(OrbitalIndex(L"a'"), OrbitalIndex(L"ρ"), C_cabs);
        auto C_ribs_space = OrbitalSpace(OrbitalIndex(L"P'"), OrbitalIndex(L"ρ"), C_ri);

        orbital_registry.add(C_cabs_space);
        orbital_registry.add(C_ribs_space);

        auto mo_time1 = mpqc_time::fenced_now(world);
        auto mo_time = mpqc_time::duration_in_s(mo_time0,mo_time1);
        utility::print_par(world,"ClosedShell CABS MO Build Time: ", mo_time, " S\n");
    }

};

} // end of namespace mpqc


#endif //MPQC_MO_BUILD_H
