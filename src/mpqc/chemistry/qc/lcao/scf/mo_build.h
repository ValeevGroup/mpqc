//
// Created by Chong Peng on 4/14/16.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_MO_BUILD_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_MO_BUILD_H_

#include <tiledarray.h>

#include "mpqc/chemistry/qc/lcao/expression/orbital_registry.h"
#include "mpqc/chemistry/qc/lcao/integrals/lcao_factory.h"
#include "mpqc/chemistry/qc/lcao/expression/trange1_engine.h"

namespace mpqc {
namespace lcao {

template <typename Tile, typename Policy>
std::shared_ptr<::mpqc::utility::TRange1Engine> closed_shell_obs_mo_build_eigen_solve(
    LCAOFactory<Tile, Policy> &lcao_factory, Eigen::VectorXd &ens, std::size_t nocc,
    const Molecule &mols, bool frozen_core, std::size_t occ_blocksize,
    std::size_t vir_blocksize) {
  auto &ao_factory = lcao_factory.ao_factory();
  auto &orbital_registry = lcao_factory.orbital_space();
  auto &world = ao_factory.world();
  using TArray = TA::DistArray<Tile, Policy>;

  auto mo_time0 = mpqc::fenced_now(world);
  utility::print_par(world, "\nBuilding ClosedShell OBS MO Orbital\n");

  // find fock matrix
  TArray F;
  if (ao_factory.registry().have(Formula(L"<μ|F|ν>"))) {
    F = ao_factory.compute(Formula(L"<μ|F|ν>"));
  } else {
    F = ao_factory.compute(Formula(L"<μ|F|ν>[df]"));
  }

  auto S = ao_factory.compute(L"<κ|λ>");

  RowMatrixXd F_eig = array_ops::array_to_eigen(F);
  RowMatrixXd S_eig = array_ops::array_to_eigen(S);

  // solve mo coefficients
  Eigen::GeneralizedSelfAdjointEigenSolver<RowMatrixXd> es(F_eig, S_eig);

  // start to solve coefficient

  std::size_t n_frozen_core = 0;
  if (frozen_core) {
    n_frozen_core = mols.core_electrons();
    utility::print_par(world, "Frozen Core: ", n_frozen_core, " electrons",
                       "\n");
    n_frozen_core = n_frozen_core / 2;
  }

  ens = es.eigenvalues();
  RowMatrixXd C_all = es.eigenvectors();
  RowMatrixXd C_occ = C_all.block(0, 0, S_eig.rows(), nocc);
  RowMatrixXd C_corr_occ =
      C_all.block(0, n_frozen_core, S_eig.rows(), nocc - n_frozen_core);
  RowMatrixXd C_vir = C_all.rightCols(S_eig.rows() - nocc);

  utility::print_par(world, "OccBlockSize: ", occ_blocksize, "\n");
  utility::print_par(world, "VirBlockSize: ", vir_blocksize, "\n");

  std::size_t all = S.trange().elements_range().extent()[0];
  using TRange1Engine = ::mpqc::utility::TRange1Engine;
  auto tre = std::make_shared<TRange1Engine>(nocc, all, occ_blocksize,
                                             vir_blocksize, n_frozen_core);

  // get all the trange1s
  auto tr_obs = S.trange().data().back();
  auto tr_corr_occ = tre->get_active_occ_tr1();
  auto tr_occ = tre->compute_range(nocc, occ_blocksize);
  auto tr_vir = tre->get_vir_tr1();
  auto tr_all = tre->get_all_tr1();

  mpqc::detail::parallel_print_range_info(world, tr_occ, "Occ");
  mpqc::detail::parallel_print_range_info(world, tr_corr_occ, "CorrOcc");
  mpqc::detail::parallel_print_range_info(world, tr_vir, "Vir");
  mpqc::detail::parallel_print_range_info(world, tr_all, "Obs");

  // convert to TA
  auto C_occ_ta = array_ops::eigen_to_array<Tile,Policy>(world, C_occ, tr_obs, tr_occ);
  auto C_corr_occ_ta =
      array_ops::eigen_to_array<Tile,Policy>(world, C_corr_occ, tr_obs, tr_corr_occ);
  auto C_vir_ta = array_ops::eigen_to_array<Tile,Policy>(world, C_vir, tr_obs, tr_vir);
  auto C_all_ta = array_ops::eigen_to_array<Tile,Policy>(world, C_all, tr_obs, tr_all);

  // insert to registry
  using OrbitalSpaceTArray = OrbitalSpace<TA::DistArray<Tile, Policy>>;
  auto occ_space =
      OrbitalSpaceTArray(OrbitalIndex(L"m"), OrbitalIndex(L"κ"), C_occ_ta);
  orbital_registry.add(occ_space);

  auto corr_occ_space =
      OrbitalSpaceTArray(OrbitalIndex(L"i"), OrbitalIndex(L"κ"), C_corr_occ_ta);
  orbital_registry.add(corr_occ_space);

  auto vir_space =
      OrbitalSpaceTArray(OrbitalIndex(L"a"), OrbitalIndex(L"κ"), C_vir_ta);
  orbital_registry.add(vir_space);

  auto obs_space =
      OrbitalSpaceTArray(OrbitalIndex(L"p"), OrbitalIndex(L"κ"), C_all_ta);
  orbital_registry.add(obs_space);

  auto mo_time1 = mpqc::fenced_now(world);
  auto mo_time = mpqc::duration_in_s(mo_time0, mo_time1);
  utility::print_par(world, "ClosedShell OBS MO Build Time: ", mo_time, " S\n");

  return tre;
}

template <typename Tile, typename Policy>
void closed_shell_cabs_mo_build_svd(
    LCAOFactory<Tile, Policy> &lcao_factory,
    const std::shared_ptr<::mpqc::utility::TRange1Engine> tre, std::size_t vir_blocksize) {
  auto &ao_factory = lcao_factory.ao_factory();
  auto &orbital_registry = lcao_factory.orbital_space();
  auto &world = ao_factory.world();
  // CABS fock build
  auto mo_time0 = mpqc::fenced_now(world);
  utility::print_par(world, "\nBuilding ClosedShell CABS MO Orbital\n");

  // build the RI basis

  const auto abs_basis =
      *ao_factory.orbital_basis_registry().retrieve(OrbitalIndex(L"α"));
  const auto obs_basis =
      *ao_factory.orbital_basis_registry().retrieve(OrbitalIndex(L"κ"));

  gaussian::Basis ri_basis;
  ri_basis = merge(obs_basis, abs_basis);

  mpqc::detail::parallel_print_range_info(world, ri_basis.create_trange1(),
                                    "RI Basis");
  ao_factory.orbital_basis_registry().add(OrbitalIndex(L"ρ"), std::make_shared<gaussian::Basis>(ri_basis));

  // integral
  auto S_ribs_inv = ao_factory.compute(L"<ρ|σ>[inv_sqr]");
  auto S_obs_ribs = ao_factory.compute(L"<μ|σ>");
  auto S_obs_inv = ao_factory.compute(L"<κ|λ>[inv_sqr]");

  // construct cabs
  TA::DistArray<Tile, Policy> C_cabs, C_ri, C_allvir;
  {

    // orthogonalize
    decltype(S_obs_inv) S_obs_ribs_ortho;
    S_obs_ribs_ortho("i,j") = S_obs_inv("i,k") * S_obs_ribs("k,l")*S_ribs_inv("l,j");
    RowMatrixXd S_obs_ribs_ortho_eigen = array_ops::array_to_eigen(S_obs_ribs_ortho);

    // SVD solve
    Eigen::JacobiSVD<RowMatrixXd> svd(S_obs_ribs_ortho_eigen,
                                      Eigen::ComputeFullV);
    RowMatrixXd V_eigen = svd.matrixV();
    size_t nbf_ribs = S_obs_ribs_ortho_eigen.cols();
    auto nbf_cabs = nbf_ribs - svd.nonzeroSingularValues();
    RowMatrixXd Vnull(nbf_ribs, nbf_cabs);
    Vnull = V_eigen.block(0, svd.nonzeroSingularValues(), nbf_ribs, nbf_cabs);

    auto tr_ribs = ri_basis.create_trange1();
    auto tr_cabs_mo = tre->compute_range(nbf_cabs, vir_blocksize);
    mpqc::detail::parallel_print_range_info(world, tr_cabs_mo, "CABS MO");

    C_cabs = array_ops::eigen_to_array<Tile,Policy>(world, Vnull, tr_ribs, tr_cabs_mo);
    C_cabs("i,j") = S_ribs_inv("i,k") * C_cabs("k, j");

    RowMatrixXd C_cabs_eigen = array_ops::array_to_eigen(C_cabs);

    // solve orbitals for all virtual

    auto n_vir = tre->get_vir();
    RowMatrixXd C_allvirtual_eigen =
        RowMatrixXd::Zero(nbf_ribs, n_vir + nbf_cabs);

    {
      auto C_vir = orbital_registry.retrieve(OrbitalIndex(L"a")).coefs();
      RowMatrixXd C_vir_eigen = array_ops::array_to_eigen(C_vir);

      auto n_obs = C_vir_eigen.rows();

      C_allvirtual_eigen.block(0, 0, n_obs, n_vir) << C_vir_eigen;
      C_allvirtual_eigen.block(0, n_vir, nbf_ribs, nbf_cabs) << C_cabs_eigen;
    }


    // reblock C_ribs
    auto tr_ribs_mo = tre->compute_range(nbf_ribs, vir_blocksize);
    mpqc::detail::parallel_print_range_info(world, tr_ribs_mo, "RIBS MO");
    auto ribs_to_mo = array_ops::create_diagonal_array_from_eigen<Tile, Policy>(world, tr_ribs, tr_ribs_mo, 1.0);
    C_ri("i,j") = S_ribs_inv("i,k")*ribs_to_mo("k,j");


    auto tr_allvir_mo = tre->compute_range(nbf_cabs + n_vir, vir_blocksize);
    mpqc::detail::parallel_print_range_info(world, tr_allvir_mo, "All Virtual MO");

    C_allvir = array_ops::eigen_to_array<Tile,Policy>(world, C_allvirtual_eigen,
                                                      tr_ribs, tr_allvir_mo);

    // insert to orbital space
    using OrbitalSpaceTArray = OrbitalSpace<TA::DistArray<Tile, Policy>>;
    auto C_cabs_space =
        OrbitalSpaceTArray(OrbitalIndex(L"a'"), OrbitalIndex(L"ρ"), C_cabs);
    auto C_ribs_space =
        OrbitalSpaceTArray(OrbitalIndex(L"P'"), OrbitalIndex(L"ρ"), C_ri);
    auto C_allvir_space =
        OrbitalSpaceTArray(OrbitalIndex(L"A'"), OrbitalIndex(L"ρ"), C_allvir);

    orbital_registry.add(C_cabs_space);
    orbital_registry.add(C_ribs_space);
    orbital_registry.add(C_allvir_space);

    auto mo_time1 = mpqc::fenced_now(world);
    auto mo_time = mpqc::duration_in_s(mo_time0, mo_time1);
    utility::print_par(world, "ClosedShell CABS MO Build Time: ", mo_time,
                       " S\n");
  }
};

template <typename Tile, typename Policy>
std::shared_ptr<::mpqc::utility::TRange1Engine> closed_shell_dualbasis_mo_build_eigen_solve_svd(
    LCAOFactory<Tile, Policy> &lcao_factory, Eigen::VectorXd &ens, std::size_t nocc,
    const Molecule &mols, bool frozen_core, std::size_t occ_blocksize,
    std::size_t vir_blocksize) {
  auto &ao_factory = lcao_factory.ao_factory();
  auto &world = ao_factory.world();
  using TArray = TA::DistArray<Tile, Policy>;

  utility::print_par(world, "\nBuilding ClosedShell Dual Basis MO Orbital\n");
  auto mo_time0 = mpqc::fenced_now(world);

  // solving occupied orbitals
  TArray F;
  if (ao_factory.registry().have(Formula(L"<μ|F|ν>"))) {
    F = ao_factory.compute(Formula(L"<μ|F|ν>"));
  } else {
    F = ao_factory.compute(Formula(L"<μ|F|ν>[df]"));
  }

  auto S = ao_factory.compute(L"<κ|λ>");

  RowMatrixXd F_eig = array_ops::array_to_eigen(F);
  RowMatrixXd S_eig = array_ops::array_to_eigen(S);

  // solve mo coefficients
  Eigen::GeneralizedSelfAdjointEigenSolver<RowMatrixXd> es(F_eig, S_eig);

  std::size_t n_frozen_core = 0;
  if (frozen_core) {
    n_frozen_core = mols.core_electrons();
    utility::print_par(world, "Frozen Core: ", n_frozen_core, " electrons",
                       "\n");
    n_frozen_core = n_frozen_core / 2;
  }

  Eigen::VectorXd ens_occ = es.eigenvalues().segment(0, nocc);
//  std::cout << "Energy of Occupied: \n" << ens_occ << std::endl;
  RowMatrixXd C_all = es.eigenvectors();
  RowMatrixXd C_occ = C_all.block(0, 0, S_eig.rows(), nocc);
  RowMatrixXd C_corr_occ =
      C_all.block(0, n_frozen_core, S_eig.rows(), nocc - n_frozen_core);

  // finished solving occupied orbitals

  // start to solve virtual orbitals
  auto S_vbs = ao_factory.compute(L"<Α|Β>");
  auto S_obs_vbs = ao_factory.compute(L"<μ|Α>");

  // construct C_vbs
  RowMatrixXd C_vbs;
  std::size_t nbf_vbs;
  std::size_t nbf_v;
  {
    // S_A^B -(1/2)
    RowMatrixXd S_vbs_eigen = array_ops::array_to_eigen(S_vbs);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es2(S_vbs_eigen);
    RowMatrixXd X_vbs_eigen_inv = es2.operatorInverseSqrt();

    // orthogonalize
    RowMatrixXd S_obs_vbs_eigen = array_ops::array_to_eigen(S_obs_vbs);
    RowMatrixXd S_obs_vbs_ortho_eigen = S_obs_vbs_eigen * X_vbs_eigen_inv;

    RowMatrixXd X_i_mu = C_occ.transpose() * S_obs_vbs_ortho_eigen;

    // SVD solve
    Eigen::JacobiSVD<RowMatrixXd> svd(X_i_mu, Eigen::ComputeFullV);
    RowMatrixXd V_eigen = svd.matrixV();
    nbf_vbs = S_obs_vbs_ortho_eigen.cols();
    nbf_v = nbf_vbs - svd.nonzeroSingularValues();
    std::cout << "Virtual Orbital Size:  " << nbf_v << std::endl;
    RowMatrixXd Vnull(nbf_vbs, nbf_v);
    Vnull = V_eigen.block(0, svd.nonzeroSingularValues(), nbf_vbs, nbf_v);
    C_vbs = X_vbs_eigen_inv.transpose() * Vnull;
  }

  utility::print_par(world, "OccBlockSize: ", occ_blocksize, "\n");
  utility::print_par(world, "VirBlockSize: ", vir_blocksize, "\n");

  using TRange1Engine = ::mpqc::utility::TRange1Engine;
  auto tre = std::make_shared<TRange1Engine>(nocc, nbf_vbs, occ_blocksize,
                                             vir_blocksize, n_frozen_core);
  auto tr_obs = S.trange().data().back();
  auto tr_vbs = S_vbs.trange().data().back();
  auto tr_occ = tre->compute_range(nocc, occ_blocksize);
  auto tr_corr_occ = tre->get_active_occ_tr1();
  auto tr_vir = tre->get_vir_tr1();

  mpqc::detail::parallel_print_range_info(world, tr_occ, "Occ");
  mpqc::detail::parallel_print_range_info(world, tr_corr_occ, "CorrOcc");
  mpqc::detail::parallel_print_range_info(world, tr_vir, "Vir");

  // convert to TA
  auto C_occ_ta = array_ops::eigen_to_array<Tile,Policy>(world, C_occ, tr_obs, tr_occ);
  auto C_corr_occ_ta =
      array_ops::eigen_to_array<Tile,Policy>(world, C_corr_occ, tr_obs, tr_corr_occ);
  auto C_vir_ta = array_ops::eigen_to_array<Tile,Policy>(world, C_vbs, tr_vbs, tr_vir);

  // insert to registry
  using OrbitalSpaceTArray = OrbitalSpace<TA::DistArray<Tile, Policy>>;
  auto occ_space =
      OrbitalSpaceTArray(OrbitalIndex(L"m"), OrbitalIndex(L"κ"), C_occ_ta);
  lcao_factory.orbital_space().add(occ_space);

  auto corr_occ_space =
      OrbitalSpaceTArray(OrbitalIndex(L"i"), OrbitalIndex(L"κ"), C_corr_occ_ta);
  lcao_factory.orbital_space().add(corr_occ_space);

  auto vir_space =
      OrbitalSpaceTArray(OrbitalIndex(L"a"), OrbitalIndex(L"Α"), C_vir_ta);
  lcao_factory.orbital_space().add(vir_space);

  // solve energy in virtual orbital
  TArray F_vbs;
  if (ao_factory.registry().have(Formula(L"<μ|F|ν>"))) {
    F_vbs = lcao_factory.compute(L"<a|F|b>");
  } else {
    F_vbs = lcao_factory.compute(L"<a|F|b>[df]");
  }
  RowMatrixXd F_vbs_mo_eigen = array_ops::array_to_eigen(F_vbs);
  //    std::cout << "F_vbs MO" << std::endl;
  //    std::cout << F_vbs_mo_eigen << std::endl;

  Eigen::SelfAdjointEigenSolver<RowMatrixXd> es3(F_vbs_mo_eigen);
  auto ens_vir = es3.eigenvalues();

  ens = Eigen::VectorXd(nbf_vbs);
  ens << ens_occ, ens_vir;

  std::cout << "Energy of Orbitals " << std::endl;
  std::cout << ens << std::endl;

  // resolve the virtual orbitals
  RowMatrixXd C_vir_rotate = es3.eigenvectors();
  C_vbs = C_vbs * C_vir_rotate;
  TArray C_vir_ta_new =
      array_ops::eigen_to_array<Tile,Policy>(world, C_vbs, tr_vbs, tr_vir);

  // remove old virtual orbitals
  lcao_factory.orbital_space().remove(OrbitalIndex(L"a"));
  lcao_factory.registry().purge_index(world, OrbitalIndex(L"a"));

  // add new virtual orbial
  vir_space =
      OrbitalSpaceTArray(OrbitalIndex(L"a"), OrbitalIndex(L"Α"), C_vir_ta_new);
  lcao_factory.orbital_space().add(vir_space);

  auto mo_time1 = mpqc::fenced_now(world);
  auto mo_time = mpqc::duration_in_s(mo_time0, mo_time1);
  utility::print_par(world, "ClosedShell Dual Basis MO Build Time: ", mo_time,
                     " S\n");

  return tre;
}

template <typename Tile, typename Policy>
void closed_shell_dualbasis_cabs_mo_build_svd(
    LCAOFactory<Tile, Policy> &lcao_factory,
    const std::shared_ptr<::mpqc::utility::TRange1Engine> tre, std::string ri_method,
    std::size_t vir_blocksize) {
  auto &ao_factory = lcao_factory.ao_factory();
  auto &orbital_registry = lcao_factory.orbital_space();
  auto &world = ao_factory.world();
  // CABS fock build
  auto mo_time0 = mpqc::fenced_now(world);
  utility::print_par(world,
                     "\nBuilding ClosedShell Dual Basis CABS MO Orbital\n");

  // build RI Basis First
  const auto abs_basis =
      *ao_factory.orbital_basis_registry().retrieve(OrbitalIndex(L"α"));
  const auto vir_basis =
      *ao_factory.orbital_basis_registry().retrieve(OrbitalIndex(L"Α"));
  const auto obs_basis =
      *ao_factory.orbital_basis_registry().retrieve(OrbitalIndex(L"κ"));

  std::shared_ptr<gaussian::Basis> ri_basis;

  if (ri_method == "VBS") {
    *ri_basis = merge(vir_basis, abs_basis);
    mpqc::detail::parallel_print_range_info(world, ri_basis->create_trange1(),
                                      "RI Basis with VBS");
  } else if (ri_method == "OBS") {
    *ri_basis = merge(obs_basis, abs_basis);
    mpqc::detail::parallel_print_range_info(world, ri_basis->create_trange1(),
                                      "RI Basis with OBS");
  } else {
    throw std::runtime_error("Invalid RI Method!");
  }
  ao_factory.orbital_basis_registry().add(OrbitalIndex(L"ρ"), ri_basis);

  // integral
  auto S_ribs = ao_factory.compute(L"<ρ|σ>");
  auto S_obs_ribs = ao_factory.compute(L"<μ|σ>");
  auto S_vbs_ribs = ao_factory.compute(L"<Α|σ>");

  // construct cabs
  RowMatrixXd C_cabs_eigen;
  RowMatrixXd C_allvir_eigen;
  std::size_t nbf_ribs_minus_occ;
  {
    RowMatrixXd S_ribs_eigen = array_ops::array_to_eigen(S_ribs);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(S_ribs_eigen);
    RowMatrixXd X_ribs_eigen_inv = es.operatorInverseSqrt();

    RowMatrixXd S_obs_ribs_eigen = array_ops::array_to_eigen(S_obs_ribs);

    // C_mu^i
    TA::DistArray<Tile, Policy> Ci =
        orbital_registry.retrieve(OrbitalIndex(L"m")).coefs();
    RowMatrixXd Ci_eigen = array_ops::array_to_eigen(Ci);

    RowMatrixXd X1 = Ci_eigen.transpose() * S_obs_ribs_eigen * X_ribs_eigen_inv;

    // SVD solve to project out occupied
    Eigen::JacobiSVD<RowMatrixXd> svd1(X1, Eigen::ComputeFullV);
    RowMatrixXd V_eigen = svd1.matrixV();
    size_t nbf_ribs = S_ribs_eigen.cols();
    nbf_ribs_minus_occ = nbf_ribs - svd1.nonzeroSingularValues();
    RowMatrixXd Vnull1(nbf_ribs, nbf_ribs_minus_occ);
    Vnull1 = V_eigen.block(0, svd1.nonzeroSingularValues(), nbf_ribs,
                           nbf_ribs_minus_occ);
    C_allvir_eigen = X_ribs_eigen_inv * Vnull1;

    RowMatrixXd S_vbs_ribs_eigen = array_ops::array_to_eigen(S_vbs_ribs);
    // C_a
    TA::DistArray<Tile, Policy> Ca =
        orbital_registry.retrieve(OrbitalIndex(L"a")).coefs();
    RowMatrixXd Ca_eigen = array_ops::array_to_eigen(Ca);
    RowMatrixXd X2 = Ca_eigen.transpose() * S_vbs_ribs_eigen * C_allvir_eigen;

    Eigen::JacobiSVD<RowMatrixXd> svd2(X2, Eigen::ComputeFullV);
    RowMatrixXd V_eigen2 = svd2.matrixV();
    auto nbf_cabs = nbf_ribs_minus_occ - svd2.nonzeroSingularValues();
    RowMatrixXd Vnull2(nbf_ribs_minus_occ, nbf_cabs);
    Vnull2 = V_eigen2.block(0, svd2.nonzeroSingularValues(), nbf_ribs_minus_occ,
                            nbf_cabs);
    C_cabs_eigen = C_allvir_eigen * Vnull2;
  }

  utility::print_par(world, "VirBlockSize: ", vir_blocksize, "\n");
  // get cabs trange
  auto tr_cabs = lcao_factory.ao_factory()
                     .orbital_basis_registry()
                     .retrieve(OrbitalIndex(L"α"))
                     ->create_trange1();
  auto tr_ribs = S_ribs.trange().data().back();
  auto tr_cabs_mo =
      tre->compute_range(tr_cabs.elements_range().second, vir_blocksize);
  auto tr_allvir_mo = tre->compute_range(nbf_ribs_minus_occ, vir_blocksize);

  mpqc::detail::parallel_print_range_info(world, tr_cabs_mo, "CABS MO");
  TA::DistArray<Tile, Policy> C_cabs = array_ops::eigen_to_array<Tile,Policy>(
      world, C_cabs_eigen, tr_ribs, tr_cabs_mo);

  // insert to orbital space
  using OrbitalSpaceTArray = OrbitalSpace<TA::DistArray<Tile, Policy>>;
  auto C_cabs_space =
      OrbitalSpaceTArray(OrbitalIndex(L"a'"), OrbitalIndex(L"ρ"), C_cabs);
  orbital_registry.add(C_cabs_space);

  mpqc::detail::parallel_print_range_info(world, tr_allvir_mo, "All Virtual MO");
  TA::DistArray<Tile, Policy> C_allvir = array_ops::eigen_to_array<Tile,Policy>(
      world, C_allvir_eigen, tr_ribs, tr_allvir_mo);

  // insert to orbital space
  auto C_allvir_space =
      OrbitalSpaceTArray(OrbitalIndex(L"A'"), OrbitalIndex(L"ρ"), C_allvir);
  orbital_registry.add(C_allvir_space);

  auto mo_time1 = mpqc::fenced_now(world);
  auto mo_time = mpqc::duration_in_s(mo_time0, mo_time1);
  utility::print_par(world, "ClosedShell Dual Basis CABS MO Build Time: ",
                     mo_time, " S\n");
};

}  // namespace lcao
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_MO_BUILD_H_
