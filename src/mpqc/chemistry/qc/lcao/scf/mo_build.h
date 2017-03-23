//
// Created by Chong Peng on 4/14/16.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_MO_BUILD_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_MO_BUILD_H_

#include <tiledarray.h>

#include "mpqc/chemistry/qc/lcao/expression/trange1_engine.h"
#include "mpqc/chemistry/qc/lcao/factory/lcao_factory.h"
#include "mpqc/chemistry/qc/lcao/factory/periodic_lcao_factory.h"
#include "mpqc/util/misc/provider.h"

namespace mpqc {
namespace lcao {

/// computes the MO-basis Fock matrix and extracts the diagonal elements
template <typename Tile, typename Policy>
std::shared_ptr<Eigen::VectorXd> make_diagonal_fpq(
    LCAOFactoryBase<Tile, Policy> &lcao_factory,
    gaussian::AOFactoryBase<Tile, Policy> &ao_factory) {
  bool df = ao_factory.registry().have(Formula(L"<μ|F|ν>[df]"));
  auto str = df ? L"<p|F|q>[df]" : L"<p|F|q>";
  auto Fpq_eig = array_ops::array_to_eigen(lcao_factory.compute(str));
  return std::make_shared<Eigen::VectorXd>(Fpq_eig.diagonal());
}

/// @brief populates the LCAOFactory object with
///        the canonical single-determinant closed-shell reference state orbital
///        spaces

/// Populates the LCAOFactory object with the occupied, active occupied,
/// unoccupied, and full orbtial spaces.
/// @param lcao the LCAOFactory object to be populated
/// @param p_space the CanonicalOrbitalSpace object for the full (occupied and
/// unoccupied) space
/// @param ndocc the number of doubly occupied orbitals in the reference
/// @param n_frozen_core the number of frozen core orbitals
/// @param occ_blksize the target block size for the occupied orbitals
/// @param unocc_blksize the target block size for the unoccupied orbitals
template <typename Tile, typename Policy>
void make_closed_shell_canonical_sdref_subspaces(
    std::shared_ptr<LCAOFactoryBase<Tile, Policy>> lcao_factory,
    std::shared_ptr<const CanonicalOrbitalSpace<TA::DistArray<Tile, Policy>>>
        p_space,
    std::size_t ndocc, std::size_t n_frozen_core, std::size_t occ_blksize,
    std::size_t unocc_blksize) {
  using TArray = TA::DistArray<Tile, Policy>;
  using COrbSpace = CanonicalOrbitalSpace<TArray>;

  const auto &eps_p = p_space->attributes();

  auto &orbital_registry = lcao_factory->orbital_registry();
  auto &world = lcao_factory->world();

  // divide the LCAO space into subspaces using Eigen .. boo
  RowMatrixXd C_p = array_ops::array_to_eigen(p_space->coefs());
  const auto n = C_p.cols();
  const auto nao = C_p.rows();
  RowMatrixXd C_m = C_p.block(0, 0, C_p.rows(), ndocc);
  RowMatrixXd C_i = C_m.block(0, n_frozen_core, nao, ndocc - n_frozen_core);
  RowMatrixXd C_a = C_p.rightCols(n - ndocc);

  using TRange1Engine = ::mpqc::utility::TRange1Engine;
  auto tre = std::make_shared<TRange1Engine>(ndocc, n, occ_blksize,
                                             unocc_blksize, n_frozen_core);

  // get all the trange1s
  auto tr_ao = p_space->coefs().trange().data()[0];
  auto tr_i = tre->get_active_occ_tr1();
  auto tr_m = utility::compute_trange1(ndocc, occ_blksize);
  auto tr_a = tre->get_vir_tr1();
  auto tr_p = tre->get_all_tr1();

  mpqc::detail::parallel_print_range_info(world, tr_m, "Occ Range");
  mpqc::detail::parallel_print_range_info(world, tr_i, "ActiveOcc Range");
  mpqc::detail::parallel_print_range_info(world, tr_a, "Unocc Range");
  mpqc::detail::parallel_print_range_info(world, tr_p, "Obs Range");

  // convert eigen arrays to TA
  auto C_m_ta =
      array_ops::eigen_to_array<Tile, Policy>(world, C_m, tr_ao, tr_m);
  auto C_i_ta =
      array_ops::eigen_to_array<Tile, Policy>(world, C_i, tr_ao, tr_i);
  auto C_a_ta =
      array_ops::eigen_to_array<Tile, Policy>(world, C_a, tr_ao, tr_a);

  // create orbital spaces and push into registry
  std::vector<double> eps_m(eps_p.begin(), eps_p.begin() + ndocc);
  auto m_space =
      COrbSpace(OrbitalIndex(L"m"), OrbitalIndex(L"κ"), C_m_ta, eps_m);
  orbital_registry.add(m_space);

  std::vector<double> eps_i(eps_m.begin() + n_frozen_core, eps_m.end());
  auto i_space =
      COrbSpace(OrbitalIndex(L"i"), OrbitalIndex(L"κ"), C_i_ta, eps_i);
  orbital_registry.add(i_space);

  std::vector<double> eps_a(eps_p.begin() + ndocc, eps_p.end());
  auto a_space =
      COrbSpace(OrbitalIndex(L"a"), OrbitalIndex(L"κ"), C_a_ta, eps_a);
  orbital_registry.add(a_space);

  // reblock the full space and add to the regisry
  auto C_p_ta =
      array_ops::eigen_to_array<Tile, Policy>(world, C_p, tr_ao, tr_p);
  auto p_space_reblocked =
      COrbSpace(OrbitalIndex(L"p"), OrbitalIndex(L"κ"), C_p_ta, eps_p);
  orbital_registry.add(p_space_reblocked);
}

/// @brief populates the LCAOFactory object with
///        the single-determinant closed-shell reference state orbital spaces.

/// Populates the LCAOFactory object with the occupied, active occupied,
/// unoccupied, and full orbital spaces. The input state is either the full
/// populated set of reference orbitals, or just the occupied subset. In the
/// latter
/// case the unoccupied orbitals are formed by svd-based projection of the
/// occupied
/// subspace of the occupied orbitals.
/// @param lcao the LCAOFactory object to be populated
/// @param wfn provides CanonicalOrbitalSpace, will be computed
/// @param target_ref_precision the precision with which to compute the
///                orbitals (also see
///                CanonicalOrbitalSpace::Provider::evaluate())
/// @param ndocc the number of doubly occupied orbitals in the reference
/// @param n_frozen_core the number of frozen core orbitals
/// @param occ_blksize the target block size for the occupied orbitals
/// @param unocc_blksize the target block size for the unoccupied orbitals
template <typename Tile, typename Policy>
void make_closed_shell_sdref_subspaces(
    std::shared_ptr<gaussian::AOFactoryBase<Tile, Policy>> ao_factory,
    std::shared_ptr<
        typename PopulatedOrbitalSpace<TA::DistArray<Tile, Policy>>::Provider>
        wfn,
    double target_ref_precision, std::size_t ndocc, std::size_t n_frozen_core,
    std::size_t occ_blksize, std::size_t unocc_blksize) {
  using TArray = TA::DistArray<Tile, Policy>;
  using POrbSpace = PopulatedOrbitalSpace<TArray>;
  using TRange1Engine = ::mpqc::utility::TRange1Engine;

  // compute a populated space
  std::shared_ptr<POrbSpace> input_space = std::make_shared<POrbSpace>();
  evaluate(*input_space, wfn, target_ref_precision);

  // receive either occupied or all orbitals
  assert(input_space->index() == OrbitalIndex(L"m") ||
         input_space->index() == OrbitalIndex(L"p"));

  auto &orbital_registry = ao_factory->orbital_registry();
  auto &world = ao_factory->world();

  if (input_space->rank() == ndocc) {
    assert(input_space->index() ==
           OrbitalIndex(
               L"m"));  // make sure ndocc is consistent with the occupied space

    const auto &occ_m = input_space->attributes();
    for (auto occ : occ_m) assert(occ == 2.0);  // all should be doubly occupied
    auto nao = input_space->ao_rank();
    auto tr_ao = input_space->ao_trange();

    //////////////////////////////////////////////////////////////////////////////////
    // rebuild occupied orbital space (in case block size differs)
    RowMatrixXd C_m_eig = array_ops::array_to_eigen(input_space->coefs());
    auto C_m = array_ops::eigen_to_array<Tile, Policy>(
        world, C_m_eig, tr_ao,
        utility::compute_trange1(C_m_eig.cols(), occ_blksize));
    auto m_space =
        POrbSpace(OrbitalIndex(L"m"), OrbitalIndex(L"κ"), C_m, occ_m);
    orbital_registry.add(m_space);

    //////////////////////////////////////////////////////////////////////////////////
    // build active occupied space
    RowMatrixXd C_i_eig =
        C_m_eig.block(0, n_frozen_core, nao, ndocc - n_frozen_core);
    auto C_i = array_ops::eigen_to_array<Tile, Policy>(
        world, C_i_eig, tr_ao,
        utility::compute_trange1(C_i_eig.cols(), occ_blksize));
    std::vector<double> occ_i(occ_m.begin() + n_frozen_core, occ_m.end());
    auto i_space =
        POrbSpace(OrbitalIndex(L"i"), OrbitalIndex(L"κ"), C_m, occ_i);
    orbital_registry.add(i_space);

    //////////////////////////////////////////////////////////////////////////////////
    // obtain unoccupieds by projecting occupieds from OBS
    // should really make PAOs here, but these should be close, although much
    // more expensive
    auto obs_basis = ao_factory->basis_registry()->retrieve(OrbitalIndex(L"κ"));

    // need some integrals
    auto S_m_obs = ao_factory->compute(L"<m|λ>");
    auto S_obs_inv = ao_factory->compute(L"<κ|λ>[inv_sqr]");

    decltype(S_m_obs) S_m_obs_ortho;
    S_m_obs_ortho("i,k") = S_m_obs("i,l") * S_obs_inv("l,k");
    RowMatrixXd S_m_obs_ortho_eig = array_ops::array_to_eigen(S_m_obs_ortho);

    // SVD to obtain the null-space basis of <m|kappa>, i.e. the unoccupied
    // orbitals
    Eigen::JacobiSVD<RowMatrixXd> svd(S_m_obs_ortho_eig, Eigen::ComputeFullV);
    RowMatrixXd V_eig = svd.matrixV();
    size_t nbf = S_m_obs_ortho_eig.cols();
    auto n_unocc = nbf - svd.nonzeroSingularValues();
    RowMatrixXd Vnull(nbf, n_unocc);
    Vnull = V_eig.block(0, svd.nonzeroSingularValues(), nbf, n_unocc);

    auto C_a = array_ops::eigen_to_array<Tile, Policy>(
        world, Vnull, obs_basis->create_trange1(),
        utility::compute_trange1(n_unocc, unocc_blksize));
    C_a("i,j") = S_obs_inv("i,k") * C_a("k, j");
    std::vector<double> occ_a(n_unocc, 0);
    auto a_space =
        POrbSpace(OrbitalIndex(L"a"), OrbitalIndex(L"κ"), C_a, occ_a);
    orbital_registry.add(a_space);

    //////////////////////////////////////////////////////////////////////////////////
    // make a union of occupied and unoccupied orbs
    RowMatrixXd C_p_eig = RowMatrixXd::Zero(nao, ndocc + n_unocc);
    TArray C_p;
    {
      RowMatrixXd C_a_eig = array_ops::array_to_eigen(C_a);
      C_p_eig.block(0, 0, nao, ndocc) << C_m_eig;
      C_p_eig.block(0, ndocc, nao, n_unocc) << C_a_eig;

      auto tr_all = utility::join_trange1(m_space.trange(), a_space.trange());
      C_p = array_ops::eigen_to_array<Tile, Policy>(
          world, C_p_eig, obs_basis->create_trange1(), tr_all);

      std::vector<double> occ_p(C_p_eig.cols(), 0);
      std::fill(occ_p.begin(), occ_p.begin() + ndocc, 2.0);
      auto p_space =
          POrbSpace(OrbitalIndex(L"p"), OrbitalIndex(L"κ"), C_p, occ_p);
      orbital_registry.add(p_space);
    }
  } else {  // received the full space
    const auto &p_space = *input_space;
    const auto &occ_p = p_space.attributes();

    // divide the LCAO space into subspaces using Eigen .. boo
    RowMatrixXd C_p = array_ops::array_to_eigen(p_space.coefs());
    const auto n = C_p.cols();
    const auto nao = C_p.rows();
    RowMatrixXd C_m = C_p.block(0, 0, C_p.rows(), ndocc);
    RowMatrixXd C_i = C_m.block(0, n_frozen_core, nao, ndocc - n_frozen_core);
    RowMatrixXd C_a = C_p.rightCols(n - ndocc);

    auto tre = std::make_shared<TRange1Engine>(ndocc, n, occ_blksize,
                                               unocc_blksize, n_frozen_core);

    // get all the trange1s
    auto tr_ao = p_space.coefs().trange().data()[0];
    auto tr_i = tre->get_active_occ_tr1();
    auto tr_m = utility::compute_trange1(ndocc, occ_blksize);
    auto tr_a = tre->get_vir_tr1();
    auto tr_p = tre->get_all_tr1();

    // convert eigen arrays to TA
    auto C_m_ta =
        array_ops::eigen_to_array<Tile, Policy>(world, C_m, tr_ao, tr_m);
    auto C_i_ta =
        array_ops::eigen_to_array<Tile, Policy>(world, C_i, tr_ao, tr_i);
    auto C_a_ta =
        array_ops::eigen_to_array<Tile, Policy>(world, C_a, tr_ao, tr_a);

    // create orbital spaces and push into registry
    std::vector<double> occ_m(occ_p.begin(), occ_p.begin() + ndocc);
    auto m_space =
        POrbSpace(OrbitalIndex(L"m"), OrbitalIndex(L"κ"), C_m_ta, occ_m);
    orbital_registry.add(m_space);

    std::vector<double> occ_i(occ_m.begin() + n_frozen_core, occ_m.end());
    auto i_space =
        POrbSpace(OrbitalIndex(L"i"), OrbitalIndex(L"κ"), C_i_ta, occ_i);
    orbital_registry.add(i_space);

    std::vector<double> occ_a(occ_p.begin() + ndocc, occ_p.end());
    auto a_space =
        POrbSpace(OrbitalIndex(L"a"), OrbitalIndex(L"κ"), C_a_ta, occ_a);
    orbital_registry.add(a_space);

    // reblock the full space and add to the regisry
    auto C_p_ta =
        array_ops::eigen_to_array<Tile, Policy>(world, C_p, tr_ao, tr_p);
    auto p_space_reblocked =
        POrbSpace(OrbitalIndex(L"p"), OrbitalIndex(L"κ"), C_p_ta, occ_p);
    orbital_registry.add(p_space_reblocked);
  }
}

// produces canonical eigenvalues and coefficients
template <typename Tile, typename Policy>
std::shared_ptr<CanonicalOrbitalSpace<TA::DistArray<Tile, Policy>>>
make_closed_shell_canonical_orbitals(
    std::shared_ptr<gaussian::AOFactoryBase<Tile, Policy>> ao_factory,
    std::size_t ndocc, std::size_t target_blocksize) {
  using TRange1Engine = ::mpqc::utility::TRange1Engine;

  auto &world = ao_factory->world();

  RowMatrixXd F_eig =
      array_ops::array_to_eigen(ao_factory->compute(L"<κ|F|λ>"));
  auto S = ao_factory->compute(L"<κ|λ>");
  RowMatrixXd S_eig = array_ops::array_to_eigen(S);

  // solve mo coefficients
  Eigen::GeneralizedSelfAdjointEigenSolver<RowMatrixXd> es(F_eig, S_eig);
  auto evals = es.eigenvalues();
  auto C = es.eigenvectors();

  // convert coeffs to TA
  auto nobs = S_eig.rows();
  using TRange1Engine = ::mpqc::utility::TRange1Engine;
  auto tre = std::make_shared<TRange1Engine>(ndocc, nobs, target_blocksize,
                                             target_blocksize, 0);
  auto tr_ao = S.trange().data().back();
  auto tr_all = tre->get_all_tr1();
  auto C_obs = array_ops::eigen_to_array<Tile, Policy>(world, C, tr_ao, tr_all);

  // convert eigenvalues to std::vec
  std::vector<double> evals_vec(evals.data(), evals.data() + evals.rows());

  return std::make_shared<CanonicalOrbitalSpace<TA::DistArray<Tile, Policy>>>(
      OrbitalIndex(L"p"), OrbitalIndex(L"κ"), C_obs, evals_vec);
}

template <typename Tile, typename Policy>
void closed_shell_cabs_mo_build_svd(
    gaussian::AOFactoryBase<Tile, Policy> &ao_factory,
    const std::shared_ptr<const ::mpqc::utility::TRange1Engine> tre,
    std::size_t vir_blocksize) {
  auto &orbital_registry = ao_factory.orbital_registry();
  auto &world = ao_factory.world();
  // CABS fock build
  auto mo_time0 = mpqc::fenced_now(world);
  utility::print_par(world, "\nBuilding ClosedShell CABS MO Orbital\n");

  // build the RI basis

  const auto abs_basis =
      *ao_factory.basis_registry()->retrieve(OrbitalIndex(L"α"));
  const auto obs_basis =
      *ao_factory.basis_registry()->retrieve(OrbitalIndex(L"κ"));

  gaussian::Basis ri_basis;
  ri_basis = merge(obs_basis, abs_basis);

  mpqc::detail::parallel_print_range_info(world, ri_basis.create_trange1(),
                                          "RI Basis");
  ao_factory.basis_registry()->add(OrbitalIndex(L"ρ"),
                                   std::make_shared<gaussian::Basis>(ri_basis));

  // integral
  auto S_ribs_inv = ao_factory.compute(L"<ρ|σ>[inv_sqr]");
  auto S_obs_ribs = ao_factory.compute(L"<μ|σ>");
  auto S_obs_inv = ao_factory.compute(L"<κ|λ>[inv_sqr]");

  // construct cabs
  TA::DistArray<Tile, Policy> C_cabs, C_ri, C_allvir;
  {
    // orthogonalize
    decltype(S_obs_inv) S_obs_ribs_ortho;
    S_obs_ribs_ortho("i,j") =
        S_obs_inv("i,k") * S_obs_ribs("k,l") * S_ribs_inv("l,j");
    RowMatrixXd S_obs_ribs_ortho_eigen =
        array_ops::array_to_eigen(S_obs_ribs_ortho);

    // SVD solve
    Eigen::JacobiSVD<RowMatrixXd> svd(S_obs_ribs_ortho_eigen,
                                      Eigen::ComputeFullV);
    RowMatrixXd V_eigen = svd.matrixV();
    size_t nbf_ribs = S_obs_ribs_ortho_eigen.cols();
    auto nbf_cabs = nbf_ribs - svd.nonzeroSingularValues();
    RowMatrixXd Vnull(nbf_ribs, nbf_cabs);
    Vnull = V_eigen.block(0, svd.nonzeroSingularValues(), nbf_ribs, nbf_cabs);

    auto tr_ribs = ri_basis.create_trange1();
    auto tr_cabs_mo = utility::compute_trange1(nbf_cabs, vir_blocksize);
    mpqc::detail::parallel_print_range_info(world, tr_cabs_mo, "CABS MO");

    C_cabs = array_ops::eigen_to_array<Tile, Policy>(world, Vnull, tr_ribs,
                                                     tr_cabs_mo);
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
    auto tr_ribs_mo = utility::compute_trange1(nbf_ribs, vir_blocksize);
    mpqc::detail::parallel_print_range_info(world, tr_ribs_mo, "RIBS MO");
    auto ribs_to_mo = array_ops::create_diagonal_array_from_eigen<Tile, Policy>(
        world, tr_ribs, tr_ribs_mo, 1.0);
    C_ri("i,j") = S_ribs_inv("i,k") * ribs_to_mo("k,j");

    auto tr_allvir_mo = utility::compute_trange1(nbf_cabs + n_vir, vir_blocksize);
    mpqc::detail::parallel_print_range_info(world, tr_allvir_mo,
                                            "All Virtual MO");

    C_allvir = array_ops::eigen_to_array<Tile, Policy>(
        world, C_allvirtual_eigen, tr_ribs, tr_allvir_mo);

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
std::shared_ptr<::mpqc::utility::TRange1Engine>
closed_shell_dualbasis_mo_build_eigen_solve_svd(
    LCAOFactoryBase<Tile, Policy> &lcao_factory, Eigen::VectorXd &ens,
    std::size_t nocc, const Molecule &mols, bool frozen_core,
    std::size_t occ_blocksize, std::size_t vir_blocksize) {
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
  auto tr_occ = utility::compute_trange1(nocc, occ_blocksize);
  auto tr_corr_occ = tre->get_active_occ_tr1();
  auto tr_vir = tre->get_vir_tr1();

  mpqc::detail::parallel_print_range_info(world, tr_occ, "Occ");
  mpqc::detail::parallel_print_range_info(world, tr_corr_occ, "CorrOcc");
  mpqc::detail::parallel_print_range_info(world, tr_vir, "Vir");

  // convert to TA
  auto C_occ_ta =
      array_ops::eigen_to_array<Tile, Policy>(world, C_occ, tr_obs, tr_occ);
  auto C_corr_occ_ta = array_ops::eigen_to_array<Tile, Policy>(
      world, C_corr_occ, tr_obs, tr_corr_occ);
  auto C_vir_ta =
      array_ops::eigen_to_array<Tile, Policy>(world, C_vbs, tr_vbs, tr_vir);

  // insert to registry
  using OrbitalSpaceTArray = OrbitalSpace<TA::DistArray<Tile, Policy>>;
  auto occ_space =
      OrbitalSpaceTArray(OrbitalIndex(L"m"), OrbitalIndex(L"κ"), C_occ_ta);
  lcao_factory.orbital_registry().add(occ_space);

  auto corr_occ_space =
      OrbitalSpaceTArray(OrbitalIndex(L"i"), OrbitalIndex(L"κ"), C_corr_occ_ta);
  lcao_factory.orbital_registry().add(corr_occ_space);

  auto vir_space =
      OrbitalSpaceTArray(OrbitalIndex(L"a"), OrbitalIndex(L"Α"), C_vir_ta);
  lcao_factory.orbital_registry().add(vir_space);

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
      array_ops::eigen_to_array<Tile, Policy>(world, C_vbs, tr_vbs, tr_vir);

  // remove old virtual orbitals
  lcao_factory.orbital_registry().remove(OrbitalIndex(L"a"));
  lcao_factory.registry().purge_index(world, OrbitalIndex(L"a"));

  // add new virtual orbial
  vir_space =
      OrbitalSpaceTArray(OrbitalIndex(L"a"), OrbitalIndex(L"Α"), C_vir_ta_new);
  lcao_factory.orbital_registry().add(vir_space);

  auto mo_time1 = mpqc::fenced_now(world);
  auto mo_time = mpqc::duration_in_s(mo_time0, mo_time1);
  utility::print_par(world, "ClosedShell Dual Basis MO Build Time: ", mo_time,
                     " S\n");

  return tre;
}

template <typename Tile, typename Policy>
void closed_shell_dualbasis_cabs_mo_build_svd(
    LCAOFactoryBase<Tile, Policy> &lcao_factory,
    const std::shared_ptr<::mpqc::utility::TRange1Engine> tre,
    std::string ri_method, std::size_t vir_blocksize) {
  auto &ao_factory = lcao_factory.ao_factory();
  auto &orbital_registry = lcao_factory.orbital_registry();
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
      utility::compute_trange1(tr_cabs.elements_range().second, vir_blocksize);
  auto tr_allvir_mo = utility::compute_trange1(nbf_ribs_minus_occ, vir_blocksize);

  mpqc::detail::parallel_print_range_info(world, tr_cabs_mo, "CABS MO");
  TA::DistArray<Tile, Policy> C_cabs = array_ops::eigen_to_array<Tile, Policy>(
      world, C_cabs_eigen, tr_ribs, tr_cabs_mo);

  // insert to orbital space
  using OrbitalSpaceTArray = OrbitalSpace<TA::DistArray<Tile, Policy>>;
  auto C_cabs_space =
      OrbitalSpaceTArray(OrbitalIndex(L"a'"), OrbitalIndex(L"ρ"), C_cabs);
  orbital_registry.add(C_cabs_space);

  mpqc::detail::parallel_print_range_info(world, tr_allvir_mo,
                                          "All Virtual MO");
  TA::DistArray<Tile, Policy> C_allvir =
      array_ops::eigen_to_array<Tile, Policy>(world, C_allvir_eigen, tr_ribs,
                                              tr_allvir_mo);

  // insert to orbital space
  auto C_allvir_space =
      OrbitalSpaceTArray(OrbitalIndex(L"A'"), OrbitalIndex(L"ρ"), C_allvir);
  orbital_registry.add(C_allvir_space);

  auto mo_time1 = mpqc::fenced_now(world);
  auto mo_time = mpqc::duration_in_s(mo_time0, mo_time1);
  utility::print_par(world, "ClosedShell Dual Basis CABS MO Build Time: ",
                     mo_time, " S\n");
}

/*!
 * \brief This inserts crystal orbitals to registry for gamma-point methods
 */
template <typename Tile, typename Policy>
std::shared_ptr<::mpqc::utility::TRange1Engine> mo_insert_gamma_point(PeriodicLCAOFactory<Tile, Policy>& plcao_factory,
                           RowMatrixXd& C_gamma_point, UnitCell& unitcell,
                           size_t occ_block, size_t vir_block) {
  using TRange1Engine = ::mpqc::utility::TRange1Engine;

  auto& orbital_registry = plcao_factory.orbital_registry();
  auto& world = plcao_factory.world();

  auto all = C_gamma_point.cols();

  // the unit cell must be electrically neutral
  const auto charge = 0;
  auto occ = (unitcell.total_atomic_number() - charge) / 2;
  auto vir = all - occ;
  std::size_t n_frozen_core = 0;  // TODO: should be determined by user

  RowMatrixXd C_occ = C_gamma_point.leftCols(occ);
  RowMatrixXd C_corr_occ =
      C_gamma_point.block(0, n_frozen_core, all, occ - n_frozen_core);
  RowMatrixXd C_vir = C_gamma_point.rightCols(vir);

  ExEnv::out0() << "OccBlockSize: " << occ_block << std::endl;
  ExEnv::out0() << "VirBlockSize: " << vir_block << std::endl;

  auto obs_basis =
      plcao_factory.pao_factory().basis_registry()->retrieve(
          OrbitalIndex(L"κ"));
  auto tre = std::make_shared<TRange1Engine>(occ, all, occ_block, vir_block, 0);

  // get all trange1s
  auto tr_obs = obs_basis->create_trange1();
  auto tr_corr_occ = tre->get_active_occ_tr1();
  auto tr_occ = utility::compute_trange1(occ, occ_block);
  auto tr_vir = tre->get_vir_tr1();
  auto tr_all = tre->get_all_tr1();

  mpqc::detail::parallel_print_range_info(world, tr_obs, "Obs");
  mpqc::detail::parallel_print_range_info(world, tr_occ, "Occ");
  mpqc::detail::parallel_print_range_info(world, tr_corr_occ, "CorrOcc");
  mpqc::detail::parallel_print_range_info(world, tr_vir, "Vir");
  mpqc::detail::parallel_print_range_info(world, tr_all, "All");

  // convert Eigen matrices to TA
  auto C_occ_ta =
      array_ops::eigen_to_array<Tile, Policy>(world, C_occ, tr_obs, tr_occ);
  auto C_corr_occ_ta = array_ops::eigen_to_array<Tile, Policy>(
      world, C_corr_occ, tr_obs, tr_corr_occ);
  auto C_vir_ta =
      array_ops::eigen_to_array<Tile, Policy>(world, C_vir, tr_obs, tr_vir);
  auto C_all_ta = array_ops::eigen_to_array<Tile, Policy>(world, C_gamma_point,
                                                          tr_obs, tr_all);

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

  auto all_space =
      OrbitalSpaceTArray(OrbitalIndex(L"p"), OrbitalIndex(L"κ"), C_all_ta);
  orbital_registry.add(all_space);

  return tre;
}



}  // namespace lcao
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_SCF_MO_BUILD_H_
