//
// Created by Chong Peng on 12/6/16.
//

#ifndef SRC_MPQC_CHEMISTRY_QC_MBPT_DBMP2_IMPL_H_
#define SRC_MPQC_CHEMISTRY_QC_MBPT_DBMP2_IMPL_H_


namespace mpqc {
namespace lcao {

namespace detail {
template <typename Tile, typename Policy>
std::shared_ptr<TRange1Engine> closed_shell_dual_basis_mo_build_steele(
    LCAOFactory<Tile, Policy> &lcao_factory,
    Eigen::VectorXd &ens,
    const Molecule &mols,
    bool frozen_core,
    std::size_t occ_blocksize,
    std::size_t vir_blocksize)
{
  auto &ao_factory = lcao_factory.ao_factory();
  auto& world = ao_factory.world();
  auto& orbital_registry = lcao_factory.orbital_space();
  using TArray = TA::DistArray<Tile, Policy>;

  auto mo_time0 = mpqc::fenced_now(world);
  auto occ = mols.occupation() / 2;
  utility::print_par(
      world,
      "\nBuilding ClosedShell Dual Basis MO Orbital (Steele Approach) \n");

  // solving occupied orbitals
  TArray F;
  if (ao_factory.registry().have(Formula(L"<μ|F|ν>"))) {
    F = ao_factory.registry().retrieve(Formula(L"<μ|F|ν>"));
  } else {
    F = ao_factory.registry().retrieve(Formula(L"<μ|F|ν>[df]"));
  }

  auto S = ao_factory.compute(L"<κ|λ>");

  RowMatrixXd F_eig = array_ops::array_to_eigen(F);
  RowMatrixXd S_eig = array_ops::array_to_eigen(S);

  // solve mo coefficients
  Eigen::GeneralizedSelfAdjointEigenSolver<RowMatrixXd> es(F_eig, S_eig);

  std::size_t n_frozen_core = 0;
  if (frozen_core) {
    n_frozen_core = mols.core_electrons();
    utility::print_par(world, "Frozen Core: ", n_frozen_core, " electrons. \n");
    n_frozen_core = n_frozen_core / 2;
  }

  std::cout << es.eigenvalues() << std::endl;
  RowMatrixXd C_occ = es.eigenvectors().leftCols(occ);

  // finished solving occupied orbitals

  // get all the sizes
  utility::print_par(world, "OccBlockSize: ", occ_blocksize, "\n");
  utility::print_par(world, "VirBlockSize: ", vir_blocksize, "\n");

  std::size_t n_obs = S.trange().elements_range().extent()[0];
  auto tre = std::make_shared<TRange1Engine>(occ, n_obs, occ_blocksize,
                                             vir_blocksize, n_frozen_core);
  auto tr_obs = S.trange().data().back();
  auto tr_occ = tre->compute_range(occ, occ_blocksize);
  auto tr_corr_occ = tre->get_active_occ_tr1();
  mpqc::detail::parallel_print_range_info(world, tr_occ, "Occ");

  // convert to TA
  auto C_occ_ta = array_ops::eigen_to_array<Tile,Policy>(world, C_occ, tr_obs, tr_occ);

  // project to large basis set

  TArray S_vbs_inv = ao_factory.compute(L"<Α|Β>[inv]");
  TArray S_vbs_obs = ao_factory.compute(L"<Α|μ>");
  TArray S_vbs_vbs = ao_factory.compute(L"<Α|Β>");
  auto n_vbs = S_vbs_inv.trange().elements_range().extent()[0];
  auto tr_vbs = S_vbs_inv.trange().data().back();

  TArray t;
  t("A,mu") = S_vbs_inv("A,B") * S_vbs_obs("B,mu");

  C_occ_ta("A,i") = t("A,mu") * C_occ_ta("mu,i");

  // insert to registry
  using OrbitalSpaceTArray = OrbitalSpace<TA::DistArray<Tile, Policy>>;
  auto occ_space =
      OrbitalSpaceTArray(OrbitalIndex(L"m"), OrbitalIndex(L"Α"), C_occ_ta);
  orbital_registry.add(occ_space);

  // solving the new orbitals in large basis set
  // find fock matrix
  TArray F_vbs;
  // if use density fitting
  if (ao_factory.registry().have(Formula(L"<μ|F|ν>[df]"))) {
    F_vbs = ao_factory.compute(Formula(L"<Α|F|Β>[df]"));
  } else {
    F_vbs = ao_factory.compute(Formula(L"<Α|F|Β>"));
  }
  tre = std::make_shared<TRange1Engine>(occ, n_vbs, occ_blocksize,
                                        vir_blocksize, n_frozen_core);

  RowMatrixXd C_vir;
  RowMatrixXd C_corr_occ;
  {
    RowMatrixXd F_vbs_eigen = array_ops::array_to_eigen(F_vbs);
    RowMatrixXd S_vbs_eigen = array_ops::array_to_eigen(S_vbs_vbs);
    Eigen::GeneralizedSelfAdjointEigenSolver<RowMatrixXd> es(F_vbs_eigen, S_vbs_eigen);

    assert(es.info() == Eigen::ComputationInfo::Success);
    std::cout << es.eigenvalues() << std::endl;
    ens = es.eigenvalues();
    auto env = es.eigenvectors();

    C_occ = env.leftCols(occ);
    C_corr_occ = C_occ.rightCols(occ - n_frozen_core);
    C_vir = env.rightCols(n_vbs - occ);
  }

  // get all the trange1s
  auto tr_vir = tre->get_vir_tr1();
  mpqc::detail::parallel_print_range_info(world, tr_vir, "Vir");

  C_occ_ta = array_ops::eigen_to_array<Tile,Policy>(world, C_occ, tr_vbs, tr_occ);

  auto C_corr_occ_ta =
      array_ops::eigen_to_array<Tile,Policy>(world, C_corr_occ, tr_vbs, tr_corr_occ);

  auto C_vir_ta = array_ops::eigen_to_array<Tile,Policy>(world, C_vir, tr_vbs, tr_vir);

  // insert to registry
  occ_space =
      OrbitalSpaceTArray(OrbitalIndex(L"m"), OrbitalIndex(L"Α"), C_occ_ta);
  orbital_registry.remove(OrbitalIndex(L"m"));
  orbital_registry.add(occ_space);

  auto corr_occ_space =
      OrbitalSpaceTArray(OrbitalIndex(L"i"), OrbitalIndex(L"Α"), C_corr_occ_ta);
  orbital_registry.add(corr_occ_space);

  auto vir_space =
      OrbitalSpaceTArray(OrbitalIndex(L"a"), OrbitalIndex(L"Α"), C_vir_ta);
  orbital_registry.add(vir_space);

  auto mo_time1 = mpqc::fenced_now(world);
  auto mo_time = mpqc::duration_in_s(mo_time0, mo_time1);
  utility::print_par(world, "ClosedShell Dual Basis MO Build Time: ", mo_time,
                     " S\n");

  return tre;
}
} // end of namespace detail


///
/// member function of DBRMP2
///

template<typename Tile, typename Policy>
void DBRMP2<Tile,Policy>::obsolete() {
  scf_correction_ = 0.0;
  RMP2<Tile,Policy>::obsolete();
}


template<typename Tile, typename Policy>
void DBRMP2<Tile,Policy>::evaluate(Energy *result) {

  if(!this->compute()){

    // call RMP2 evaluate function
    RMP2<Tile,Policy>::evaluate(result);

    double mp2_energy = this->get_value(result).derivs(0)[0];

    // compute scf correction
    double scf_correction = compute_scf_correction();

    this->computed_ = true;
    this->set_value(result, mp2_energy + scf_correction);
  }
}

template<typename Tile, typename Policy>
void DBRMP2<Tile,Policy>::init() {
  // if not initialized
  if (this->trange1_engine() == nullptr || this->orbital_energy() == nullptr) {
    auto mol = this->wfn_world()->atoms();
    Eigen::VectorXd orbital_energy;

    if (method_ == "valeev") {
      this->trange1_engine_ = closed_shell_dualbasis_mo_build_eigen_solve_svd(
          this->lcao_factory(), orbital_energy, *mol, this->is_frozen_core(),
          this->occ_block(), this->unocc_block());
    } else if (method_ == "steele") {
      this->trange1_engine_ = detail::closed_shell_dual_basis_mo_build_steele(
          this->lcao_factory(), orbital_energy, *mol, this->is_frozen_core(),
          this->occ_block(), this->unocc_block());
    }
    this->orbital_energy_ = std::make_shared<Eigen::VectorXd>(orbital_energy);
  }
}

template<typename Tile, typename Policy>
double DBRMP2<Tile,Policy>::compute_scf_correction() {
  init();

  int occ = this->trange1_engine()->get_occ();
  TArray F_ma;

  F_ma = this->lcao_factory().compute(L"<m|F|a>");
  //    } else if (method == "df") {
  //      F_ma = this->lcao_factory().compute(L"<m|F|a>[df]");

  double scf_correction = 2 *
      F_ma("m,a").reduce(detail::ScfCorrection<Tile>(
          this->orbital_energy(), occ));

  if (F_ma.world().rank() == 0) {
    std::cout << "SCF Correction: " << scf_correction << std::endl;
  }

  return scf_correction;
}

///
/// member function of RIDBRMP2
///
template<typename Tile, typename Policy>
/// \return
double RIDBRMP2<Tile,Policy>::compute() {
  return detail::compute_mp2(this->lcao_factory(), this->orbital_energy(),
                             this->trange1_engine(), true);
}

template<typename Tile, typename Policy>
double RIDBRMP2<Tile,Policy>::compute_scf_correction() {
  this->init();

  int occ = this->trange1_engine()->get_occ();
  TA::DistArray<Tile,Policy> F_ma;

  F_ma = this->lcao_factory().compute(L"<m|F|a>[df]");

  double scf_correction = 2 *
      F_ma("m,a").reduce(detail::ScfCorrection<Tile>(
          this->orbital_energy(), occ));

  if (F_ma.world().rank() == 0) {
    std::cout << "SCF Correction: " << scf_correction << std::endl;
  }

  return scf_correction;
}

}  // namespace lcao
}  // namespace mpqc


#endif //SRC_MPQC_CHEMISTRY_QC_MBPT_DBMP2_IMPL_H_
