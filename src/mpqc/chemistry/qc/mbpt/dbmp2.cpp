//
// Created by Chong Peng on 11/2/16.
//

#include "mpqc/chemistry/qc/mbpt/dbmp2.h"
#include "mpqc/util/keyval/forcelink.h"

MPQC_CLASS_EXPORT2("DBRMP2", mpqc::lcao::DBRMP2);
MPQC_CLASS_EXPORT2("RI-DBRMP2", mpqc::lcao::RIDBRMP2);

namespace mpqc {
namespace lcao {
///
/// member function of DBRMP2
///
DBRMP2::~DBRMP2() = default;

double DBRMP2::value() {
  if (this->energy_ == 0.0) {
    double mp2_energy = RMP2::value();

    double scf_correction = compute_scf_correction();
    this->energy_ = scf_correction + mp2_energy;
  }
  return this->energy_;
}

void DBRMP2::init() {
  // if not initialized
  if (this->trange1_engine() == nullptr || this->orbital_energy() == nullptr) {
    auto mol = this->wfn_world()->molecule();
    Eigen::VectorXd orbital_energy;

    if (method_ == "valeev") {
      this->trange1_engine_ = closed_shell_dualbasis_mo_build_eigen_solve_svd(
          this->lcao_factory(), orbital_energy, mol, this->is_frozen_core(),
          this->occ_block(), this->unocc_block());
    } else if (method_ == "steele") {
      this->trange1_engine_ = detail::closed_shell_dual_basis_mo_build_steele(
          this->lcao_factory(), orbital_energy, mol, this->is_frozen_core(),
          this->occ_block(), this->unocc_block());
    }
    this->orbital_energy_ = std::make_shared<Eigen::VectorXd>(orbital_energy);
  }
}

double DBRMP2::compute_scf_correction() {
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

double RIDBRMP2::compute() {
  return mbpt::detail::compute_mp2(this->lcao_factory(), this->orbital_energy(),
                                   this->trange1_engine(), true);
}

double RIDBRMP2::compute_scf_correction() {
  init();

  int occ = this->trange1_engine()->get_occ();
  TArray F_ma;

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
