//
// Created by Chong Peng on 7/12/16.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_F12_DBCCSD_F12_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_F12_DBCCSD_F12_H_

#include "mpqc/chemistry/qc/f12/ccsd_f12.h"
#include "mpqc/chemistry/qc/f12/db_f12_intermediates.h"
#include "mpqc/mpqc_config.h"

namespace mpqc {
namespace lcao {

/***
 * \class DBCCSDF12 Dual basis CCSD F12 class
 *
 * takes all options from CCSDF12 and the following options:
 *
 * @param RIMethod: OBS or VBS, if OBS, RI Basis is union of OBS+CABS; if VBS,
 * RI Basis is union of VBS+CABS, default VBS
 *
 */

template <typename Tile>
class DBCCSD_F12 : public CCSD_F12<Tile> {
 public:
  using Policy = TA::SparsePolicy;
  using TArray = TA::DistArray<Tile, Policy>;
  using LCAOFactoryType = LCAOFactory<Tile, Policy>;

  using real_t = typename Tile::scalar_type;
  using Matrix = RowMatrix<real_t>;

  /**
   * KeyVal constructor
   * @param kv
   *
   * keywords: takes all keywords from CCSDF12 class
   *
   * invalid keywords: approximation, vt_couple
   */
  DBCCSD_F12(const KeyVal& kv)
      : CCSD<Tile, TA::SparsePolicy>(kv), CCSD_F12<Tile>(kv) {}
  virtual ~DBCCSD_F12() = default;

 private:
  /// overide initialization of CCSD
  void init() override {
    if (this->orbital_energy() == nullptr ||
        this->trange1_engine() == nullptr) {
      auto& lcao_factory = this->lcao_factory();
      auto mol = lcao_factory.ao_factory().molecule();
      Eigen::VectorXd orbital_energy;
      this->trange1_engine_ = closed_shell_dualbasis_mo_build_eigen_solve_svd(
          lcao_factory, orbital_energy, mol, this->is_frozen_core(),
          this->occ_block(), this->unocc_block());
      this->orbital_energy_ = std::make_shared<Eigen::VectorXd>(orbital_energy);
    }
  }

  /// override initialization of CABS in CCSDF12
  void init_cabs() override {
    closed_shell_dualbasis_cabs_mo_build_svd(this->lcao_factory(),
                                             this->trange1_engine(), "VBS",
                                             this->unocc_block());
  }

  /// override cabs singles in CCSDF12
  void compute_cabs_singles() override;

  /// override compute_f12 in CCSDF12
  void compute_f12() override;

  Matrix compute_db_ccsd_f12_df();
};

template <typename Tile>
void DBCCSD_F12<Tile>::compute_cabs_singles() {
  auto& world = this->wfn_world()->world();

  mpqc::utility::print_par(world, " CABS Singles \n");
  auto single_time0 = mpqc::fenced_now(world);

  // non-canonical, don't include F_m^a
  CABSSingles<Tile> cabs_singles(this->lcao_factory());
  this->singles_energy_ = cabs_singles.compute(true, false, false);
  utility::print_par(world, "E_S: ", this->singles_energy_, "\n");
  auto single_time1 = mpqc::fenced_now(world);
  auto single_time = mpqc::duration_in_s(single_time0, single_time1);
  mpqc::utility::print_par(world, "Total CABS Singles Time:  ", single_time,
                           "\n");
}

template <typename Tile>
void DBCCSD_F12<Tile>::compute_f12() {
  auto& world = this->wfn_world()->world();

  Matrix Eij_F12;
  if (this->method_ != "df") {
    utility::print_par(world, "\n Warning! DBCCSD only support method==df! \n");
  }
  Eij_F12 = compute_db_ccsd_f12_df();

  this->f12_energy_ = Eij_F12.sum();
}

template <typename Tile>
typename DBCCSD_F12<Tile>::Matrix DBCCSD_F12<Tile>::compute_db_ccsd_f12_df() {
  auto& lcao_factory = this->lcao_factory();
  auto& world = lcao_factory.world();
  Matrix Eij_F12;

  auto nocc = this->trange1_engine()->get_active_occ();

  // create shape
  auto occ_tr1 = this->trange1_engine()->get_active_occ_tr1();
  TiledArray::TiledRange occ4_trange({occ_tr1, occ_tr1, occ_tr1, occ_tr1});
  auto ijij_ijji_shape = f12::make_ijij_ijji_shape(occ4_trange);

  // compute V_ijij_ijji
  TArray V_ijij_ijji = f12::compute_V_ijij_ijji_db_df(lcao_factory, ijij_ijji_shape);

  // VT2 contribution
  TArray tmp =
      mpqc::lcao::f12::compute_VT2_ijij_ijji_db_df(lcao_factory, this->t2(), ijij_ijji_shape);
  V_ijij_ijji("i1,j1,i2,j2") += tmp("i1,j1,i2,j2");

  // VT1 contribution
  {
    TArray tmp =
        f12::compute_VT1_ijij_ijji_db_df(lcao_factory, this->t1(), ijij_ijji_shape);
    V_ijij_ijji("i1,j1,i2,j2") += tmp("i1,j1,i2,j2");
  }

  // V contribution to energy
  Eij_F12 = V_ijij_ijji("i1,j1,i2,j2")
                .reduce(f12::F12PairEnergyReductor<Tile>(
                    2 * f12::C_ijij_bar, 2 * f12::C_ijji_bar, nocc));
  utility::print_par(world, "E_V: ", Eij_F12.sum(), "\n");

  // compute X term
  TArray X_ijij_ijji = f12::compute_X_ijij_ijji_db_df(lcao_factory, ijij_ijji_shape);

  auto Fij = this->lcao_factory().compute(L"(i|F|j)[df]");
  auto Fij_eigen = array_ops::array_to_eigen(Fij);
  f12::convert_X_ijkl(X_ijij_ijji, Fij_eigen);
  {
    Matrix eij = X_ijij_ijji("i1,j1,i2,j2")
                     .reduce(f12::F12PairEnergyReductor<Tile>(
                         f12::CC_ijij_bar, f12::CC_ijji_bar, nocc));
    eij *= -1;
    utility::print_par(world, "E_X: ", eij.sum(), "\n");
    Eij_F12 += eij;
  }

  // compute B term
  {
    TArray B_ijij_ijji =
        f12::compute_B_ijij_ijji_db_df(lcao_factory, ijij_ijji_shape);
    Matrix eij = B_ijij_ijji("i1,j1,i2,j2")
                     .reduce(f12::F12PairEnergyReductor<Tile>(
                         f12::CC_ijij_bar, f12::CC_ijji_bar, nocc));
    utility::print_par(world, "E_B: ", eij.sum(), "\n");
    Eij_F12 += eij;
  }

  return Eij_F12;
}

#if TA_DEFAULT_POLICY == 1
extern template class DBCCSD_F12<TA::TensorD>;
#endif

}  // namespace lcao
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_F12_DBCCSD_F12_H_
