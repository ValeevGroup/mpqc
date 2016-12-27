//
// Created by Chong Peng on 4/12/16.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_F12_CCSD_F12_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_F12_CCSD_F12_H_

#include <tiledarray.h>

#include "mpqc/chemistry/qc/cc/ccsd.h"
#include "mpqc/chemistry/qc/f12/f12_intermediates.h"
#include "mpqc/chemistry/qc/wfn/trange1_engine.h"
#include "mpqc/chemistry/qc/f12/cabs_singles.h"
#include "mpqc/mpqc_config.h"

namespace mpqc {
namespace lcao {

/**
 *  CCSD(2)F12 class
 *
 */

template <typename Tile>
class CCSD_F12 : virtual public CCSD<Tile, TA::SparsePolicy> {
 public:
  using Policy = TA::SparsePolicy;
  using TArray = TA::DistArray<Tile, Policy>;
  using DirectArray =
      typename CCSD<Tile, Policy>::DirectAOIntegral::DirectTArray;
  using LCAOFactoryType = LCAOFactory<Tile, Policy>;

  using real_t = typename Tile::scalar_type;
  using Matrix = RowMatrix<real_t>;

  /**
   * KeyVal constructor
   *
   * @param kv
   *
   * keywords: takes all keywords from CCSD class
   *
   * | KeyWord | Type | Default| Description |
   * |---------|------|--------|-------------|
   * | approx | char | C | approaximation to compute F12 (C or D) |
   * | cabs_singles | bool | true | if do CABSSingles calculation |
   * | vt_couple | bool | true | if couple last two term in VT2 and VT1 term |
   *
   */
  CCSD_F12(const KeyVal& kv) : CCSD<Tile, Policy>(kv) {
    vt_couple_ = kv.value<bool>("vt_couple", true);
    cabs_singles_ = kv.value<bool>("cabs_singles", true);

    approximation_ = kv.value<char>("approx", 'C');
    if (approximation_ != 'C' && approximation_ != 'D') {
      throw std::runtime_error("Wrong CCSD_F12 Approach");
    }

    method_ = kv.value<std::string>("method", "df");
    if (method_ != "standard" && method_ != "df" && method_ != "direct") {
      throw std::invalid_argument("Invalid Method For CCSD_F12");
    }

    f12_energy_ = 0.0;
    singles_energy_ = 0.0;
  }

  virtual ~CCSD_F12() = default;

  double value() override {
    if (this->energy_ == 0.0) {
      auto& world = this->wfn_world()->world();

      auto ccsd_time0 = mpqc::fenced_now(world);

      // compute ccsd
      real_t ccsd = CCSD<Tile, Policy>::value();

      auto ccsd_time1 = mpqc::fenced_now(world);
      auto ccsd_time = mpqc::duration_in_s(ccsd_time0, ccsd_time1);
      mpqc::utility::print_par(world, "Total CCSD Time:  ", ccsd_time, "\n");

      // initialize CABS orbitals
      init_cabs();

      utility::print_par(world, "VTCouple: ", vt_couple_, "\n");

      // clean LCAO Integrals
      this->lcao_factory().registry().purge(world);

      // compute, this will set f12_energy_
      compute_f12();

      if (debug()) utility::print_par(world, "E_F12: ", f12_energy_, "\n");

      // compute cabs singles, this will set singles_energy_
      if (cabs_singles_) {
        compute_cabs_singles();
      }

      auto f12_time0 = mpqc::fenced_now(world);
      auto f12_time = mpqc::duration_in_s(ccsd_time1, f12_time0);
      mpqc::utility::print_par(world, "Total F12 Time:  ", f12_time, "\n");

      this->energy_ = ccsd + f12_energy_ + singles_energy_;
    }
    return this->energy_;
  }

  void obsolete() override {
    f12_energy_ = 0.0;
    singles_energy_ = 0.0;
    CCSD<Tile, Policy>::obsolete();
  }

 private:
  virtual void init_cabs() {
    closed_shell_cabs_mo_build_svd(this->lcao_factory(), this->trange1_engine(),
                                   this->unocc_block());
  }

  virtual void compute_cabs_singles();

  virtual void compute_f12();

  /// standard approach
  Matrix compute_ccsd_f12_df(const DirectArray& darray, const char approach);

  Matrix compute_ccsd_f12(const DirectArray& darray);

 protected:
  bool cabs_singles_;
  double f12_energy_;
  double singles_energy_;
  std::string method_;
  int debug() const { return 1; }

 private:
  bool vt_couple_;
  char approximation_;
};

template <typename Tile>
void CCSD_F12<Tile>::compute_cabs_singles() {
  auto& world = this->wfn_world()->world();

  mpqc::utility::print_par(world, " CABS Singles \n");

  auto single_time0 = mpqc::fenced_now(world);

  bool df = this->is_df();

  if (approximation_ == 'D') {
    CABSSingles<Tile> cabs_singles(this->lcao_factory());
    singles_energy_ = cabs_singles.compute(df, true, true);
  } else {
    CABSSingles<Tile> cabs_singles(this->lcao_factory());
    singles_energy_ = cabs_singles.compute(df, false, true);
  }
  if (debug()) {
    utility::print_par(world, "E_S: ", singles_energy_, "\n");
  }
  auto single_time1 = mpqc::fenced_now(world);
  auto single_time = mpqc::duration_in_s(single_time0, single_time1);
  mpqc::utility::print_par(world, "Total CABS Singles Time:  ", single_time,
                           "\n");
}

template <typename Tile>
void CCSD_F12<Tile>::compute_f12() {
  auto lazy_two_electron_int = this->get_direct_ao_integral();
  Matrix Eij_F12;
  if (method_ == "standard") {
    Eij_F12 = compute_ccsd_f12(lazy_two_electron_int);
  } else if (method_ == "df") {
    Eij_F12 = compute_ccsd_f12_df(lazy_two_electron_int, approximation_);
  } else if (method_ == "direct") {
    Eij_F12 = compute_ccsd_f12_df(lazy_two_electron_int, approximation_);
  }

  f12_energy_ = Eij_F12.sum();
}

template <typename Tile>
typename CCSD_F12<Tile>::Matrix CCSD_F12<Tile>::compute_ccsd_f12_df(
    const DirectArray& darray, const char approach) {
  auto& lcao_factory = this->lcao_factory();
  auto& world = lcao_factory.world();
  Matrix Eij_F12;

  utility::print_par(world, "\n Computing CCSD_F12 ", approach, " Approach \n");

  auto n_active_occ = this->trange1_engine()->get_active_occ();

  // create shape
  auto occ_tr1 = this->trange1_engine()->get_active_occ_tr1();
  TiledArray::TiledRange occ4_trange({occ_tr1, occ_tr1, occ_tr1, occ_tr1});
  auto ijij_ijji_shape = f12::make_ijij_ijji_shape(occ4_trange);

  // compute B term
  {
    TArray B_ijij_ijji;

    if (approach == 'C') {
      B_ijij_ijji = f12::compute_B_ijij_ijji_C_df(lcao_factory, ijij_ijji_shape);
    } else if (approach == 'D') {
      B_ijij_ijji = f12::compute_B_ijij_ijji_D_df(lcao_factory, ijij_ijji_shape);
    }

    Matrix eij = B_ijij_ijji("i1,j1,i2,j2")
                     .reduce(f12::F12PairEnergyReductor<Tile>(
                         f12::CC_ijij_bar, f12::CC_ijji_bar, n_active_occ));
    if (debug()) utility::print_par(world, "E_B: ", eij.sum(), "\n");
    Eij_F12 = eij;
  }

  lcao_factory.ao_factory().registry().purge_operator(world, L"R");

  // compute X term
  TArray X_ijij_ijji = f12::compute_X_ijij_ijji_df(lcao_factory, ijij_ijji_shape);

  lcao_factory.purge_operator(world, L"R2");

  auto Fij = lcao_factory.compute(L"<i|F|j>[df]");
  auto Fij_eigen = array_ops::array_to_eigen(Fij);
  f12::convert_X_ijkl(X_ijij_ijji, Fij_eigen);
  {
    Matrix eij = X_ijij_ijji("i1,j1,i2,j2")
                     .reduce(f12::F12PairEnergyReductor<Tile>(
                         f12::CC_ijij_bar, f12::CC_ijji_bar, n_active_occ));
    eij *= -1;
    if (debug()) utility::print_par(world, "E_X: ", eij.sum(), "\n");
    Eij_F12 += eij;
  }

  // compute V_ijij_ijji
  TArray V_ijij_ijji = f12::compute_V_ijij_ijji_df(lcao_factory, ijij_ijji_shape);

  // VT2 contribution
  if (darray.array().is_initialized()) {
    TArray tmp = f12::compute_VT2_ijij_ijji_df_direct(lcao_factory, this->t2(),
                                                 ijij_ijji_shape, darray);
    V_ijij_ijji("i1,j1,i2,j2") += tmp("i1,j1,i2,j2");
  } else {
    TArray tmp = f12::compute_VT2_ijij_ijji_df(lcao_factory, this->t2(),
                                          ijij_ijji_shape, vt_couple_);
    V_ijij_ijji("i1,j1,i2,j2") += tmp("i1,j1,i2,j2");
  }

  // VT1 contribution
  {
    TArray tmp = f12::compute_VT1_ijij_ijji_df(lcao_factory, this->t1(),
                                          ijij_ijji_shape, vt_couple_);
    V_ijij_ijji("i1,j1,i2,j2") += tmp("i1,j1,i2,j2");
  }

  // V contribution to energy
  Matrix e_ij = V_ijij_ijji("i1,j1,i2,j2")
                    .reduce(f12::F12PairEnergyReductor<Tile>(
                        2 * f12::C_ijij_bar, 2 * f12::C_ijji_bar, n_active_occ));
  Eij_F12 += e_ij;
  if (debug()) utility::print_par(world, "E_V: ", e_ij.sum(), "\n");

  return Eij_F12;
}

template <typename Tile>
typename CCSD_F12<Tile>::Matrix CCSD_F12<Tile>::compute_ccsd_f12(
    const DirectArray& darray) {
  auto& lcao_factory = this->lcao_factory();
  auto& world = lcao_factory.world();
  Matrix Eij_F12;

  utility::print_par(world, "\n Computing CCSD_F12 C Approach \n");

  auto n_active_occ = this->trange1_engine()->get_active_occ();

  // create shape
  auto occ_tr1 = this->trange1_engine()->get_active_occ_tr1();
  TiledArray::TiledRange occ4_trange({occ_tr1, occ_tr1, occ_tr1, occ_tr1});
  auto ijij_ijji_shape = f12::make_ijij_ijji_shape(occ4_trange);

  {
    // compute B term
    TArray B_ijij_ijji = f12::compute_B_ijij_ijji_C(lcao_factory, ijij_ijji_shape);
    Matrix eij = B_ijij_ijji("i1,j1,i2,j2")
                     .reduce(f12::F12PairEnergyReductor<Tile>(
                         f12::CC_ijij_bar, f12::CC_ijji_bar, n_active_occ));
    if (debug()) utility::print_par(world, "E_B: ", eij.sum(), "\n");
    Eij_F12 = eij;
  }

  // compute X term
  TArray X_ijij_ijji = f12::compute_X_ijij_ijji(lcao_factory, ijij_ijji_shape);
  //    std::cout << "X_ijij_ijji" << std::endl;
  //    std::cout << X_ijij_ijji << std::endl;
  lcao_factory.purge_operator(world, L"R2");

  auto Fij = lcao_factory.compute(L"<i|F|j>");
  auto Fij_eigen = array_ops::array_to_eigen(Fij);
  f12::convert_X_ijkl(X_ijij_ijji, Fij_eigen);
  {
    Matrix eij = X_ijij_ijji("i1,j1,i2,j2")
                     .reduce(f12::F12PairEnergyReductor<Tile>(
                         f12::CC_ijij_bar, f12::CC_ijji_bar, n_active_occ));
    eij *= -1;
    if (debug()) utility::print_par(world, "E_X: ", eij.sum(), "\n");
    Eij_F12 += eij;
  }

  // compute V_ijij_ijji
  TArray V_ijij_ijji = f12::compute_V_ijij_ijji(lcao_factory, ijij_ijji_shape);

  //    std::cout << "V_ijij_ijji" << std::endl;
  //    std::cout << V_ijij_ijji << std::endl;

  // VT2 contribution
  //    if(darray.is_initialized()){
  //        TArray tmp = f12::compute_VT2_ijij_ijji_df_direct(lcao_factory,
  //        this->t2(), ijij_ijji_shape, darray);
  //        V_ijij_ijji("i1,j1,i2,j2") += tmp("i1,j1,i2,j2");
  //    }else
  {
    TArray tmp = f12::compute_VT2_ijij_ijji(lcao_factory, this->t2(),
                                       ijij_ijji_shape, vt_couple_);
    V_ijij_ijji("i1,j1,i2,j2") += tmp("i1,j1,i2,j2");
  }

  // VT1 contribution
  {
    TArray tmp = f12::compute_VT1_ijij_ijji(lcao_factory, this->t1(),
                                       ijij_ijji_shape, vt_couple_);
    V_ijij_ijji("i1,j1,i2,j2") += tmp("i1,j1,i2,j2");
  }

  // V contribution to energy
  Matrix eij = V_ijij_ijji("i1,j1,i2,j2")
                   .reduce(f12::F12PairEnergyReductor<Tile>(
                       2 * f12::C_ijij_bar, 2 * f12::C_ijji_bar, n_active_occ));
  if (debug()) utility::print_par(world, "E_V: ", eij.sum(), "\n");
  Eij_F12 += eij;

  return Eij_F12;
}

#if TA_DEFAULT_POLICY == 1
extern template class CCSD_F12<TA::TensorD>;
#endif

}  // namespace lcao
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_F12_CCSD_F12_H_
