//
// Created by Chong Peng on 4/12/16.
//

#ifndef MPQC_CCSDF12_H
#define MPQC_CCSDF12_H

#include "../../../../../common/namespaces.h"
#include "../../../../../include/eigen.h"
#include "../../../../../include/tiledarray.h"
#include "../../../../../utility/trange1_engine.h"
#include <mpqc/chemistry/qc/cc/ccsd.h>
#include <mpqc/chemistry/qc/f12/f12_intermediates.h>

#include <mpqc/chemistry/qc/f12/cabs_singles.h>

namespace mpqc {
namespace f12 {

/**
 *  CCSD(2)F12 Takes all options from CCSD
 *
// *  @param ///MP2F12: bool, default false
 *  @param Singles, bool, if compute cabs singles correction, default true
 *  @param VTCouple, bool, if include <p q |G| m a'> <p q|R|m a'> term in VT coupling , default true
 *  @param Approach = string, use C or D approach, default is C
 */

template <typename Tile>
class CCSDF12 {
 public:
  using Policy = TA::SparsePolicy;
  using TArray = TA::DistArray<Tile, Policy>;
  using LCAOFactoryType = integrals::LCAOFactory<Tile, Policy>;

  using real_t = typename Tile::scalar_type;
  using Matrix = RowMatrix<real_t>;

  CCSDF12(integrals::LCAOFactory<Tile, Policy>& lcao_factory,
          rapidjson::Document& options)
      : lcao_factory_(lcao_factory) {
    ccsd_ = std::make_shared<cc::CCSD<Tile, Policy>>(lcao_factory, options);
    vt_couple_ = ccsd_->options().HasMember("VTCouple") ? ccsd_->options()["VTCouple"].GetBool() : true;
  }

  CCSDF12(std::shared_ptr<cc::CCSD<Tile, Policy>> ccsd)
      : lcao_factory_(ccsd->intermediate()->lcao_factory()), ccsd_(ccsd) {}

  virtual real_t compute() {
    auto& world = lcao_factory_.get_world();
    // compute ccsd
    real_t ccsd = ccsd_->compute();

    auto& option = ccsd_->options();
    // initialize CABS orbitals
    closed_shell_cabs_mo_build_svd(this->lcao_factory_, option,
                                   this->ccsd_->trange1_engine());

    auto lazy_two_electron_int = ccsd_->intermediate()->direct_ao();
    Matrix Eij_F12;

    std::string method =
        option.HasMember("Method") ? option["Method"].GetString() : "df";

    std::string approach =
        option.HasMember("Approach") ? option["Approach"].GetString() : "C";

    if (approach != "C" && approach != "D") {
      throw std::runtime_error("Wrong CCSDF12 Approach");
    }

    utility::print_par(lcao_factory_.get_world(), "VTCouple: ", vt_couple_);

    if (method == "four center") {
      Eij_F12 = compute_ccsd_f12(lazy_two_electron_int);
    } else if (method == "df") {
      Eij_F12 = compute_ccsd_f12_df(lazy_two_electron_int, approach);
    } else {
      throw std::runtime_error("Wrong CCSDF12 Method");
    }

    real_t e_f12 = Eij_F12.sum();
    if (debug()) utility::print_par(world, "E_F12: ", e_f12, "\n");

    // compute cabs singles
    real_t e_s = 0.0;
    bool singles =
        option.HasMember("Singles") ? option["Singles"].GetBool() : true;
    if (singles) {
      auto single_time0 = mpqc_time::fenced_now(world);

      CABSSingles<Tile> cabs_singles(lcao_factory_);
      e_s = cabs_singles.compute();
      if (debug()) {
        utility::print_par(world, "E_S: ", e_s, "\n");
      }
      auto single_time1 = mpqc_time::fenced_now(world);
      auto single_time = mpqc_time::duration_in_s(single_time0, single_time1);
      mpqc::utility::print_par(world, "Total CABS Singles Time:  ", single_time,
                               "\n");
    }

    return ccsd + e_f12 + e_s;
  }

 private:
  /// standard approach
  template <typename DirectArray>
  Matrix compute_ccsd_f12_df(const DirectArray& darray,
                             const std::string& approach);

  template <typename DirectArray>
  Matrix compute_ccsd_f12(const DirectArray& darray);

 protected:
  LCAOFactoryType& lcao_factory_;
  std::shared_ptr<cc::CCSD<Tile, Policy>> ccsd_;
  bool vt_couple_;

  int debug() const { return 1; }
};

template <typename Tile>
template <typename DirectArray>
typename CCSDF12<Tile>::Matrix CCSDF12<Tile>::compute_ccsd_f12_df(
    const DirectArray& darray, const std::string& approach) {
  auto& lcao_factory = lcao_factory_;
  auto& world = lcao_factory.get_world();
  Matrix Eij_F12;

  utility::print_par(world, "\n Computing CCSDF12 ", approach, " Approach \n");

  // clean LCAO Integrals
  lcao_factory.registry().clear();

  auto n_active_occ = ccsd_->trange1_engine()->get_active_occ();

  // create shape
  auto occ_tr1 = ccsd_->trange1_engine()->get_occ_tr1();
  TiledArray::TiledRange occ4_trange({occ_tr1, occ_tr1, occ_tr1, occ_tr1});
  auto ijij_ijji_shape = f12::make_ijij_ijji_shape(occ4_trange);

  // compute V_ijij_ijji
  TArray V_ijij_ijji = compute_V_ijij_ijji_df(lcao_factory, ijij_ijji_shape);

  // VT2 contribution
  if (darray.is_initialized()) {
    TArray tmp = compute_VT2_ijij_ijji_df_direct(lcao_factory, ccsd_->t2(),
                                                 ijij_ijji_shape, darray);
    V_ijij_ijji("i1,j1,i2,j2") += tmp("i1,j1,i2,j2");
  } else {
    TArray tmp =
        compute_VT2_ijij_ijji_df(lcao_factory, ccsd_->t2(), ijij_ijji_shape, vt_couple_);
    V_ijij_ijji("i1,j1,i2,j2") += tmp("i1,j1,i2,j2");
  }

  // VT1 contribution
  {
    TArray tmp =
        compute_VT1_ijij_ijji_df(lcao_factory, ccsd_->t1(), ijij_ijji_shape,vt_couple_);
    V_ijij_ijji("i1,j1,i2,j2") += tmp("i1,j1,i2,j2");
  }

  // V contribution to energy
  Eij_F12 = V_ijij_ijji("i1,j1,i2,j2")
                .reduce(f12::F12PairEnergyReductor<Tile>(
                    2 * C_ijij_bar, 2 * C_ijji_bar, n_active_occ));
  if (debug()) utility::print_par(world, "E_V: ", Eij_F12.sum(), "\n");

  // compute X term
  TArray X_ijij_ijji = compute_X_ijij_ijji_df(lcao_factory, ijij_ijji_shape);
  // R_ipjq not needed
  lcao_factory_.registry().purge_formula(world, L"(i1 p|R|j1 q)[df]");

  auto Fij = lcao_factory_.compute(L"(i|F|j)[df]");
  auto Fij_eigen = array_ops::array_to_eigen(Fij);
  f12::convert_X_ijkl(X_ijij_ijji, Fij_eigen);
  {
    Matrix eij = X_ijij_ijji("i1,j1,i2,j2")
                     .reduce(f12::F12PairEnergyReductor<Tile>(
                         CC_ijij_bar, CC_ijji_bar, n_active_occ));
    eij *= -1;
    if (debug()) utility::print_par(world, "E_X: ", eij.sum(), "\n");
    Eij_F12 += eij;
  }

  // compute B term
  {
    TArray B_ijij_ijji;

    if (approach == "C") {
      B_ijij_ijji = compute_B_ijij_ijji_df(lcao_factory, ijij_ijji_shape);
    } else if (approach == "D") {
      B_ijij_ijji = compute_B_ijij_ijji_D_df(lcao_factory, ijij_ijji_shape);
    }

    Matrix eij = B_ijij_ijji("i1,j1,i2,j2")
                     .reduce(F12PairEnergyReductor<Tile>(
                         CC_ijij_bar, CC_ijji_bar, n_active_occ));
    if (debug()) utility::print_par(world, "E_B: ", eij.sum(), "\n");
    Eij_F12 += eij;
  }

  return Eij_F12;
}

template <typename Tile>
template <typename DirectArray>
typename CCSDF12<Tile>::Matrix CCSDF12<Tile>::compute_ccsd_f12(
    const DirectArray& darray) {
  auto& lcao_factory = lcao_factory_;
  auto& world = lcao_factory.get_world();
  Matrix Eij_F12;

  utility::print_par(world, "\n Computing CCSDF12 C Approach \n");

  // clean LCAO Integrals
  lcao_factory.registry().clear();

  auto n_active_occ = ccsd_->trange1_engine()->get_active_occ();

  // create shape
  auto occ_tr1 = ccsd_->trange1_engine()->get_occ_tr1();
  TiledArray::TiledRange occ4_trange({occ_tr1, occ_tr1, occ_tr1, occ_tr1});
  auto ijij_ijji_shape = f12::make_ijij_ijji_shape(occ4_trange);

  // compute V_ijij_ijji
  TArray V_ijij_ijji = compute_V_ijij_ijji(lcao_factory, ijij_ijji_shape);

  //    std::cout << "V_ijij_ijji" << std::endl;
  //    std::cout << V_ijij_ijji << std::endl;

  // VT2 contribution
  //    if(darray.is_initialized()){
  //        TArray tmp = compute_VT2_ijij_ijji_df_direct(lcao_factory,
  //        ccsd_->t2(), ijij_ijji_shape, darray);
  //        V_ijij_ijji("i1,j1,i2,j2") += tmp("i1,j1,i2,j2");
  //    }else
  {
    TArray tmp =
        compute_VT2_ijij_ijji(lcao_factory, ccsd_->t2(), ijij_ijji_shape, vt_couple_);
    V_ijij_ijji("i1,j1,i2,j2") += tmp("i1,j1,i2,j2");
  }

  // VT1 contribution
  {
    TArray tmp =
        compute_VT1_ijij_ijji(lcao_factory, ccsd_->t1(), ijij_ijji_shape,vt_couple_);
    V_ijij_ijji("i1,j1,i2,j2") += tmp("i1,j1,i2,j2");
  }

  // V contribution to energy
  Eij_F12 = V_ijij_ijji("i1,j1,i2,j2")
                .reduce(f12::F12PairEnergyReductor<Tile>(
                    2 * C_ijij_bar, 2 * C_ijji_bar, n_active_occ));
  if (debug()) utility::print_par(world, "E_V: ", Eij_F12.sum(), "\n");

  //    {
  //        utility::print_par(world, "Compute CC Term Without DF \n");
  //        auto C_ijab = compute_C_ijab(lcao_factory);
  //        auto C_bar_ijab = f12::convert_C_ijab(C_ijab, occ,
  //        *orbital_energy_);
  //        V_ijij_ijji("i1,j1,i2,j2") =
  //        (C_ijab("i1,j1,a,b")*C_bar_ijab("i2,j2,a,b")).set_shape(ijij_ijji_shape);
  //
  //        double E_cc =
  //        V_ijij_ijji("i1,j1,i2,j2").reduce(f12::CLF12Energy<Tile>(CC_ijij_bar,CC_ijji_bar));
  //        utility::print_par(world, "E_CC: ", E_cc, "\n");
  //        E += E_cc;
  //    }
  // compute X term
  TArray X_ijij_ijji = compute_X_ijij_ijji(lcao_factory, ijij_ijji_shape);
  //    std::cout << "X_ijij_ijji" << std::endl;
  //    std::cout << X_ijij_ijji << std::endl;

  // R_ipjq not needed
  lcao_factory_.registry().purge_formula(world, L"(i1 p|R|j1 q)");

  auto Fij = lcao_factory_.compute(L"(i|F|j)");
  auto Fij_eigen = array_ops::array_to_eigen(Fij);
  f12::convert_X_ijkl(X_ijij_ijji, Fij_eigen);
  {
    Matrix eij = X_ijij_ijji("i1,j1,i2,j2")
                     .reduce(f12::F12PairEnergyReductor<Tile>(
                         CC_ijij_bar, CC_ijji_bar, n_active_occ));
    eij *= -1;
    if (debug()) utility::print_par(world, "E_X: ", eij.sum(), "\n");
    Eij_F12 += eij;
  }

  {
    // compute B term
    TArray B_ijij_ijji = compute_B_ijij_ijji(lcao_factory, ijij_ijji_shape);
    Matrix eij = B_ijij_ijji("i1,j1,i2,j2")
                     .reduce(F12PairEnergyReductor<Tile>(
                         CC_ijij_bar, CC_ijji_bar, n_active_occ));
    if (debug()) utility::print_par(world, "E_B: ", eij.sum(), "\n");
    Eij_F12 += eij;
  }

  return Eij_F12;
}

}  // end of namespace f12
}  // end of namespace mpqc

#endif  // MPQC_CCSDF12_H
