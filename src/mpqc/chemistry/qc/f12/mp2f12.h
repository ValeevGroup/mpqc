//
// Created by Chong Peng on 3/31/16.
//

#ifndef MPQC_MP2F12_H
#define MPQC_MP2F12_H

#include <string>

#include "../../../../../include/eigen.h"
#include "../../../../../utility/cc_utility.h"
#include "../../../../../utility/trange1_engine.h"
#include <mpqc/chemistry/qc/f12/cabs_singles.h>
#include <mpqc/chemistry/qc/f12/f12_intermediates.h>
#include <mpqc/chemistry/qc/f12/f12_utility.h>
#include <mpqc/chemistry/qc/mbpt/mp2.h>

namespace mpqc {
namespace f12 {

/**
 *  MP2F12 Class
 *  Take all the options from MP2
 *
 *  Other options:
 *  @param Singles = bool, if perform cabs_singles calculation, default is true
 *  @param Approach = string, use C or D approach, default is C
 *
 */

template <typename Tile>
class MP2F12 {
 public:
  using Policy = TA::SparsePolicy;
  using TArray = TA::DistArray<Tile, Policy>;
  using LCAOFactoryType = integrals::LCAOFactory<Tile, Policy>;

  using real_t = typename Tile::scalar_type;
  using Matrix = RowMatrix<real_t>;

  MP2F12() = default;

  MP2F12(LCAOFactoryType& lcao_factory) {
    mp2_ = std::make_shared<mbpt::MP2<Tile, Policy>>(lcao_factory);
  }

  MP2F12(std::shared_ptr<mbpt::MP2<Tile, Policy>> mp2) : mp2_(mp2) {}

  LCAOFactoryType& lcao_factory() const { return mp2_->lcao_factory(); }

  const std::shared_ptr<TRange1Engine> trange1_engine() const {
    return mp2_->trange1_engine();
  }

  const std::shared_ptr<Eigen::VectorXd> orbital_energy() const {
    return mp2_->orbital_energy();
  }

  std::tuple<Matrix, Matrix> compute_mp2_f12_df(const std::string& approach);

  std::tuple<Matrix, Matrix> compute_mp2_f12();

  real_t compute(const rapidjson::Document& in) {
    auto& world = lcao_factory().get_world();
    auto f12_time0 = mpqc_time::fenced_now(world);
    mp2_->init(in);

    // solve cabs orbitals
    auto orbital_registry = lcao_factory().orbital_space();
    closed_shell_cabs_mo_build_svd(lcao_factory(), in,
                                   this->mp2_->trange1_engine());

    std::string method =
        in.HasMember("Method") ? in["Method"].GetString() : "df";

    std::string approach = in.HasMember("Approach") ? in["Approach"].GetString() : "C";

    if(approach!="C" && approach != "D"){
      throw std::runtime_error("Wrong MP2F12 Approach");
    }

    Matrix mp2_eij, f12_eij;

    lcao_factory().registry().purge(world);

    bool df;
    if (method == "four center") {
      std::tie(mp2_eij, f12_eij) = compute_mp2_f12();
      df = false;
    } else if (method == "df") {
      std::tie(mp2_eij, f12_eij) = compute_mp2_f12_df(approach);
      df = true;
    } else {
      throw std::runtime_error("Wrong MP2F12 Method");
    }

    if (world.rank() == 0) {
      auto nocc = mp2_->trange1_engine()->get_active_occ();
      printf(
          "  i0     i1       eij(mp2)        eij(f12)      eij(mp2-f12) \n"
          "====== ====== =============== =============== ===============\n");
      for (int i = 0; i != nocc; ++i)
        for (int j = i; j != nocc; ++j)
          printf("%4d   %4d   %15.12lf %15.12lf %15.12lf\n", i, j,
                 mp2_eij(i, j), f12_eij(i, j), mp2_eij(i, j) + f12_eij(i, j));
    }

    auto emp2 = mp2_eij.sum();
    auto ef12 = f12_eij.sum();
    if (debug()) {
      utility::print_par(lcao_factory().get_world(), "E_MP2: ", emp2, "\n");
      utility::print_par(lcao_factory().get_world(), "E_F12: ", ef12, "\n");
    }

    // compute cabs singles
    real_t e_s = 0.0;
    bool singles = in.HasMember("Singles") ? in["Singles"].GetBool() : true;
    if (singles) {
      auto single_time0 = mpqc_time::fenced_now(world);

      if(approach=="D"){
        CABSSingles<Tile> cabs_singles(lcao_factory());
        e_s = cabs_singles.compute(df,true,true);
      }else{
        CABSSingles<Tile> cabs_singles(lcao_factory());
        e_s = cabs_singles.compute(df,false,true);
      }

      if (debug()) {
        utility::print_par(lcao_factory().get_world(), "E_S: ", e_s, "\n");
      }
      auto single_time1 = mpqc_time::fenced_now(world);
      auto single_time = mpqc_time::duration_in_s(single_time0, single_time1);
      mpqc::utility::print_par(world, "Total CABS Singles Time:  ", single_time,
                               "\n");
    }

    auto f12_time1 = mpqc_time::fenced_now(world);
    auto f12_time = mpqc_time::duration_in_s(f12_time0, f12_time1);
    mpqc::utility::print_par(world, "Total MP2F12 Time:  ", f12_time, "\n");

    return emp2 + ef12 + e_s;
  }

 protected:
  std::shared_ptr<mbpt::MP2<Tile, Policy>> mp2_;

  int debug() const { return 1; }
};

template <typename Tile>
std::tuple<typename MP2F12<Tile>::Matrix, typename MP2F12<Tile>::Matrix>
MP2F12<Tile>::compute_mp2_f12_df(const std::string& approach) {

  auto& world = lcao_factory().get_world();

  utility::print_par(world, "\n Computing MP2F12 ", approach ," Approach \n");

  Matrix Eij_MP2, Eij_F12;

  auto n_active_occ = mp2_->trange1_engine()->get_active_occ();
  auto n_occ = mp2_->trange1_engine()->get_occ();
  auto n_frozen = mp2_->trange1_engine()->get_nfrozen();

  // create shape
  auto occ_tr1 = mp2_->trange1_engine()->get_occ_tr1();
  TiledArray::TiledRange occ4_trange({occ_tr1, occ_tr1, occ_tr1, occ_tr1});
  auto ijij_ijji_shape = f12::make_ijij_ijji_shape(occ4_trange);

  // compute B term
  TArray B_ijij_ijji;
  if(approach=="C"){
    B_ijij_ijji = compute_B_ijij_ijji_C_df(lcao_factory(), ijij_ijji_shape);
  }
  else if(approach=="D"){
    B_ijij_ijji = compute_B_ijij_ijji_D_df(lcao_factory(), ijij_ijji_shape);
  }

  {
    Matrix Eij_b = B_ijij_ijji("i1,j1,i2,j2")
        .reduce(F12PairEnergyReductor<Tile>(
            CC_ijij_bar, CC_ijji_bar, n_active_occ));
    if (debug()) utility::print_par(world, "E_B: ", Eij_b.sum(), "\n");
    Eij_F12 = Eij_b;
  }

  // compute X term
  TArray X_ijij_ijji = compute_X_ijij_ijji_df(lcao_factory(), ijij_ijji_shape);
  lcao_factory().purge_operator(world,L"R2");

  {
    auto Fij = lcao_factory().compute(L"<i|F|j>[df]");
    auto Fij_eigen = array_ops::array_to_eigen(Fij);
    f12::convert_X_ijkl(X_ijij_ijji, Fij_eigen);

    Matrix Eij_x = X_ijij_ijji("i1,j1,i2,j2")
        .reduce(F12PairEnergyReductor<Tile>(
            CC_ijij_bar, CC_ijji_bar, n_active_occ));
    Eij_x *= -1.0;
    if (debug()) utility::print_par(world, "E_X: ", Eij_x.sum(), "\n");
    Eij_F12 += Eij_x;
  }

  // compute V term
  TArray V_ijij_ijji = compute_V_ijij_ijji_df(lcao_factory(), ijij_ijji_shape);
  {
    // G integral in MO not needed, still need G integral in AO to compute F, K,
    // hJ
    lcao_factory().registry().purge_operator(world, L"G");
    lcao_factory().purge_operator(world, L"GR");

    // contribution from V_ijij_ijji
    // NB factor of 2 from the Hylleraas functional
    Matrix e_ij = V_ijij_ijji("i1,j1,i2,j2")
        .reduce(F12PairEnergyReductor<Tile>(
            2 * C_ijij_bar, 2 * C_ijji_bar, n_active_occ));
    Eij_F12 += e_ij;
    if (debug()) utility::print_par(world, "E_V: ", e_ij.sum(), "\n");
  }

  TArray t2;
  {
    utility::print_par(world, "Compute T_abij With DF \n");

    TArray g_abij;
    g_abij("a,b,i,j") = lcao_factory().compute(L"<i j|G|a b>[df]")("i,j,a,b");
    t2 = mpqc::cc::d_abij(g_abij, *(mp2_->orbital_energy()), n_occ, n_frozen);

    // compute MP2 energy and pair energies
    TArray TG_ijij_ijji;
    TG_ijij_ijji("i1,j1,i2,j2") =
        (t2("a,b,i1,j1") * g_abij("a,b,i2,j2")).set_shape(ijij_ijji_shape);
    Eij_MP2 = TG_ijij_ijji("i1,j1,i2,j2")
                  .reduce(F12PairEnergyReductor<Tile>(2, -1, n_active_occ));
  }

  // compute C term
  TArray C_ijab = compute_C_ijab_df(lcao_factory());

  {
    utility::print_par(world, "Compute CT With DF \n");
    V_ijij_ijji("i1,j1,i2,j2") =
        (C_ijab("i1,j1,a,b") * t2("a,b,i2,j2")).set_shape(ijij_ijji_shape);

    // NB factor of 2 from the Hylleraas functional
    Matrix Eij_ct = V_ijij_ijji("i1,j1,i2,j2")
                        .reduce(F12PairEnergyReductor<Tile>(
                            2 * C_ijij_bar, 2 * C_ijji_bar, n_active_occ));
    if (debug()) utility::print_par(world, "E_CT: ", Eij_ct.sum(), "\n");
    Eij_F12 += Eij_ct;
  }

  {
    utility::print_par(world, "Compute CC Term With DF \n");
    auto C_bar_ijab =
        f12::convert_C_ijab(C_ijab, n_occ, n_frozen, *(mp2_->orbital_energy()));
    B_ijij_ijji("i1,j1,i2,j2") = (C_ijab("i1,j1,a,b") * C_bar_ijab("i2,j2,a,b"))
        .set_shape(ijij_ijji_shape);

    Matrix Eij_cc = B_ijij_ijji("i1,j1,i2,j2")
        .reduce(F12PairEnergyReductor<Tile>(
            CC_ijij_bar, CC_ijji_bar, n_active_occ));
    if (debug()) utility::print_par(world, "E_CC: ", Eij_cc.sum(), "\n");
    Eij_F12 += Eij_cc;
  }


  return std::make_tuple(Eij_MP2, Eij_F12);
}

template <typename Tile>
std::tuple<typename MP2F12<Tile>::Matrix, typename MP2F12<Tile>::Matrix>
MP2F12<Tile>::compute_mp2_f12() {

  auto& world = lcao_factory().get_world();

  utility::print_par(world, "\n Computing MP2F12 C Approach \n");

  Matrix Eij_MP2, Eij_F12;

  auto n_active_occ = mp2_->trange1_engine()->get_active_occ();
  auto n_occ = mp2_->trange1_engine()->get_occ();
  auto n_frozen = mp2_->trange1_engine()->get_nfrozen();

  // create shape
  auto occ_tr1 = mp2_->trange1_engine()->get_occ_tr1();
  TiledArray::TiledRange occ4_trange({occ_tr1, occ_tr1, occ_tr1, occ_tr1});
  auto ijij_ijji_shape = f12::make_ijij_ijji_shape(occ4_trange);

  TArray B_ijij_ijji_nodf =
      compute_B_ijij_ijji_C(lcao_factory(), ijij_ijji_shape);
  {
    Matrix Eij_b = B_ijij_ijji_nodf("i1,j1,i2,j2")
        .reduce(F12PairEnergyReductor<Tile>(
            CC_ijij_bar, CC_ijji_bar, n_active_occ));
    if (debug()) utility::print_par(world, "E_B: ", Eij_b.sum(), "\n");
    Eij_F12 = Eij_b;
  }

  TArray X_ijij_ijji_nodf =
      compute_X_ijij_ijji(lcao_factory(), ijij_ijji_shape);
  lcao_factory().purge_operator(world,L"R2");
  {
    // compute energy contribution
    auto Fij = lcao_factory().compute(L"<i|F|j>");
    auto Fij_eigen = array_ops::array_to_eigen(Fij);
    f12::convert_X_ijkl(X_ijij_ijji_nodf, Fij_eigen);

    Matrix Eij_x = X_ijij_ijji_nodf("i1,j1,i2,j2")
        .reduce(F12PairEnergyReductor<Tile>(
            CC_ijij_bar, CC_ijji_bar, n_active_occ));
    Eij_x *= -1;
    if (debug()) utility::print_par(world, "E_X: ", Eij_x.sum(), "\n");
    Eij_F12 += Eij_x;
  }


  TArray V_ijij_ijji_nodf =
      compute_V_ijij_ijji(lcao_factory(), ijij_ijji_shape);
  lcao_factory().registry().purge_operator(world, L"G");
  lcao_factory().purge_operator(world, L"GR");
  {
    Matrix e_ij = V_ijij_ijji_nodf("i1,j1,i2,j2")
                  .reduce(F12PairEnergyReductor<Tile>(
                      2 * C_ijij_bar, 2 * C_ijji_bar, n_active_occ));
    Eij_F12 += e_ij;
    if (debug()) utility::print_par(world, "E_V: ", e_ij.sum(), "\n");
  }

  TArray t2_nodf;  // t2_abij
  {
    utility::print_par(world, "Compute T_abij Without DF \n");
    TArray g_abij;
    g_abij("a,b,i,j") = lcao_factory().compute(L"<i j|G|a b>")("i,j,a,b");
    t2_nodf =
        mpqc::cc::d_abij(g_abij, *(mp2_->orbital_energy()), n_occ, n_frozen);
    TArray TG_ijij_ijji_nodf;
    TG_ijij_ijji_nodf("i1,j1,i2,j2") =
        (t2_nodf("a,b,i1,j1") * g_abij("a,b,i2,j2")).set_shape(ijij_ijji_shape);
    Eij_MP2 = TG_ijij_ijji_nodf("i1,j1,i2,j2")
        .reduce(F12PairEnergyReductor<Tile>(2, -1, n_active_occ));
  }

  TArray C_ijab_nodf = compute_C_ijab(lcao_factory());

  {
    utility::print_par(world, "Compute CT Without DF \n");
    V_ijij_ijji_nodf("i1,j1,i2,j2") =
        (C_ijab_nodf("i1,j1,a,b") * t2_nodf("a,b,i2,j2"))
            .set_shape(ijij_ijji_shape);

    Matrix Eij_ct = V_ijij_ijji_nodf("i1,j1,i2,j2")
                        .reduce(F12PairEnergyReductor<Tile>(
                            2 * C_ijij_bar, 2 * C_ijji_bar, n_active_occ));
    if (debug()) utility::print_par(world, "E_CT: ", Eij_ct.sum(), "\n");
    Eij_F12 += Eij_ct;
  }

  {
    utility::print_par(world, "Compute CC Term Without DF \n");
    auto C_bar_ijab = f12::convert_C_ijab(C_ijab_nodf, n_occ, n_frozen,
                                          *(mp2_->orbital_energy()));
    B_ijij_ijji_nodf("i1,j1,i2,j2") =
        (C_ijab_nodf("i1,j1,a,b") * C_bar_ijab("i2,j2,a,b"))
            .set_shape(ijij_ijji_shape);

    Matrix Eij_cc = B_ijij_ijji_nodf("i1,j1,i2,j2")
                        .reduce(F12PairEnergyReductor<Tile>(
                            CC_ijij_bar, CC_ijji_bar, n_active_occ));
    if (debug()) utility::print_par(world, "E_CC: ", Eij_cc.sum(), "\n");
    Eij_F12 += Eij_cc;
  }

  return std::make_tuple(Eij_MP2, Eij_F12);
}

class RMP2F12 : public qc::LCAOWavefunction{

public:
  using TArray = qc::LCAOWavefunction::ArrayType;
  using Matrix = RowMatrix<double>;

  RMP2F12(const KeyVal& kv);
  ~RMP2F12() = default;

  double value() override;
  std::tuple<Matrix,Matrix> compute();
  void compute(qc::PropertyBase* pb) override;
  void obsolete() override;

private:

  virtual TArray compute_B();
  virtual TArray compute_V();
  virtual TArray compute_X();
  virtual TArray compute_C();
  virtual std::tuple<TArray, TArray> compute_T();
  virtual double compute_cabs_singles();

protected:
  char approximation_;
  TA::SparseShape<float> ijij_ijji_shape_;

private:
  std::shared_ptr<qc::Wavefunction> ref_wfn_;
  double rmp2f12_energy_;
  bool cabs_singles_;
};

class RIRMP2F12 : public RMP2F12{

public:
  RIRMP2F12(const KeyVal& kv);
  ~RIRMP2F12() = default;

private:
  TArray compute_B() override ;
  TArray compute_V() override ;
  TArray compute_X() override ;
  TArray compute_C() override ;
  std::tuple<TArray, TArray> compute_T() override ;
  double compute_cabs_singles() override;

};

}  // end of namespace f12
}  // mpqc

#endif  // MPQC_MP2F12_H
