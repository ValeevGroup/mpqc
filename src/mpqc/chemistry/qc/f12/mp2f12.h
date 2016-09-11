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

/*
 *  MP2F12 Class
 *  Take all the options from MP2
 *
 *  Other options:
 *  Singles = bool, if perform cabs_singles calculation, default is true
 *  Approach = string, use C or D approach, default is C
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

    if (method == "four center") {
      std::tie(mp2_eij, f12_eij) = compute_mp2_f12();
    } else if (method == "df") {
      std::tie(mp2_eij, f12_eij) = compute_mp2_f12_df(approach);
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

      CABSSingles<Tile> cabs_singles(lcao_factory());
      e_s = cabs_singles.compute();
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

  // compute V term
  TArray V_ijij_ijji = compute_V_ijij_ijji_df(lcao_factory(), ijij_ijji_shape);
  {
    // G integral in MO not needed, still need G integral in AO to compute F, K,
    // hJ
    lcao_factory().registry().purge_operator(world, L"G");

    // contribution from V_ijij_ijji
    // NB factor of 2 from the Hylleraas functional
    Eij_F12 = V_ijij_ijji("i1,j1,i2,j2")
                  .reduce(F12PairEnergyReductor<Tile>(
                      2 * C_ijij_bar, 2 * C_ijji_bar, n_active_occ));
    if (debug()) utility::print_par(world, "E_V: ", Eij_F12.sum(), "\n");
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

  // compute X term
  TArray X_ijij_ijji = compute_X_ijij_ijji_df(lcao_factory(), ijij_ijji_shape);
  {
    // R_ipjq not needed
    lcao_factory().registry().purge_formula(world, L"<i1 j1|R|p q>[df]");

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

  // compute B term
  TArray B_ijij_ijji;
  if(approach=="C"){
    B_ijij_ijji = compute_B_ijij_ijji_df(lcao_factory(), ijij_ijji_shape);
  }
  else if(approach=="D"){
    B_ijij_ijji = compute_B_ijij_ijji_D_df(lcao_factory(), ijij_ijji_shape);
  }

  {
    Matrix Eij_b = B_ijij_ijji("i1,j1,i2,j2")
                       .reduce(F12PairEnergyReductor<Tile>(
                           CC_ijij_bar, CC_ijji_bar, n_active_occ));
    if (debug()) utility::print_par(world, "E_B: ", Eij_b.sum(), "\n");
    Eij_F12 += Eij_b;
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

  TArray V_ijij_ijji_nodf =
      compute_V_ijij_ijji(lcao_factory(), ijij_ijji_shape);
  {
    Eij_F12 = V_ijij_ijji_nodf("i1,j1,i2,j2")
                  .reduce(F12PairEnergyReductor<Tile>(
                      2 * C_ijij_bar, 2 * C_ijji_bar, n_active_occ));
    if (debug()) utility::print_par(world, "E_V: ", Eij_F12.sum(), "\n");
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

  TArray X_ijij_ijji_nodf =
      compute_X_ijij_ijji(lcao_factory(), ijij_ijji_shape);
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

  TArray B_ijij_ijji_nodf =
      compute_B_ijij_ijji(lcao_factory(), ijij_ijji_shape);
  {
    Matrix Eij_b = B_ijij_ijji_nodf("i1,j1,i2,j2")
                       .reduce(F12PairEnergyReductor<Tile>(
                           CC_ijij_bar, CC_ijji_bar, n_active_occ));
    if (debug()) utility::print_par(world, "E_B: ", Eij_b.sum(), "\n");
    Eij_F12 += Eij_b;
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

}  // end of namespace f12

namespace dyson {

enum class Denominator {
  pqrE,  // ( e(p) + e(q) - e(r) - E)
  rEpq   // (-e(p) - e(q) + e(r) + E)
};

template <Denominator Denom, typename Tile, typename Policy>
TA::Array<double, 4, Tile, Policy> d_pqrE(
    TA::Array<double, 4, Tile, Policy>& pqrs, const Eigen::VectorXd& evals_p,
    const Eigen::VectorXd& evals_q, const Eigen::VectorXd& evals_r,
    typename Tile::scalar_type E) {
  auto convert = [evals_p, evals_q, evals_r, E](Tile& result_tile,
                                                const Tile& arg_tile) {

    result_tile = Tile(arg_tile.range());

    // compute index
    const auto p0 = result_tile.range().lobound()[0];
    const auto p1 = result_tile.range().upbound()[0];
    const auto q0 = result_tile.range().lobound()[1];
    const auto q1 = result_tile.range().upbound()[1];
    const auto r0 = result_tile.range().lobound()[2];
    const auto r1 = result_tile.range().upbound()[2];
    const auto s0 = result_tile.range().lobound()[3];
    const auto s1 = result_tile.range().upbound()[3];

    auto tile_idx = 0;
    typename Tile::value_type norm2 = 0.0;
    for (auto p = p0; p != p1; ++p) {
      const auto e_p = evals_p[p];
      for (auto q = q0; q != q1; ++q) {
        const auto e_q = evals_q[q];
        for (auto r = r0; r < r1; ++r) {
          const auto e_r = evals_r[r];
          const auto denom =
              1 / ((Denom == Denominator::pqrE) ? -e_r - E + e_p + e_q
                                                : e_r + E - e_p - e_q);
          for (auto s = s0; s < s1; ++s, ++tile_idx) {
            const auto v_pqrs = arg_tile[tile_idx];
            const auto vv_pqrs = v_pqrs * denom;
            norm2 += vv_pqrs * vv_pqrs;
            result_tile[tile_idx] = vv_pqrs;
          }
        }
      }
    }
    return std::sqrt(norm2);
  };

  return TA::foreach (pqrs, convert);
}

}  // namespace dyson

namespace f12 {

template <typename Tile>
class GF2F12 {
 public:
  using Policy = TA::SparsePolicy;
  using TArray = TA::DistArray<Tile, Policy>;
  using LCAOFactoryType = integrals::LCAOFactory<Tile, Policy>;

  using real_t = typename Tile::scalar_type;
  using Matrix = RowMatrix<real_t>;

  GF2F12() = default;

  GF2F12(LCAOFactoryType& lcao_factory)
      : mp2_(std::make_shared<mbpt::MP2<Tile, Policy>>(lcao_factory)) {}

  LCAOFactoryType& lcao_factory() const { return mp2_->lcao_factory(); }

  const std::shared_ptr<TRange1Engine> trange1_engine() const {
    return mp2_->trange1_engine();
  }

  const std::shared_ptr<Eigen::VectorXd> orbital_energy() const {
    return mp2_->orbital_energy();
  }

  virtual real_t compute(const rapidjson::Document& in) {
    using mpqc::utility::print_par;

    auto& world = this->lcao_factory().get_world();

    this->mp2_->init(in);

    // compute cabs
    auto orbital_registry = lcao_factory().orbital_space();
    closed_shell_cabs_mo_build_svd(lcao_factory(), in,
                                   this->mp2_->trange1_engine());

    auto time0 = mpqc_time::fenced_now(world);

    orbital_ = in.HasMember("Orbital") ? in["Orbital"].GetInt() : -1;
    use_cabs_ = in.HasMember("UseCABS") ? in["UseCABS"].GetBool() : true;
    if (orbital_ == 0)
      throw std::runtime_error(
          "GF2F12::Orbital must be positive (for particles) or negative (for "
          "holes)");
    if (orbital_ < 0 &&
        (abs(orbital_) - 1 >= trange1_engine()->get_active_occ()))
      throw std::runtime_error(
          "GF2F12::orbital is invalid (the number of holes exceeded)");
    if (orbital_ > 0 && (orbital_ - 1 >= trange1_engine()->get_vir()))
      throw std::runtime_error(
          "GF2F12::orbital is invalid (the number of particles exceeded)");

    std::string method_str = in.HasMember("DysonMethod")
                                 ? in["DysonMethod"].GetString()
                                 : "diagonal";
    enum class Method { diag, nondiag };
    Method method;
    if (method_str.find("diagonal") == 0)
      method = Method::diag;
    else if (method_str == "nondiagonal")
      method = Method::nondiag;
    else
      throw std::runtime_error("GF2F12: unknown value for keyword \"method\"");

    print_par(world, "orbital = ", orbital_, " method = ", method_str,
              " cabs = ", std::to_string(use_cabs_), "\n");

    if (method == Method::diag)
      compute_diagonal();
    else
      compute_nondiagonal();

    auto time1 = mpqc_time::fenced_now(world);
    auto time = mpqc_time::duration_in_s(time0, time1);

    print_par(world, "Total GF2F12 Time:  ", time, "\n");

    return 0.0;
  }

 private:
  std::shared_ptr<mbpt::MP2<Tile, Policy>> mp2_;
  int orbital_;
  bool use_cabs_;

  /// use self-energy in diagonal representation
  void compute_diagonal(int max_niter = 100);
  /// use non-diagonal self-energy
  void compute_nondiagonal(int max_niter = 100);
};

template <typename Tile>
void GF2F12<Tile>::compute_diagonal(int max_niter) {
  auto& world = lcao_factory().get_world();

  auto nfzc = this->mp2_->trange1_engine()->get_nfrozen();
  auto nocc = this->mp2_->trange1_engine()->get_active_occ();
  auto nuocc = this->mp2_->trange1_engine()->get_vir();
  // map orbital_ (index relative to Fermi level) to orbital index in active
  // orbital space
  const auto orbital = nfzc + nocc + ((orbital_ < 0) ? orbital_ : orbital_ - 1);
  auto SE = orbital_energy()->operator()(orbital);

  Eigen::VectorXd occ_evals = orbital_energy()->segment(nfzc, nocc + nfzc);
  Eigen::VectorXd uocc_evals = orbital_energy()->segment(nfzc + nocc, nuocc);

  // will use only the target orbital to transform ints
  // create an OrbitalSpace here
  {
    auto& world = this->lcao_factory().get_world();
    auto orbital_registry = this->lcao_factory().orbital_space();
    auto p_space = orbital_registry->retrieve(OrbitalIndex(L"p"));
    auto C_p = array_ops::array_to_eigen(p_space.array());
    auto C_x = C_p.block(0, orbital, C_p.rows(), 1);
    auto tr_obs = p_space.array().trange().data().front();
    TA::TiledRange1 tr_x{0, 1};
    auto C_x_ta = array_ops::eigen_to_array<Tile>(world, C_x, tr_obs, tr_x);

    using OrbitalSpaceTArray = OrbitalSpace<TA::DistArray<Tile, Policy>>;
    auto x_space =
        OrbitalSpaceTArray(OrbitalIndex(L"x"), OrbitalIndex(L"κ"), C_x_ta);
    orbital_registry->add(x_space);
  }

  lcao_factory().keep_partial_transforms(true);

  TArray g_vvog = lcao_factory().compute(L"<a b|G|i x>[df]");
  TArray g_oovg = lcao_factory().compute(L"<i j|G|a x>[df]");

  if (world.rank() == 0) {
    printf("Iter     SE2(in)     SE2(out)   SE2(delta)\n");
    printf("==== =========== =========== ===========\n");
  }
  size_t iter = 0;
  decltype(SE) SE_diff;
  do {
    TArray Sigma_pph;
    {
      TArray dg_vvog = mpqc::dyson::d_pqrE<dyson::Denominator::rEpq>(
          g_vvog, uocc_evals, uocc_evals, occ_evals, SE);
      Sigma_pph("x,y") = 0.5 * (4 * g_vvog("a,b,i,x") - 2 * g_vvog("b,a,i,x")) *
                         dg_vvog("a,b,i,y");
    }
    // std::cout << "Sigma_pph:\n" << Sigma_pph << std::endl;

    TArray Sigma_hhp;
    {
      TArray dg_oovg = mpqc::dyson::d_pqrE<dyson::Denominator::rEpq>(
          g_oovg, occ_evals, occ_evals, uocc_evals, SE);
      Sigma_hhp("x,y") = 0.5 * (4 * g_oovg("i,j,a,x") - 2 * g_oovg("j,i,a,x")) *
                         dg_oovg("i,j,a,y");
    }
    // std::cout << "Sigma_hhp:\n" << Sigma_hhp << std::endl;

    RowMatrixXd Sigma;
    {
      TArray Sigma_ta;
      Sigma_ta("x,y") = Sigma_pph("x,y") + Sigma_hhp("x,y");
      Sigma = array_ops::array_to_eigen(Sigma_ta);
    }

    // std::cout << "Sigma = " << Sigma << std::endl;

    auto SE_updated = Sigma(0, 0) + orbital_energy()->operator()(orbital);
    SE_diff = SE_updated - SE;

    if (world.rank() == 0)
      printf(" %3ld %10.4lf %10.4lf %10.4lf\n", iter, SE, SE_updated, SE_diff);

    SE = SE_updated;
    ++iter;

  } while ((fabs(SE_diff) > 1e-6) && (iter <= max_niter));

  lcao_factory().purge_formula(world, L"<a b|G|i x>[df]");
  lcao_factory().purge_formula(world, L"<i j|G|a x>[df]");

  // for now the f12 contribution is energy-independent
  TArray Sigma_pph_f12;
  {
    TArray V_ixix, V_ixxi;
    const bool df = true;  // always do DF
    std::tie(V_ixix, V_ixxi) = mpqc::f12::VX_pqrs_pqsr(
        "V", lcao_factory(), "i", "x", "j", "y", df, use_cabs_);
    Sigma_pph_f12("x,y") =
        0.5 * (1.25 * V_ixix("i,x,j,y") - 0.25 * V_ixxi("i,x,y,j")) *
        lcao_factory()(L"<i|I|j>");
  }
  // std::cout << "Sigma_pph_f12:\n" << Sigma_pph_f12 << std::endl;
  Matrix Sigma_f12 = array_ops::array_to_eigen(Sigma_pph_f12);

  // done with F12 ... remove all geminal ints and CABS indices
  lcao_factory().purge_operator(world, L"R");
  lcao_factory().purge_index(world, L"a'");
  lcao_factory().purge_index(world, L"ρ");

  lcao_factory().keep_partial_transforms(false);

  if (world.rank() == 0) {
    auto SE_F12 = SE + Sigma_f12(0, 0);
    auto Hartree2eV = 27.21138602;
    std::string orblabel =
        std::string(orbital_ < 0 ? "IP" : "EA") + std::to_string(abs(orbital_));
    printf("final       GF2 %6s = %11.3lf eV (%10.4lf a.u.)\n",
           orblabel.c_str(), SE * Hartree2eV, SE);
    printf("final GF2-F12-V %6s = %11.3lf eV (%10.4lf a.u.)\n",
           orblabel.c_str(), SE_F12 * Hartree2eV, SE_F12);
    if (orbital_ > 0)
      printf(
          "WARNING: non-strongly-orthogonal F12 projector is used for the F12 "
          "correction to EA!!!");
  }
}

template <typename Tile>
void GF2F12<Tile>::compute_nondiagonal(int max_niter) {
  auto& world = lcao_factory().get_world();

  auto nfzc = this->mp2_->trange1_engine()->get_nfrozen();
  auto nocc = this->mp2_->trange1_engine()->get_active_occ();
  auto nuocc = this->mp2_->trange1_engine()->get_vir();
  // map orbital_ (index relative to Fermi level) to orbital index in active
  // orbital space
  const auto orbital = nfzc + nocc + ((orbital_ < 0) ? orbital_ : orbital_ - 1);
  auto SE = orbital_energy()->operator()(orbital);

  Eigen::VectorXd occ_evals = orbital_energy()->segment(nfzc, nocc + nfzc);
  Eigen::VectorXd uocc_evals = orbital_energy()->segment(nfzc + nocc, nuocc);

  lcao_factory().keep_partial_transforms(true);

  auto qp_str = L"p";
  // auto qp_str = (orbital_ < 0) ? L"i" : L"a";
  using mpqc::utility::wconcat;
  TArray g_vvog = lcao_factory().compute(wconcat("<a b|G|i ", qp_str, ">[df]"));
  TArray g_oovg = lcao_factory().compute(wconcat("<i j|G|a ", qp_str, ">[df]"));

  if (world.rank() == 0) {
    printf("Iter     SE2(in)     SE2(out)   SE2(delta)\n");
    printf("==== =========== =========== ===========\n");
  }
  size_t iter = 0;
  decltype(SE) SE_diff;
  RowMatrixXd C_dyson;  // Dyson orbital coefficients
  do {
    TArray Sigma_pph;
    {
      TArray dg_vvog = mpqc::dyson::d_pqrE<dyson::Denominator::rEpq>(
          g_vvog, uocc_evals, uocc_evals, occ_evals, SE);
      Sigma_pph("x,y") = 0.5 * (4 * g_vvog("a,b,i,x") - 2 * g_vvog("b,a,i,x")) *
                         dg_vvog("a,b,i,y");
    }
    // std::cout << "Sigma_pph:\n" << Sigma_pph << std::endl;

    TArray Sigma_hhp;
    {
      TArray dg_oovg = mpqc::dyson::d_pqrE<dyson::Denominator::rEpq>(
          g_oovg, occ_evals, occ_evals, uocc_evals, SE);
      Sigma_hhp("x,y") = 0.5 * (4 * g_oovg("i,j,a,x") - 2 * g_oovg("j,i,a,x")) *
                         dg_oovg("i,j,a,y");
    }
    // std::cout << "Sigma_hhp:\n" << Sigma_hhp << std::endl;

    RowMatrixXd Sigma;
    {
      TArray Sigma_ta;
      Sigma_ta("x,y") = Sigma_pph("x,y") + Sigma_hhp("x,y");
      Sigma = array_ops::array_to_eigen(Sigma_ta);
    }

    // std::cout << "Sigma = " << Sigma << std::endl;

    // to keep numerics consistent, diagonalize on 1 node and propagate to the
    // rest
    decltype(SE) SE_updated = 0.0;
    if (world.rank() == 0) {
      RowMatrixXd F_dyson = Sigma + RowMatrixXd(orbital_energy()->asDiagonal());
      Eigen::SelfAdjointEigenSolver<RowMatrixXd> eig_solver(F_dyson);
      auto eps_dyson = eig_solver.eigenvalues();
      SE_updated = eps_dyson(orbital);
      C_dyson = eig_solver.eigenvectors();
      assert(world.size() == 1);
    }

    SE_diff = SE_updated - SE;

    if (world.rank() == 0)
      printf(" %3ld %10.4lf %10.4lf %10.4lf\n", iter, SE, SE_updated, SE_diff);

    SE = SE_updated;
    ++iter;

  } while ((fabs(SE_diff) > 1e-6) && (iter <= max_niter));

  lcao_factory().purge_index(world, qp_str);

  // will use only the target orbital to transform ints
  // create an OrbitalSpace here
  {
    auto& world = this->lcao_factory().get_world();
    auto orbital_registry = this->lcao_factory().orbital_space();
    auto qp_space = orbital_registry->retrieve(OrbitalIndex(qp_str));
    auto C_qp = array_ops::array_to_eigen(qp_space.array());
    auto C_qp_dyson = C_qp * C_dyson;
    auto C_x = C_qp_dyson.block(0, orbital, C_qp_dyson.rows(), 1);
    auto tr_obs = qp_space.array().trange().data().front();
    TA::TiledRange1 tr_x{0, 1};
    auto C_x_ta = array_ops::eigen_to_array<Tile>(world, C_x, tr_obs, tr_x);

    using OrbitalSpaceTArray = OrbitalSpace<TA::DistArray<Tile, Policy>>;
    auto x_space =
        OrbitalSpaceTArray(OrbitalIndex(L"x"), OrbitalIndex(L"κ"), C_x_ta);
    orbital_registry->add(x_space);
  }

  // for now the f12 contribution is energy-independent
  TArray Sigma_pph_f12;
  {
    TArray V_ixix, V_ixxi;
    const bool df = true;  // always do DF
    std::tie(V_ixix, V_ixxi) = mpqc::f12::VX_pqrs_pqsr(
        "V", lcao_factory(), "i", "x", "j", "y", df, use_cabs_);
    Sigma_pph_f12("x,y") =
        0.5 * (1.25 * V_ixix("i,x,j,y") - 0.25 * V_ixxi("i,x,y,j")) *
        lcao_factory()(L"<i|I|j>");
  }
  // std::cout << "Sigma_pph_f12:\n" << Sigma_pph_f12 << std::endl;
  Matrix Sigma_f12 = array_ops::array_to_eigen(Sigma_pph_f12);

  // done with F12 ... remove all geminal ints and CABS indices
  lcao_factory().purge_operator(world, L"R");
  lcao_factory().purge_index(world, L"a'");
  lcao_factory().purge_index(world, L"ρ");

  lcao_factory().keep_partial_transforms(false);

  if (world.rank() == 0) {
    auto SE_F12 = SE + Sigma_f12(0, 0);
    auto Hartree2eV = 27.21138602;
    std::string orblabel =
        std::string(orbital_ < 0 ? "IP" : "EA") + std::to_string(abs(orbital_));
    printf("final       GF2 %6s = %11.3lf eV (%10.4lf a.u.)\n",
           orblabel.c_str(), SE * Hartree2eV, SE);
    printf("final GF2-F12-V %6s = %11.3lf eV (%10.4lf a.u.)\n",
           orblabel.c_str(), SE_F12 * Hartree2eV, SE_F12);
    if (orbital_ > 0)
      printf(
          "WARNING: non-strongly-orthogonal F12 projector is used for the F12 "
          "correction to EA!!!");
  }
}

}  // namespace f12

}  // mpqc

#endif  // MPQC_MP2F12_H
