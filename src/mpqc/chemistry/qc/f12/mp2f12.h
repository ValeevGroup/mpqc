//
// Created by Chong Peng on 3/31/16.
//

#ifndef MPQC_MP2F12_H
#define MPQC_MP2F12_H

#include <string>

#include "../../../../../include/eigen.h"
#include "../../../../../utility/cc_utility.h"
#include "../../../../../utility/trange1_engine.h"
#include <mpqc/chemistry/qc/f12/f12_intermediates.h>
#include <mpqc/chemistry/qc/f12/f12_utility.h>
#include <mpqc/chemistry/qc/mbpt/mp2.h>

namespace mpqc {
namespace f12 {

template <typename Tile>
class MP2F12 {
 public:
  using Policy = TA::SparsePolicy;
  using TArray = TA::DistArray<Tile, Policy>;
  using MolecularIntegralClass = integrals::MolecularIntegral<Tile, Policy>;

  using real_t = typename Tile::scalar_type;
  using Matrix = RowMatrix<real_t>;

  MP2F12() = default;

  MP2F12(MolecularIntegralClass& mo_int) {
    mp2_ = std::make_shared<mbpt::MP2<Tile, Policy>>(mo_int);
  }

  MP2F12(std::shared_ptr<mbpt::MP2<Tile, Policy>> mp2)
      : mp2_(mp2) {}

  MolecularIntegralClass &mo_integral() const { return mp2_->mo_integral(); }

  const std::shared_ptr<TRange1Engine> trange1_engine() const {
    return mp2_->trange1_engine();
  }

  const std::shared_ptr<Eigen::VectorXd> orbital_energy() const {
    return mp2_->orbital_energy();
  }

  std::tuple<Matrix, Matrix> compute_mp2_f12_c_df();

  std::tuple<Matrix, Matrix> compute_mp2_f12_c();

  real_t compute(const rapidjson::Document& in) {
    auto& world = mo_integral().get_world();
    auto f12_time0 = mpqc_time::fenced_now(world);
    mp2_->init(in);

    // solve cabs orbitals
    auto orbital_registry = mo_integral().orbital_space();
    closed_shell_cabs_mo_build_svd(mo_integral(), in,
                                   this->mp2_->trange1_engine());

    std::string method =
        in.HasMember("Method") ? in["Method"].GetString() : "df";

    Matrix mp2_eij, f12_eij;

    if (method == "four center") {
      std::tie(mp2_eij, f12_eij) = compute_mp2_f12_c();
    } else if (method == "df") {
      std::tie(mp2_eij, f12_eij) = compute_mp2_f12_c_df();
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
      utility::print_par(mo_integral().get_world(), "E_MP2: ", emp2, "\n");
      utility::print_par(mo_integral().get_world(), "E_F12: ", ef12, "\n");
    }

    auto f12_time1 = mpqc_time::fenced_now(world);
    auto f12_time = mpqc_time::duration_in_s(f12_time0, f12_time1);
    mpqc::utility::print_par(world, "Total MP2F12 Time:  ", f12_time, "\n");

    return emp2 + ef12;
  }

 protected:
  std::shared_ptr<mbpt::MP2<Tile, Policy>> mp2_;

  int debug() const { return 1; }
};

template <typename Tile>
std::tuple<typename MP2F12<Tile>::Matrix, typename MP2F12<Tile>::Matrix>
MP2F12<Tile>::compute_mp2_f12_c_df() {
  auto& world = mo_integral().get_world();

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
    g_abij("a,b,i,j") = mo_integral().compute(L"<i j|G|a b>[df]")("i,j,a,b");
    t2 = mpqc::cc::d_abij(g_abij, *(mp2_->orbital_energy()), n_occ, n_frozen);

    // compute MP2 energy and pair energies
    TArray TG_ijij_ijji;
    TG_ijij_ijji("i1,j1,i2,j2") =
        (t2("a,b,i1,j1") * g_abij("a,b,i2,j2")).set_shape(ijij_ijji_shape);
    Eij_MP2 = TG_ijij_ijji("i1,j1,i2,j2")
                  .reduce(F12PairEnergyReductor<Tile>(2, -1, n_active_occ));
  }

  // compute V term
  TArray V_ijij_ijji = compute_V_ijij_ijji_df(mo_integral(), ijij_ijji_shape);
  {
    // G integral in MO not needed, still need G integral in AO to compute F, K,
    // hJ
    mo_integral().registry().remove_operation(world, L"G");

    // contribution from V_ijij_ijji
    // NB factor of 2 from the Hylleraas functional
    Eij_F12 = V_ijij_ijji("i1,j1,i2,j2")
                  .reduce(F12PairEnergyReductor<Tile>(2 * C_ijij_bar,
                                                      2 * C_ijji_bar, n_active_occ));
    if (debug()) utility::print_par(world, "E_V: ", Eij_F12.sum(), "\n");
  }

  // compute C term
  TArray C_ijab = compute_C_ijab_df(mo_integral());

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
  TArray X_ijij_ijji = compute_X_ijij_ijji_df(mo_integral(), ijij_ijji_shape);
  {
    // R_ipjq not needed
    mo_integral().registry().remove_formula(world, L"<i1 j1|R|p q>[df]");

    auto Fij = mo_integral().compute(L"<i|F|j>[df]");
    auto Fij_eigen = array_ops::array_to_eigen(Fij);
    f12::convert_X_ijkl(X_ijij_ijji, Fij_eigen);

    Matrix Eij_x = X_ijij_ijji("i1,j1,i2,j2")
                       .reduce(F12PairEnergyReductor<Tile>(CC_ijij_bar,
                                                           CC_ijji_bar, n_active_occ));
    Eij_x *= -1.0;
    if (debug()) utility::print_par(world, "E_X: ", Eij_x.sum(), "\n");
    Eij_F12 += Eij_x;
  }

  // compute B term
  TArray B_ijij_ijji = compute_B_ijij_ijji_df(mo_integral(), ijij_ijji_shape);
  {
    Matrix Eij_b = B_ijij_ijji("i1,j1,i2,j2")
                       .reduce(F12PairEnergyReductor<Tile>(CC_ijij_bar,
                                                           CC_ijji_bar, n_active_occ));
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
                        .reduce(F12PairEnergyReductor<Tile>(CC_ijij_bar,
                                                            CC_ijji_bar, n_active_occ));
    if (debug()) utility::print_par(world, "E_CC: ", Eij_cc.sum(), "\n");
    Eij_F12 += Eij_cc;
  }

  return std::make_tuple(Eij_MP2, Eij_F12);
}

template <typename Tile>
std::tuple<typename MP2F12<Tile>::Matrix, typename MP2F12<Tile>::Matrix>
MP2F12<Tile>::compute_mp2_f12_c() {
  auto& world = mo_integral().get_world();

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
    g_abij("a,b,i,j") = mo_integral().compute(L"<i j|G|a b>")("i,j,a,b");
    t2_nodf = mpqc::cc::d_abij(g_abij, *(mp2_->orbital_energy()), n_occ, n_frozen);
    TArray TG_ijij_ijji_nodf;
    TG_ijij_ijji_nodf("i1,j1,i2,j2") =
        (t2_nodf("a,b,i1,j1") * g_abij("a,b,i2,j2")).set_shape(ijij_ijji_shape);
    Eij_MP2 = TG_ijij_ijji_nodf("i1,j1,i2,j2")
                  .reduce(F12PairEnergyReductor<Tile>(2, -1, n_active_occ));
  }

  TArray V_ijij_ijji_nodf = compute_V_ijij_ijji(mo_integral(), ijij_ijji_shape);
  {
    Eij_F12 = V_ijij_ijji_nodf("i1,j1,i2,j2")
                  .reduce(F12PairEnergyReductor<Tile>(2 * C_ijij_bar,
                                                      2 * C_ijji_bar, n_active_occ));
    if (debug()) utility::print_par(world, "E_V: ", Eij_F12.sum(), "\n");
  }

  TArray C_ijab_nodf = compute_C_ijab(mo_integral());

  {
    utility::print_par(world, "Compute CT Without DF \n");
    V_ijij_ijji_nodf("i1,j1,i2,j2") =
        (C_ijab_nodf("i1,j1,a,b") * t2_nodf("a,b,i2,j2"))
            .set_shape(ijij_ijji_shape);

    Matrix Eij_ct = V_ijij_ijji_nodf("i1,j1,i2,j2")
                        .reduce(F12PairEnergyReductor<Tile>(
                            2 * C_ijij_bar, 2 * C_ijji_bar, n_active_occ));
    if (debug()) utility::print_par(world, "E_CT: ", Eij_ct, "\n");
    Eij_F12 += Eij_ct;
  }

  TArray X_ijij_ijji_nodf = compute_X_ijij_ijji(mo_integral(), ijij_ijji_shape);
  {
    // compute energy contribution
    auto Fij = mo_integral().compute(L"<i|F|j>");
    auto Fij_eigen = array_ops::array_to_eigen(Fij);
    f12::convert_X_ijkl(X_ijij_ijji_nodf, Fij_eigen);

    Matrix Eij_x =
        X_ijij_ijji_nodf("i1,j1,i2,j2")
            .reduce(F12PairEnergyReductor<Tile>(CC_ijij_bar, CC_ijji_bar, n_active_occ));
    Eij_x *= -1;
    if (debug()) utility::print_par(world, "E_X: ", Eij_x.sum(), "\n");
    Eij_F12 += Eij_x;
  }

  TArray B_ijij_ijji_nodf = compute_B_ijij_ijji(mo_integral(), ijij_ijji_shape);
  {
    Matrix Eij_b =
        B_ijij_ijji_nodf("i1,j1,i2,j2")
            .reduce(F12PairEnergyReductor<Tile>(CC_ijij_bar, CC_ijji_bar, n_active_occ));
    if (debug()) utility::print_par(world, "E_B: ", Eij_b.sum(), "\n");
    Eij_F12 += Eij_b;
  }

  {
    utility::print_par(world, "Compute CC Term Without DF \n");
    auto C_bar_ijab =
        f12::convert_C_ijab(C_ijab_nodf, n_occ, n_frozen, *(mp2_->orbital_energy()));
    B_ijij_ijji_nodf("i1,j1,i2,j2") =
        (C_ijab_nodf("i1,j1,a,b") * C_bar_ijab("i2,j2,a,b"))
            .set_shape(ijij_ijji_shape);

    Matrix Eij_cc =
        B_ijij_ijji_nodf("i1,j1,i2,j2")
            .reduce(F12PairEnergyReductor<Tile>(CC_ijij_bar, CC_ijji_bar, n_active_occ));
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
    TA::Array<double, 4, Tile, Policy>& pqrs,
    const Eigen::VectorXd& evals_p,
    const Eigen::VectorXd& evals_q,
    const Eigen::VectorXd& evals_r,
    typename Tile::scalar_type E) {
  auto convert = [evals_p, evals_q, evals_r, E](Tile& result_tile, const Tile& arg_tile) {

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
  using MolecularIntegralClass = integrals::MolecularIntegral<Tile,Policy>;

  using real_t = typename Tile::scalar_type;
  using Matrix = RowMatrix<real_t>;

  GF2F12() = default;

  GF2F12(MolecularIntegralClass& mo_int) : mp2f12_(std::make_shared<mbpt::MP2<Tile, Policy>>(mo_int)) {}

  MolecularIntegralClass &mo_integral() const { return mp2f12_->mo_integral(); }

  const std::shared_ptr<TRange1Engine> trange1_engine() const {
    return mp2f12_->trange1_engine();
  }

  const std::shared_ptr<Eigen::VectorXd> orbital_energy() const {
    return mp2f12_->orbital_energy();
  }

  virtual real_t compute(const rapidjson::Document& in){
    auto& world = this->mo_integral().get_world();

    this->mp2f12_->compute(in);

    auto time0 = mpqc_time::fenced_now(world);

    orbital_ = in.HasMember("Orbital") ? in["Orbital"].GetInt() : -1;
    if (orbital_ == 0)
      throw std::runtime_error("GF2F12::Orbital must be positive (for particles) or negative (for holes)");
    if (orbital_ < 0 && (abs(orbital_)-1 >= trange1_engine()->get_active_occ()))
      throw std::runtime_error("GF2F12::orbital is invalid (the number of holes exceeded)");
    if (orbital_ > 0 && (orbital_-1 >= trange1_engine()->get_vir()))
      throw std::runtime_error("GF2F12::orbital is invalid (the number of particles exceeded)");

    std::string method = in.HasMember("DysonMethod") ? in["DysonMethod"].GetString() : "diagonal-fixed";
    TA_USER_ASSERT(method == "diagonal-fixed" || method == "diagonal-iterative",
                   "GF2F12: unknown value for keyword \"method\"");

    std::cout << "orbital = " << orbital_ << " method = " << method << std::endl;

    compute_diagonal(method == "diagonal-fixed" ? 0 : 100);

    auto time1 = mpqc_time::fenced_now(world);
    auto time = mpqc_time::duration_in_s(time0, time1);

    mpqc::utility::print_par(world, "Total GF2F12 Time:  ", time, "\n");

    return 0.0;
  }

private:

  std::shared_ptr<mbpt::MP2<Tile, Policy>> mp2f12_;
  int orbital_;

  /// use self-energy in diagonal representation
  void compute_diagonal(int max_niter = 0);
};

template <typename Tile>
void GF2F12<Tile>::compute_diagonal(int max_niter) {

  auto nfzc = this->mp2f12_->trange1_engine()->get_nfrozen();
  auto nocc = this->mp2f12_->trange1_engine()->get_active_occ();
  auto nuocc = this->mp2f12_->trange1_engine()->get_vir();
  // map orbital_ (index relative to Fermi level) to orbital index in active orbital space
  const auto act_orbital = (orbital_ < 0) ? nocc + orbital_ : nocc + orbital_ - 1;
  auto SE = orbital_energy()->operator()(act_orbital);

  Eigen::VectorXd occ_evals = orbital_energy()->segment(0,nocc);
  Eigen::VectorXd uocc_evals = orbital_energy()->segment(nocc, nuocc);

  // will use only the target orbital to transform ints
  // create an OrbitalSpace here
  {
    auto& world = this->mo_integral().get_world();
    auto orbital_registry = this->mo_integral().orbital_space();
    auto p_space = orbital_registry->retrieve(OrbitalIndex(L"p"));
    auto C_p = array_ops::array_to_eigen(p_space.array());
    auto orbital = nfzc + act_orbital;
    auto C_x = C_p.block(0, orbital, C_p.rows(), 1);
    auto tr_obs = p_space.array().trange().data().front();
    TA::TiledRange1 tr_x{0,1};
    auto C_x_ta = array_ops::eigen_to_array<Tile>(world,
                                                  C_x,
                                                  tr_obs,
                                                  tr_x);

    using OrbitalSpaceTArray = OrbitalSpace<TA::DistArray<Tile,Policy>>;
    auto x_space = OrbitalSpaceTArray(OrbitalIndex(L"x"),OrbitalIndex(L"κ"), C_x_ta);
    orbital_registry->add(x_space);
  }

  printf("Iter     SE2(in)     SE2(out)   SE2(delta)\n");
  printf("==== =========== =========== ===========\n");

  size_t iter = 0;
  decltype(SE) SE_diff;
  do {

    TArray Sigma_pph;
    {
      TArray& g_vvog = mo_integral().compute(L"<a b|G|i x>[df]");
      TArray dg_vvog = mpqc::dyson::d_pqrE<dyson::Denominator::rEpq>(
          g_vvog, uocc_evals, uocc_evals, occ_evals, SE);
      Sigma_pph("x,y") = 0.5 * (4*g_vvog("a,b,i,x") - 2*g_vvog("b,a,i,x")) * dg_vvog("a,b,i,y");
    }
    //std::cout << "Sigma_pph:\n" << Sigma_pph << std::endl;

    TArray Sigma_hhp;
    {
      TArray& g_oovg = mo_integral().compute(L"<i j|G|a x>[df]");
      TArray dg_oovg = mpqc::dyson::d_pqrE<dyson::Denominator::rEpq>(
          g_oovg, occ_evals, occ_evals, uocc_evals, SE);
      Sigma_hhp("x,y") = 0.5 * (4*g_oovg("i,j,a,x") - 2*g_oovg("j,i,a,x")) * dg_oovg("i,j,a,y");
    }
    //std::cout << "Sigma_hhp:\n" << Sigma_hhp << std::endl;

    RowMatrixXd Sigma;
    {
      TArray Sigma_ta;
      Sigma_ta("x,y") = Sigma_pph("x,y") + Sigma_hhp("x,y");
      Sigma = array_ops::array_to_eigen(Sigma_ta);
    }

    //std::cout << "Sigma = " << Sigma << std::endl;

    auto SE_updated = Sigma(0, 0) + orbital_energy()->operator()(act_orbital);
    SE_diff = SE_updated - SE;

    printf(" %3ld %11.3lf %11.3lf %11.3lf", iter, SE, SE_updated, SE_diff);

    SE = SE_updated;
    ++iter;

  } while ((fabs(SE_diff) > 1e-6) && (iter <= max_niter));
}

}  // namespace f12

}  // mpqc

#endif  // MPQC_MP2F12_H
