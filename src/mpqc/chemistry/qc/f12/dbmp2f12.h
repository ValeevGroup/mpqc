//
// Created by Chong Peng on 7/11/16.
//

#ifndef MPQC_DBMP2F12_H
#define MPQC_DBMP2F12_H

#include <mpqc/chemistry/qc/f12/mp2f12.h>
#include <mpqc/chemistry/qc/mbpt/dbmp2.h>
#include <mpqc/chemistry/qc/f12/db_f12_intermediates.h>


namespace mpqc{
namespace f12{


template <typename Tile>
class DBMP2F12 : protected MP2F12<Tile>{
public:

  using Policy = TA::SparsePolicy;
  using TArray = TA::DistArray<Tile, Policy>;
  using MolecularIntegralClass = integrals::MolecularIntegral<Tile,Policy>;

  using real_t = typename MP2F12<Tile>::real_t;
  using Matrix = typename MP2F12<Tile>::Matrix;
  using MP2F12<Tile>::debug;
  using MP2F12<Tile>::mo_integral;
  using MP2F12<Tile>::trange1_engine;
  using MP2F12<Tile>::orbital_energy;

  DBMP2F12() = default;

  DBMP2F12(MolecularIntegralClass& mo_int) : MP2F12<Tile>(construct_db_mp2(mo_int)) {}

  virtual real_t compute(const rapidjson::Document& in){
    auto& world = this->mo_integral().get_world();
    // dbmp2 time
    auto mp2_time0 = mpqc_time::fenced_now(world);

    this->mp2_->compute(in);

    auto mp2_time1 = mpqc_time::fenced_now(world);
    auto mp2_time = mpqc_time::duration_in_s(mp2_time0, mp2_time1);
    mpqc::utility::print_par(world, "Total DBMP2 Time:  ", mp2_time, "\n");

    // start f12
    auto f12_time0 = mpqc_time::fenced_now(world);

    // solve cabs orbitals
    auto orbital_registry = this->mo_integral().orbital_space();
    closed_shell_dualbasis_cabs_mo_build_svd(this->mo_integral(),in,this->mp2_->trange1_engine());

    std::string method = in.HasMember("Method") ? in["Method"].GetString() : "df";

    Matrix mp2_eij, f12_eij;

    if(method == "df"){
      std::tie(mp2_eij, f12_eij) = compute_db_mp2_f12_c_df();
    }
    else{
      throw std::runtime_error("Wrong DBMP2F12 Method");
    }

    if (world.rank() == 0) {
      auto nocc = this->mp2_->trange1_engine()->get_active_occ();
      printf(
          "  i0     i1       eij(mp2)        eij(f12)      eij(mp2-f12) \n"
          "====== ====== =============== =============== ===============\n");
      for (int i = 0; i != nocc; ++i)
        for (int j = i; j != nocc; ++j)
          printf("%4d   %4d   %15.12lf %15.12lf %15.12lf\n", i, j, mp2_eij(i, j),
                 f12_eij(i, j), mp2_eij(i, j) + f12_eij(i, j));
    }

    auto emp2 = mp2_eij.sum();
    auto ef12 = f12_eij.sum();
    if (debug()) {
      utility::print_par(this->mo_integral().get_world(), "E_MP2: ", emp2, "\n");
      utility::print_par(this->mo_integral().get_world(), "E_F12: ", ef12, "\n");
    }

    auto f12_time1 = mpqc_time::fenced_now(world);
    auto f12_time = mpqc_time::duration_in_s(f12_time0, f12_time1);

    mpqc::utility::print_par(world, "Total DBF12 Time:  ", f12_time, "\n");

    return emp2 + ef12;
  }

private:

  std::shared_ptr<mbpt::MP2<Tile,Policy>> construct_db_mp2(MolecularIntegralClass& mo_int){
    std::shared_ptr<mbpt::MP2<Tile,Policy>> dbmp2 = std::make_shared<mbpt::DBMP2<Tile,Policy>>(mo_int);
    return dbmp2;
  }

  /// DBMP2-F12 C approach with density fitting
  std::tuple<Matrix,Matrix> compute_db_mp2_f12_c_df();

};


template <typename Tile>
std::tuple<typename DBMP2F12<Tile>::Matrix, typename DBMP2F12<Tile>::Matrix>
DBMP2F12<Tile>::compute_db_mp2_f12_c_df() {

  auto& world = this->mo_integral().get_world();

  Matrix Eij_MP2, Eij_F12;

  auto nocc = this->mp2_->trange1_engine()->get_active_occ();

  // create shape
  auto occ_tr1 = this->mp2_->trange1_engine()->get_occ_tr1();
  TiledArray::TiledRange occ4_trange({occ_tr1,occ_tr1,occ_tr1,occ_tr1});
  auto ijij_ijji_shape = f12::make_ijij_ijji_shape(occ4_trange);

  // T Term
  TArray t2;
  {
    utility::print_par(world, "Compute T_abij With DF \n" );

    TArray g_abij;
    g_abij("a,b,i,j") = this->mo_integral().compute(L"<i j|G|a b>[df]")("i,j,a,b");
    t2 = mpqc::cc::d_abij(g_abij, *(this->mp2_->orbital_energy()), nocc);

    // compute MP2 energy and pair energies
    TArray TG_ijij_ijji;
    TG_ijij_ijji("i1,j1,i2,j2") =
        (t2("a,b,i1,j1") * g_abij("a,b,i2,j2"))
            .set_shape(ijij_ijji_shape);
    Eij_MP2 = TG_ijij_ijji("i1,j1,i2,j2").reduce(F12PairEnergyReductor<Tile>(2, -1, nocc));
  }

  //compute V term
  TArray V_ijij_ijji = compute_V_ijij_ijji_db_df(mo_integral(), ijij_ijji_shape);
  {
    // G integral in MO not needed, still need G integral in AO to compute F, K, hJ
    this->mo_integral().registry().remove_operation(world, L"G");

    //contribution from V_ijij_ijji
    Matrix eij = V_ijij_ijji("i1,j1,i2,j2").reduce(F12PairEnergyReductor<Tile>(2 * C_ijij_bar,2 * C_ijji_bar,nocc));
    if (debug()) utility::print_par(world, "E_V: ", eij.sum(), "\n");
    Eij_F12 = eij;
  }

  // compute C term
  TArray C_ijab = compute_C_ijab_df(mo_integral());

  {
    utility::print_par(world, "Compute CT With DF \n" );
    V_ijij_ijji("i1,j1,i2,j2") = (C_ijab("i1,j1,a,b")*t2("a,b,i2,j2")).set_shape(ijij_ijji_shape);

    Matrix eij = V_ijij_ijji("i1,j1,i2,j2").reduce(F12PairEnergyReductor<Tile>(2 * C_ijij_bar,2 * C_ijji_bar,nocc));
    if (debug()) utility::print_par(world, "E_CT: ", eij.sum(), "\n");
    Eij_F12 += eij;
  }

  // compute X term
  TArray X_ijij_ijji = compute_X_ijij_ijji_db_df(mo_integral(), ijij_ijji_shape);
  {

    auto Fij = this->mo_integral().compute(L"<i|F|j>[df]");
    auto Fij_eigen = array_ops::array_to_eigen(Fij);
    f12::convert_X_ijkl(X_ijij_ijji, Fij_eigen);

    Matrix eij = X_ijij_ijji("i1,j1,i2,j2").reduce(F12PairEnergyReductor<Tile>(CC_ijij_bar,CC_ijji_bar,nocc));
    eij *= -1;
    if (debug()) utility::print_par(world, "E_X: ", eij.sum(), "\n");
    Eij_F12 += eij;

  }

  // compute B term
  TArray B_ijij_ijji = compute_B_ijij_ijji_db_df(mo_integral(), ijij_ijji_shape);
  {
    Matrix eij = B_ijij_ijji("i1,j1,i2,j2").reduce(F12PairEnergyReductor<Tile>(CC_ijij_bar,CC_ijji_bar,nocc));
    if (debug()) utility::print_par(world, "E_B: ", eij.sum(), "\n");
    Eij_F12 += eij;
  }

  {
    utility::print_par(world, "Compute CC Term With DF \n");
    auto C_bar_ijab = f12::convert_C_ijab(C_ijab, nocc, *(this->mp2_->orbital_energy()));
    B_ijij_ijji("i1,j1,i2,j2") = (C_ijab("i1,j1,a,b")*C_bar_ijab("i2,j2,a,b")).set_shape(ijij_ijji_shape);

    Matrix eij = B_ijij_ijji("i1,j1,i2,j2").reduce(F12PairEnergyReductor<Tile>(CC_ijij_bar,CC_ijji_bar,nocc));
    if (debug()) utility::print_par(world, "E_CC: ", eij.sum(), "\n");
    Eij_F12 += eij;
  }

  return std::make_tuple(Eij_MP2,Eij_F12);
}

}  // namespace f12

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
class DBGF2F12 {
public:

  using Policy = TA::SparsePolicy;
  using TArray = TA::DistArray<Tile, Policy>;
  using MolecularIntegralClass = integrals::MolecularIntegral<Tile,Policy>;

  using real_t = typename DBMP2F12<Tile>::real_t;
  using Matrix = typename DBMP2F12<Tile>::Matrix;

  DBGF2F12() = default;

  DBGF2F12(MolecularIntegralClass& mo_int) : mp2f12_(std::make_shared<DBMP2F12<Tile>>(mo_int)) {}

  MolecularIntegralClass &mo_integral() const { return mp2f12_->mo_integral(); }

  const std::shared_ptr<TRange1Engine> trange1_engine() const {
    return mp2f12_->trange1_engine_;
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
      throw std::runtime_error("DBGF2F12::Orbital must be positive (for particles) or negative (for holes)");

    std::string method = in.HasMember("DysonMethod") ? in["DysonMethod"].GetString() : "diagonal-fixed";
    TA_USER_ASSERT(method == "diagonal-fixed" || method == "diagonal-iterative",
                   "DBGF2F12: unknown value for keyword \"method\"");

    std::cout << "orbital = " << orbital_ << " method = " << method << std::endl;

    if(method == "diagonal-fixed") {
      compute_diagonal_fixed();
    }
    else if (method == "diagonal-iterative") {
      compute_diagonal_iterative();
    }

    auto time1 = mpqc_time::fenced_now(world);
    auto time = mpqc_time::duration_in_s(time0, time1);

    mpqc::utility::print_par(world, "Total DBGF2F12 Time:  ", time, "\n");

    return 0.0;
  }

private:

  std::shared_ptr<DBMP2F12<Tile>> mp2f12_;
  int orbital_;

  void compute_diagonal_fixed();
  void compute_diagonal_iterative();
};

template <typename Tile>
void DBGF2F12<Tile>::compute_diagonal_iterative() {

  auto nocc = this->mp2f12_->trange1_engine()->get_active_occ();
  auto nuocc = this->mp2f12_->trange1_engine()->get_vir();
  const auto abs_orbital = nocc + orbital_;
  auto SE = orbital_energy()->operator()(abs_orbital);

  double ediff = 0.0, elast = 0.0;
  size_t iter = 0;

  Eigen::VectorXd occ_evals = orbital_energy()->segment(0,nocc);
  Eigen::VectorXd uocc_evals = orbital_energy()->segment(nocc, nuocc);

  do {
    iter++;

    TArray Sigma_pph;
    {
      TArray g_vvog = mo_integral().compute(L"<a b|G|i j>[df]");
      TArray dg_vvog = mpqc::dyson::d_pqrE<dyson::Denominator::rEpq>(
          g_vvog, uocc_evals, uocc_evals, occ_evals, SE);
      Sigma_pph("p,q") = 0.5 * (4*g_vvog("a,b,i,p") - 2*g_vvog("b,a,i,p")) * dg_vvog("a,b,i,q");
    }
    std::cout << "Sigma_pph:\n" << Sigma_pph << std::endl;

    TArray Sigma_hhp;
    {
      TArray g_oovg = mo_integral().compute(L"<i j|G|a j>[df]");
      TArray dg_oovg = mpqc::dyson::d_pqrE<dyson::Denominator::rEpq>(
          g_oovg, occ_evals, occ_evals, uocc_evals, SE);
      Sigma_hhp("p,q") = 0.5 * (4*g_oovg("i,j,a,p") - 2*g_oovg("j,i,a,p")) * dg_oovg("i,j,a,q");
    }
    std::cout << "Sigma_hhp:\n" << Sigma_hhp << std::endl;

    RowMatrixXd Sigma;
    {
      TArray Sigma_ta;
      Sigma_ta("p,q") = Sigma_pph("p,q") + Sigma_hhp("p,q");
      Sigma = array_ops::array_to_eigen(Sigma_ta);
    }

    std::cout << "orbital_energy = " << *orbital_energy() << std::endl;
    std::cout << "sigma2 = " << Sigma << std::endl;

    SE = Sigma(abs_orbital, abs_orbital) + orbital_energy()->operator()(abs_orbital);

    ediff = elast - SE;
    elast = SE;

  } while ((fabs(ediff) > 1e-12) && (iter < 100));
}

template <typename Tile>
void DBGF2F12<Tile>::compute_diagonal_fixed() {
  assert(false && "not yet implemented");
}

} // end of namespace f12
} // end of namespace mpqc



#endif //MPQC_DBMP2F12_H
