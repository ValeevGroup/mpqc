//
// Created by Chong Peng on 10/11/16.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_F12_GF2F12_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_F12_GF2F12_H_

#include "mpqc/chemistry/qc/f12/f12_intermediates.h"
#include "mpqc/chemistry/qc/scf/mo_build.h"
#include "mpqc/chemistry/qc/wfn/lcao_wfn.h"
#include "mpqc/mpqc_config.h"
#include <mpqc/chemistry/units/units.h>

namespace mpqc {
namespace lcao {

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

/**
 *  \brief GF2F12 class
 *  keyval name for this class GF2F12
 *
 */

template <typename Tile>
class GF2F12 : public LCAOWavefunction<Tile, TA::SparsePolicy> {
 public:
  using Policy = TA::SparsePolicy;
  using TArray = TA::DistArray<Tile, Policy>;
  using LCAOFactoryType = LCAOFactory<Tile, Policy>;

  using real_t = typename Tile::scalar_type;
  using Matrix = RowMatrix<real_t>;

  GF2F12() = default;
  virtual ~GF2F12() = default;

  // clang-format off
  /**
   * KeyVal constructor
   * @param kv
   *
   * keywords: all keywords for LCAOWavefuction
   *
   * | KeyWord | Type | Default| Description |
   * |---------|------|--------|-------------|
   * | ref | Wavefunction | none | reference Wavefunction, RHF for example |
   * | orbital | int | -1 | target orbital, by default is HOMO |
   * | use_cabs | bool | true | if includes cabs in F12-V term |
   * | dyson_method | string | diagonal | dyson_method to use, (diagonal or nondiagonal ) |
   * | max_iter | int | 100 | maximum iteration |
   */
  // clang-format on

  GF2F12(const KeyVal& kv) : LCAOWavefunction<Tile, Policy>(kv) {
    if (kv.exists("ref")) {
      ref_wfn_ = kv.keyval("ref").class_ptr<Wavefunction>();
    } else {
      throw std::invalid_argument(
          "Default Ref Wfn in GF2F12 is not support! \n");
    }
    orbital_ = kv.value<int>("orbital", -1);
    use_cabs_ = kv.value<bool>("use_cabs", true);
    dyson_method_ = kv.value<std::string>("dyson_method", "diagonal");
    max_iter_ = kv.value<int>("max_iter", 100);
  }

  void compute(PropertyBase* pb) override {
    throw std::runtime_error("Not Implemented!! \n");
  }

  real_t value() override {
    using mpqc::utility::print_par;

    auto& world = this->lcao_factory().world();

    this->energy_ = ref_wfn_->value();

    // init
    init();

    auto time0 = mpqc::fenced_now(world);

    if (orbital_ == 0)
      throw std::runtime_error(
          "GF2F12::Orbital must be positive (for particles) or negative (for "
          "holes)");
    if (orbital_ < 0 &&
        (abs(orbital_) - 1 >= this->trange1_engine()->get_active_occ()))
      throw std::runtime_error(
          "GF2F12::orbital is invalid (the number of holes exceeded)");
    if (orbital_ > 0 && (orbital_ - 1 >= this->trange1_engine()->get_vir()))
      throw std::runtime_error(
          "GF2F12::orbital is invalid (the number of particles exceeded)");

    std::string method_str = dyson_method_;
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
      compute_diagonal(max_iter_);
    else
      compute_nondiagonal(max_iter_);

    auto time1 = mpqc::fenced_now(world);
    auto time = mpqc::duration_in_s(time0, time1);

    print_par(world, "Total GF2F12 Time:  ", time, "\n");

    return this->energy_;
  }

  bool use_cabs() const { return use_cabs_; }
  int orbital() const { return orbital_; }

 private:
  /// initialize obs and cabs orbitals
  virtual void init();

  /// initialize target orbital in compute_diagonal function
  virtual void init_target_orbital_diagonal();

  /// compute V_ixjy and V_ixyj term in compute_diagonal and compute_nondiagonal
  virtual std::tuple<TArray, TArray> compute_V(){
    return f12::VX_pqrs_pqsr("V", this->lcao_factory(), "i", "x", "j", "y", true, use_cabs_);
  }

  /// use self-energy in diagonal representation
  void compute_diagonal(int max_niter = 100);

  /// use non-diagonal self-energy
  void compute_nondiagonal(int max_niter = 100);

 private:
  std::shared_ptr<Wavefunction> ref_wfn_;
  int orbital_;
  bool use_cabs_;
  std::string dyson_method_;
  std::size_t max_iter_;
};

template <typename Tile>
void GF2F12<Tile>::init() {
  // init obs
  auto mol = this->lcao_factory().ao_factory().molecule();
  Eigen::VectorXd orbital_energy;
  this->trange1_engine_ = closed_shell_obs_mo_build_eigen_solve(
      this->lcao_factory(), orbital_energy, mol, this->is_frozen_core(),
      this->occ_block(), this->unocc_block());
  this->orbital_energy_ = std::make_shared<Eigen::VectorXd>(orbital_energy);

  if (use_cabs_) {
    // compute cabs
    closed_shell_cabs_mo_build_svd(this->lcao_factory(), this->trange1_engine(),
                                   this->unocc_block());
  }
}

template <typename Tile>
void GF2F12<Tile>::init_target_orbital_diagonal() {
  auto nfzc = this->trange1_engine()->get_nfrozen();
  auto nocc = this->trange1_engine()->get_active_occ();
  const auto orbital = nfzc + nocc + ((orbital_ < 0) ? orbital_ : orbital_ - 1);

  auto& world = this->wfn_world()->world();
  auto& orbital_registry = this->lcao_factory().orbital_space();
  auto p_space = orbital_registry.retrieve(OrbitalIndex(L"p"));
  auto C_p = array_ops::array_to_eigen(p_space.coefs());
  auto C_x = C_p.block(0, orbital, C_p.rows(), 1);
  auto tr_obs = p_space.coefs().trange().data().front();
  TA::TiledRange1 tr_x{0, 1};
  auto C_x_ta = array_ops::eigen_to_array<Tile, TA::SparsePolicy>(world, C_x,
                                                                  tr_obs, tr_x);

  using OrbitalSpaceTArray = OrbitalSpace<TA::DistArray<Tile, Policy>>;
  auto x_space =
      OrbitalSpaceTArray(OrbitalIndex(L"x"), OrbitalIndex(L"κ"), C_x_ta);
}

template <typename Tile>
void GF2F12<Tile>::compute_diagonal(int max_niter) {
  auto& world = this->wfn_world()->world();

  auto nfzc = this->trange1_engine()->get_nfrozen();
  auto nocc = this->trange1_engine()->get_active_occ();
  auto nuocc = this->trange1_engine()->get_vir();
  // map orbital_ (index relative to Fermi level) to orbital index in active
  // orbital space
  const auto orbital = nfzc + nocc + ((orbital_ < 0) ? orbital_ : orbital_ - 1);
  auto SE = this->orbital_energy()->operator()(orbital);

  Eigen::VectorXd occ_evals =
      this->orbital_energy()->segment(nfzc, nocc + nfzc);
  Eigen::VectorXd uocc_evals =
      this->orbital_energy()->segment(nfzc + nocc, nuocc);

  // will use only the target orbital to transform ints
  // create an OrbitalSpace here
  { init_target_orbital_diagonal(); }

  this->lcao_factory().keep_partial_transforms(true);

  TArray g_vvog = this->lcao_factory().compute(L"<a b|G|i x>[df]");
  TArray g_oovg = this->lcao_factory().compute(L"<i j|G|a x>[df]");

  if (world.rank() == 0) {
    std::printf("Iter     SE2(in)     SE2(out)   SE2(delta)\n");
    std::printf("==== =========== =========== ===========\n");
  }
  size_t iter = 0;
  decltype(SE) SE_diff;
  do {
    TArray Sigma_pph;
    {
      TArray dg_vvog = dyson::d_pqrE<dyson::Denominator::rEpq>(
          g_vvog, uocc_evals, uocc_evals, occ_evals, SE);
      Sigma_pph("x,y") = 0.5 * (4 * g_vvog("a,b,i,x") - 2 * g_vvog("b,a,i,x")) *
                         dg_vvog("a,b,i,y");
    }
    // std::cout << "Sigma_pph:\n" << Sigma_pph << std::endl;

    TArray Sigma_hhp;
    {
      TArray dg_oovg = dyson::d_pqrE<dyson::Denominator::rEpq>(
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

    auto SE_updated = Sigma(0, 0) + this->orbital_energy()->operator()(orbital);
    SE_diff = SE_updated - SE;

    if (world.rank() == 0)
      std::printf(" %3ld %10.4lf %10.4lf %10.4lf\n", iter, SE, SE_updated,
                  SE_diff);

    SE = SE_updated;
    ++iter;

  } while ((fabs(SE_diff) > 1e-6) && (iter <= max_niter));

  this->lcao_factory().purge_formula(world, L"<a b|G|i x>[df]");
  this->lcao_factory().purge_formula(world, L"<i j|G|a x>[df]");

  // for now the f12 contribution is energy-independent
  TArray Sigma_pph_f12;
  {
    TArray V_ixix, V_ixxi;
    // always do DF
    std::tie(V_ixix, V_ixxi) = compute_V();

    Sigma_pph_f12("x,y") =
        0.5 * (1.25 * V_ixix("i,x,j,y") - 0.25 * V_ixxi("i,x,y,j")) *
        this->lcao_factory()(L"<i|I|j>");
  }
  // std::cout << "Sigma_pph_f12:\n" << Sigma_pph_f12 << std::endl;
  Matrix Sigma_f12 = array_ops::array_to_eigen(Sigma_pph_f12);

  // done with F12 ... remove all geminal ints and CABS indices
  this->lcao_factory().purge_operator(world, L"R");
  this->lcao_factory().purge_index(world, L"a'");
  this->lcao_factory().purge_index(world, L"ρ");

  this->lcao_factory().keep_partial_transforms(false);

  if (world.rank() == 0) {
    auto SE_F12 = SE + Sigma_f12(0, 0);
    auto Hartree2eV = 27.21138602;
    std::string orblabel =
        std::string(orbital_ < 0 ? "IP" : "EA") + std::to_string(abs(orbital_));
    std::printf("final       GF2 %6s = %11.3lf eV (%10.4lf a.u.)\n",
                orblabel.c_str(), SE * Hartree2eV, SE);
    std::printf("final GF2-F12-V %6s = %11.3lf eV (%10.4lf a.u.)\n",
                orblabel.c_str(), SE_F12 * Hartree2eV, SE_F12);
    if (orbital_ > 0)
      std::printf(
          "WARNING: non-strongly-orthogonal F12 projector is used for the F12 "
          "correction to EA!!!");
  }
}

template <typename Tile>
void GF2F12<Tile>::compute_nondiagonal(int max_niter) {
  auto& world = this->lcao_factory().world();

  auto nfzc = this->trange1_engine()->get_nfrozen();
  auto nocc = this->trange1_engine()->get_active_occ();
  auto nuocc = this->trange1_engine()->get_vir();
  // map orbital_ (index relative to Fermi level) to orbital index in active
  // orbital space
  const auto orbital = nfzc + nocc + ((orbital_ < 0) ? orbital_ : orbital_ - 1);
  auto SE = this->orbital_energy()->operator()(orbital);

  Eigen::VectorXd occ_evals =
      this->orbital_energy()->segment(nfzc, nocc + nfzc);
  Eigen::VectorXd uocc_evals =
      this->orbital_energy()->segment(nfzc + nocc, nuocc);

  this->lcao_factory().keep_partial_transforms(true);

  auto qp_str = L"p";
  // auto qp_str = (orbital_ < 0) ? L"i" : L"a";
  using mpqc::utility::wconcat;
  TArray g_vvog =
      this->lcao_factory().compute(wconcat("<a b|G|i ", qp_str, ">[df]"));
  TArray g_oovg =
      this->lcao_factory().compute(wconcat("<i j|G|a ", qp_str, ">[df]"));

  if (world.rank() == 0) {
    std::printf("Iter     SE2(in)     SE2(out)   SE2(delta)\n");
    std::printf("==== =========== =========== ===========\n");
  }
  size_t iter = 0;
  decltype(SE) SE_diff;
  RowMatrixXd C_dyson;  // Dyson orbital coefficients
  do {
    TArray Sigma_pph;
    {
      TArray dg_vvog = dyson::d_pqrE<dyson::Denominator::rEpq>(
          g_vvog, uocc_evals, uocc_evals, occ_evals, SE);
      Sigma_pph("x,y") = 0.5 * (4 * g_vvog("a,b,i,x") - 2 * g_vvog("b,a,i,x")) *
                         dg_vvog("a,b,i,y");
    }
    // std::cout << "Sigma_pph:\n" << Sigma_pph << std::endl;

    TArray Sigma_hhp;
    {
      TArray dg_oovg = dyson::d_pqrE<dyson::Denominator::rEpq>(
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
      RowMatrixXd F_dyson =
          Sigma + RowMatrixXd(this->orbital_energy()->asDiagonal());
      Eigen::SelfAdjointEigenSolver<RowMatrixXd> eig_solver(F_dyson);
      auto eps_dyson = eig_solver.eigenvalues();
      SE_updated = eps_dyson(orbital);
      C_dyson = eig_solver.eigenvectors();
      assert(world.size() == 1);
    }

    SE_diff = SE_updated - SE;

    if (world.rank() == 0)
      std::printf(" %3ld %10.4lf %10.4lf %10.4lf\n", iter, SE, SE_updated,
                  SE_diff);

    SE = SE_updated;
    ++iter;

  } while ((fabs(SE_diff) > 1e-6) && (iter <= max_niter));

  this->lcao_factory().purge_index(world, qp_str);

  // will use only the target orbital(x) to transform ints
  // create an OrbitalSpace here
  {
    auto& world = this->lcao_factory().world();
    auto& orbital_registry = this->lcao_factory().orbital_space();
    auto qp_space = orbital_registry.retrieve(OrbitalIndex(qp_str));
    auto C_qp = array_ops::array_to_eigen(qp_space.coefs());
    auto C_qp_dyson = C_qp * C_dyson;
    auto C_x = C_qp_dyson.block(0, orbital, C_qp_dyson.rows(), 1);
    auto tr_obs = qp_space.coefs().trange().data().front();
    TA::TiledRange1 tr_x{0, 1};
    auto C_x_ta = array_ops::eigen_to_array<Tile, TA::SparsePolicy>(
        world, C_x, tr_obs, tr_x);

    using OrbitalSpaceTArray = OrbitalSpace<TA::DistArray<Tile, Policy>>;
    auto x_space =
        OrbitalSpaceTArray(OrbitalIndex(L"x"), OrbitalIndex(L"κ"), C_x_ta);
    orbital_registry.add(x_space);
  }

  // for now the f12 contribution is energy-independent
  TArray Sigma_pph_f12;
  {
    TArray V_ixix, V_ixxi;
    // always do DF
    std::tie(V_ixix, V_ixxi) = compute_V();

    Sigma_pph_f12("x,y") =
        0.5 * (1.25 * V_ixix("i,x,j,y") - 0.25 * V_ixxi("i,x,y,j")) *
        this->lcao_factory()(L"<i|I|j>");
  }
  // std::cout << "Sigma_pph_f12:\n" << Sigma_pph_f12 << std::endl;
  Matrix Sigma_f12 = array_ops::array_to_eigen(Sigma_pph_f12);

  // done with F12 ... remove all geminal ints and CABS indices
  this->lcao_factory().purge_operator(world, L"R");
  this->lcao_factory().purge_index(world, L"a'");
  this->lcao_factory().purge_index(world, L"ρ");

  this->lcao_factory().keep_partial_transforms(false);

  if (world.rank() == 0) {
    auto SE_F12 = SE + Sigma_f12(0, 0);
    auto unit_factory = UnitFactory::get_default();
    auto Hartree2eV = unit_factory->make_unit("eV").from_atomic_units();
    std::string orblabel =
        std::string(orbital_ < 0 ? "IP" : "EA") + std::to_string(abs(orbital_));
    std::printf("final       GF2 %6s = %11.3lf eV (%10.4lf a.u.)\n",
                orblabel.c_str(), SE * Hartree2eV, SE);
    std::printf("final GF2-F12-V %6s = %11.3lf eV (%10.4lf a.u.)\n",
                orblabel.c_str(), SE_F12 * Hartree2eV, SE_F12);
    if (orbital_ > 0)
      std::printf(
          "WARNING: non-strongly-orthogonal F12 projector is used for the F12 "
          "correction to EA!!!");
  }
}

#if TA_DEFAULT_POLICY == 1
extern template class GF2F12<TA::TensorD>;
#endif

}  // namespace lcao
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_F12_GF2F12_H_
