//
// Created by Chong Peng on 10/13/16.
//

#include "mp2f12.h"

MPQC_CLASS_EXPORT_KEY2(mpqc::f12::RMP2F12, "RMP2F12");
MPQC_CLASS_EXPORT_KEY2(mpqc::f12::RIRMP2F12, "RI-RMP2F12");

namespace mpqc {
namespace f12 {

using TArray = RMP2F12::TArray;
using Matrix = RMP2F12::Matrix;

RMP2F12::RMP2F12(const KeyVal& kv) : LCAOWavefunction(kv) {
  if (kv.exists("ref")) {
    ref_wfn_ = kv.keyval("ref").class_ptr<qc::Wavefunction>();
  } else {
    throw std::invalid_argument(
        "Default Ref Wfn in RMP2F12 is not support! \n");
  }

  rmp2f12_energy_ = 0.0;

  approximation_ = kv.value<char>("approaximation", 'C');
  if (approximation_ != 'C' && approximation_ != 'D') {
    throw std::invalid_argument("Only approaximation C or D is supported!");
  }

  cabs_singles_ = kv.value<bool>("cabs_singles", true);

}

double RMP2F12::value() {
  if (rmp2f12_energy_ == 0.0) {
    auto& world = this->wfn_world()->world();

    double time;
    auto time0 = mpqc_time::fenced_now(world);

    double ref_energy = ref_wfn_->value();

    auto time1 = mpqc_time::fenced_now(world);
    time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world,"Total Ref Time: ", time, " S \n");

    // initialize
    auto mol = this->lcao_factory().atomic_integral().molecule();
    Eigen::VectorXd orbital_energy;
    this->trange1_engine_ = closed_shell_obs_mo_build_eigen_solve(
        this->lcao_factory(), orbital_energy, mol, is_frozen_core(),
        occ_block(), unocc_block());

    this->orbital_energy_ = std::make_shared<Eigen::VectorXd>(orbital_energy);

    // create shap
    auto occ_tr1 = this->trange1_engine()->get_occ_tr1();
    TiledArray::TiledRange occ4_trange({occ_tr1, occ_tr1, occ_tr1, occ_tr1});
    ijij_ijji_shape_ = f12::make_ijij_ijji_shape(occ4_trange);

    closed_shell_cabs_mo_build_svd(this->lcao_factory(), this->trange1_engine(),
                                   unocc_block());

    // compute
    Matrix mp2_eij, f12_eij;
    std::tie(mp2_eij, f12_eij) = compute();

    if (world.rank() == 0) {
      auto nocc = this->trange1_engine()->get_active_occ();
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
    utility::print_par(lcao_factory().world(), "E_MP2: ", emp2, "\n");
    utility::print_par(lcao_factory().world(), "E_F12: ", ef12, "\n");

    double e_s;
    if (cabs_singles_) {
      e_s = compute_cabs_singles();
    }

    utility::print_par(world, "E_S: ", e_s, "\n");

    rmp2f12_energy_ = ref_energy + emp2 + ef12 + e_s;

    auto time2 = mpqc_time::fenced_now(world);
    time = mpqc_time::duration_in_s(time1, time2);
    utility::print_par(world,"Total F12 Time: ", time, " S \n");

    time = mpqc_time::duration_in_s(time0, time2);
    utility::print_par(world,"Total MP2F12 Time: ", time, " S \n");
  }

  return rmp2f12_energy_;
}

void RMP2F12::obsolete() {
  rmp2f12_energy_ = 0.0;
  qc::LCAOWavefunction::obsolete();
  ref_wfn_->obsolete();
}

void RMP2F12::compute(qc::PropertyBase* pb) {}

std::tuple<Matrix, Matrix> RMP2F12::compute() {
  auto& world = lcao_factory().world();

  utility::print_par(world, "\n Computing MP2F12 ", approximation_,
                     " Approximation \n");

  Matrix Eij_MP2, Eij_F12;

  auto n_active_occ = this->trange1_engine()->get_active_occ();
  auto n_occ = this->trange1_engine()->get_occ();
  auto n_frozen = this->trange1_engine()->get_nfrozen();

  auto ijij_ijji_shape = ijij_ijji_shape_;

  // compute B term
  {
    TArray B_ijij_ijji = compute_B();
    Matrix Eij_b = B_ijij_ijji("i1,j1,i2,j2")
                       .reduce(F12PairEnergyReductor<TA::TensorD>(
                           CC_ijij_bar, CC_ijji_bar, n_active_occ));
    utility::print_par(world, "E_B: ", Eij_b.sum(), "\n");
    Eij_F12 = Eij_b;
  }

  // compute X term
  {
    TArray X_ijij_ijji = compute_X();
    lcao_factory().purge_operator(world, L"R2");

    Matrix Eij_x = X_ijij_ijji("i1,j1,i2,j2")
                       .reduce(F12PairEnergyReductor<TA::TensorD>(
                           CC_ijij_bar, CC_ijji_bar, n_active_occ));
    Eij_x *= -1.0;
    utility::print_par(world, "E_X: ", Eij_x.sum(), "\n");
    Eij_F12 += Eij_x;
  }

  // compute V term
  {
    TArray V_ijij_ijji = compute_V();
    // G integral in MO not needed, still need G integral in AO to compute F, K,
    // hJ
    lcao_factory().registry().purge_operator(world, L"G");
    lcao_factory().purge_operator(world, L"GR");

    // contribution from V_ijij_ijji
    // NB factor of 2 from the Hylleraas functional
    Matrix e_ij = V_ijij_ijji("i1,j1,i2,j2")
                      .reduce(F12PairEnergyReductor<TA::TensorD>(
                          2 * C_ijij_bar, 2 * C_ijji_bar, n_active_occ));
    Eij_F12 += e_ij;
    utility::print_par(world, "E_V: ", e_ij.sum(), "\n");
  }

  TArray t2;
  {
    TArray g_abij;
    std::tie(t2, g_abij) = compute_T();

    // compute MP2 energy and pair energies
    TArray TG_ijij_ijji;
    TG_ijij_ijji("i1,j1,i2,j2") =
        (t2("a,b,i1,j1") * g_abij("a,b,i2,j2")).set_shape(ijij_ijji_shape);
    Eij_MP2 =
        TG_ijij_ijji("i1,j1,i2,j2")
            .reduce(F12PairEnergyReductor<TA::TensorD>(2, -1, n_active_occ));
  }

  // compute C term
  TArray C_ijab = compute_C();

  {
    utility::print_par(world, "Compute CT \n");
    TArray CT;
    CT("i1,j1,i2,j2") =
        (C_ijab("i1,j1,a,b") * t2("a,b,i2,j2")).set_shape(ijij_ijji_shape);

    // NB factor of 2 from the Hylleraas functional
    Matrix Eij_ct = CT("i1,j1,i2,j2")
                        .reduce(F12PairEnergyReductor<TA::TensorD>(
                            2 * C_ijij_bar, 2 * C_ijji_bar, n_active_occ));
    utility::print_par(world, "E_CT: ", Eij_ct.sum(), "\n");
    Eij_F12 += Eij_ct;
  }

  {
    utility::print_par(world, "Compute CC Term \n");
    TArray CC;
    auto C_bar_ijab =
        f12::convert_C_ijab(C_ijab, n_occ, n_frozen, *(this->orbital_energy()));
    CC("i1,j1,i2,j2") = (C_ijab("i1,j1,a,b") * C_bar_ijab("i2,j2,a,b"))
                            .set_shape(ijij_ijji_shape);

    Matrix Eij_cc = CC("i1,j1,i2,j2")
                        .reduce(F12PairEnergyReductor<TA::TensorD>(
                            CC_ijij_bar, CC_ijji_bar, n_active_occ));
    utility::print_par(world, "E_CC: ", Eij_cc.sum(), "\n");
    Eij_F12 += Eij_cc;
  }

  return std::make_tuple(Eij_MP2, Eij_F12);
}

TArray RMP2F12::compute_B() {
  TArray B;
  if (approximation_ == 'C') {
    B = f12::compute_B_ijij_ijji_C(this->lcao_factory(), ijij_ijji_shape_);
  } else if (approximation_ == 'D') {
    B = f12::compute_B_ijij_ijji_D(this->lcao_factory(), ijij_ijji_shape_);
  }
  return B;
}

TArray RMP2F12::compute_C() {
  return f12::compute_C_ijab(this->lcao_factory());
}

TArray RMP2F12::compute_V() {
  return f12::compute_V_ijij_ijji(this->lcao_factory(), ijij_ijji_shape_);
}

TArray RMP2F12::compute_X() {
  TArray X = f12::compute_X_ijij_ijji(this->lcao_factory(), ijij_ijji_shape_);
  auto Fij = this->lcao_factory().compute(L"(i|F|j)");
  auto Fij_eigen = array_ops::array_to_eigen(Fij);
  f12::convert_X_ijkl(X, Fij_eigen);
  return X;
}

std::tuple<TArray, TArray> RMP2F12::compute_T() {
  TArray g_abij, t2;
  g_abij("a,b,i,j") = lcao_factory().compute(L"<i j|G|a b>")("i,j,a,b");
  t2 = mpqc::cc::d_abij(g_abij, *(this->orbital_energy()),
                        this->trange1_engine()->get_occ(),
                        this->trange1_engine()->get_nfrozen());

  return std::tuple<TArray, TArray>(t2, g_abij);
}

double RMP2F12::compute_cabs_singles() {
  double es;

  if (approximation_ == 'D') {
    CABSSingles<TA::TensorD> cabs_singles(lcao_factory());
    es = cabs_singles.compute(false, true, true);
  } else {
    CABSSingles<TA::TensorD> cabs_singles(lcao_factory());
    es = cabs_singles.compute(false, false, true);
  }

  return es;
}


RIRMP2F12::RIRMP2F12(const KeyVal &kv) : RMP2F12(kv){}

TArray RIRMP2F12::compute_B() {
  TArray B;
  if (approximation_ == 'C') {
    B = f12::compute_B_ijij_ijji_C_df(this->lcao_factory(), ijij_ijji_shape_);
  } else if (approximation_ == 'D') {
    B = f12::compute_B_ijij_ijji_D_df(this->lcao_factory(), ijij_ijji_shape_);
  }
  return B;
}

TArray RIRMP2F12::compute_C() {
  return f12::compute_C_ijab_df(this->lcao_factory());
}

TArray RIRMP2F12::compute_V() {
  return f12::compute_V_ijij_ijji_df(this->lcao_factory(), ijij_ijji_shape_);
}

TArray RIRMP2F12::compute_X() {
  TArray X = f12::compute_X_ijij_ijji_df(this->lcao_factory(), ijij_ijji_shape_);
  auto Fij = this->lcao_factory().compute(L"(i|F|j)[df]");
  auto Fij_eigen = array_ops::array_to_eigen(Fij);
  f12::convert_X_ijkl(X, Fij_eigen);
  return X;
}

std::tuple<TArray, TArray> RIRMP2F12::compute_T() {
  TArray g_abij, t2;
  g_abij("a,b,i,j") = lcao_factory().compute(L"<i j|G|a b>[df]")("i,j,a,b");
  t2 = mpqc::cc::d_abij(g_abij, *(this->orbital_energy()),
                        this->trange1_engine()->get_occ(),
                        this->trange1_engine()->get_nfrozen());

  return std::tuple<TArray, TArray>(t2, g_abij);
}

double RIRMP2F12::compute_cabs_singles() {
  double es;

  if (approximation_ == 'D') {
    CABSSingles<TA::TensorD> cabs_singles(lcao_factory());
    es = cabs_singles.compute(true, true, true);
  } else {
    CABSSingles<TA::TensorD> cabs_singles(lcao_factory());
    es = cabs_singles.compute(true, false, true);
  }

  return es;
}

} // end of namespace f12
} // end of namespace mpqc
