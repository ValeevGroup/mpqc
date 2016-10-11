//
// Created by Chong Peng on 10/10/16.
//

#include "dbmp2f12.h"

MPQC_CLASS_EXPORT_KEY2(mpqc::f12::RIDBRMP2F12, "RI-DBRMP2F12");

namespace mpqc {
namespace f12 {

using Matrix = RIDBRMP2F12::Matrix;

RIDBRMP2F12::RIDBRMP2F12(const KeyVal& kv) : LCAOWavefunction(kv) {
  db_rmp2f12_energy_ = 0.0;

  if (kv.exists("ref")) {
    ref_wfn_ = kv.keyval("ref").class_ptr<qc::Wavefunction>();
  } else {
    throw std::invalid_argument("Default Ref Wfn in RMP2 is not support! \n");
  }
}

double RIDBRMP2F12::value() {
  if (db_rmp2f12_energy_ == 0.0) {
    double ref_energy = ref_wfn_->value();

    // initialize
    auto mol = this->lcao_factory().atomic_integral().molecule();
    Eigen::VectorXd orbital_energy;
    this->trange1_engine_ = closed_shell_dualbasis_mo_build_eigen_solve_svd(
        this->lcao_factory(), orbital_energy, mol, is_frozen_core(),
        occ_block(), unocc_block());
    this->orbital_energy_ = std::make_shared<Eigen::VectorXd>(orbital_energy);

    closed_shell_dualbasis_cabs_mo_build_svd(
        this->lcao_factory(), this->trange1_engine(), "VBS", unocc_block());

    // compute
    double mp2_f12_energy = compute();

    // clear registry and recompute scf
    obsolete();
    this->wfn_world()->ao_integrals().registry().purge_index(wfn_world()->world(), OrbitalIndex(L"Α"));

    // recompute scf
    auto vbs = wfn_world()->ao_integrals().orbital_basis_registry().retrieve(OrbitalIndex(L"Α"));
    wfn_world()->ao_integrals().orbital_basis_registry().remove(OrbitalIndex(L"κ"));
    wfn_world()->ao_integrals().orbital_basis_registry().insert(OrbitalIndex(L"κ"),vbs);
    wfn_world()->ao_integrals().orbital_basis_registry().remove(OrbitalIndex(L"Α"));

    double new_ref_energy = ref_wfn_->value();

    db_rmp2f12_energy_ = ref_energy + mp2_f12_energy;
  }
  return db_rmp2f12_energy_;
}


void RIDBRMP2F12::obsolete() {
  db_rmp2f12_energy_ = 0.0;
  qc::LCAOWavefunction::obsolete();
  ref_wfn_->obsolete();
}


double RIDBRMP2F12::compute() {
  auto& world = this->wfn_world()->world();

  // start f12
  auto f12_time0 = mpqc_time::fenced_now(world);

  RowMatrix<double> mp2_eij, f12_eij;

  std::tie(mp2_eij, f12_eij) = compute_db_mp2_f12_c();

  if (world.rank() == 0) {
    auto nocc = this->trange1_engine()->get_active_occ();
    printf(
        "  i0     i1       eij(mp2)        eij(f12)      eij(mp2-f12) \n"
        "====== ====== =============== =============== ===============\n");
    for (int i = 0; i != nocc; ++i)
      for (int j = i; j != nocc; ++j)
        printf("%4d   %4d   %15.12lf %15.12lf %15.12lf\n", i, j, mp2_eij(i, j),
               f12_eij(i, j), mp2_eij(i, j) + f12_eij(i, j));
  }

  auto ef12 = f12_eij.sum();
  auto emp2 = mp2_eij.sum();
  utility::print_par(world, "E_DBMP2: ", emp2, "\n");
  utility::print_par(world, "E_DBMP2F12: ", ef12, "\n");

  // compute cabs singles
  double e_s = 0.0;
  auto single_time0 = mpqc_time::fenced_now(world);

  CABSSingles<TA::TensorD> cabs_singles(this->lcao_factory());
  e_s = cabs_singles.compute();
  utility::print_par(lcao_factory().get_world(), "E_S: ", e_s, "\n");
  auto single_time1 = mpqc_time::fenced_now(world);
  auto single_time = mpqc_time::duration_in_s(single_time0, single_time1);
  mpqc::utility::print_par(world, "Total CABS Singles Time:  ", single_time,
                           "\n");

  auto f12_time1 = mpqc_time::fenced_now(world);
  auto f12_time = mpqc_time::duration_in_s(f12_time0, f12_time1);

  mpqc::utility::print_par(world, "Total DBF12 Time:  ", f12_time, "\n");

  return emp2 + ef12 + e_s;
}

std::tuple<Matrix, Matrix> RIDBRMP2F12::compute_db_mp2_f12_c() {
  using TArray = TA::DistArray<TA::TensorD, TA::SparsePolicy>;

  auto& world = this->wfn_world()->world();

  Matrix Eij_MP2, Eij_F12;

  auto n_active_occ = this->trange1_engine()->get_active_occ();
  auto n_occ = this->trange1_engine()->get_occ();
  auto n_frozen = this->trange1_engine()->get_nfrozen();

  // create shape
  auto occ_tr1 = this->trange1_engine()->get_occ_tr1();
  TiledArray::TiledRange occ4_trange({occ_tr1, occ_tr1, occ_tr1, occ_tr1});
  auto ijij_ijji_shape = f12::make_ijij_ijji_shape(occ4_trange);

  // T Term
  TArray t2;
  {
    utility::print_par(world, "Compute T_abij With DF \n");

    TArray g_abij;
    g_abij("a,b,i,j") =
        this->lcao_factory().compute(L"<i j|G|a b>[df]")("i,j,a,b");
    t2 = mpqc::cc::d_abij(g_abij, *(this->orbital_energy()), n_occ, n_frozen);

    // compute MP2 energy and pair energies
    TArray TG_ijij_ijji;
    TG_ijij_ijji("i1,j1,i2,j2") =
        (t2("a,b,i1,j1") * g_abij("a,b,i2,j2")).set_shape(ijij_ijji_shape);
    Eij_MP2 =
        TG_ijij_ijji("i1,j1,i2,j2")
            .reduce(F12PairEnergyReductor<TA::TensorD>(2, -1, n_active_occ));
  }

  // compute V term
  TArray V_ijij_ijji =
      compute_V_ijij_ijji_db_df(lcao_factory(), ijij_ijji_shape);
  {
    // G integral in MO not needed, still need G integral in AO to compute F, K,
    // hJ
    this->lcao_factory().registry().purge_operator(world, L"G");

    // contribution from V_ijij_ijji
    Matrix eij = V_ijij_ijji("i1,j1,i2,j2")
                     .reduce(F12PairEnergyReductor<TA::TensorD>(
                         2 * C_ijij_bar, 2 * C_ijji_bar, n_active_occ));
    utility::print_par(world, "E_V: ", eij.sum(), "\n");
    Eij_F12 = eij;
  }

  // compute C term
  TArray C_ijab = compute_C_ijab_df(lcao_factory());

  {
    utility::print_par(world, "Compute CT With DF \n");
    V_ijij_ijji("i1,j1,i2,j2") =
        (C_ijab("i1,j1,a,b") * t2("a,b,i2,j2")).set_shape(ijij_ijji_shape);

    Matrix eij = V_ijij_ijji("i1,j1,i2,j2")
                     .reduce(F12PairEnergyReductor<TA::TensorD>(
                         2 * C_ijij_bar, 2 * C_ijji_bar, n_active_occ));
    utility::print_par(world, "E_CT: ", eij.sum(), "\n");
    Eij_F12 += eij;
  }

  // compute X term
  TArray X_ijij_ijji =
      compute_X_ijij_ijji_db_df(lcao_factory(), ijij_ijji_shape);
  {
    auto Fij = this->lcao_factory().compute(L"<i|F|j>[df]");
    auto Fij_eigen = array_ops::array_to_eigen(Fij);
    f12::convert_X_ijkl(X_ijij_ijji, Fij_eigen);

    Matrix eij = X_ijij_ijji("i1,j1,i2,j2")
                     .reduce(F12PairEnergyReductor<TA::TensorD>(
                         CC_ijij_bar, CC_ijji_bar, n_active_occ));
    eij *= -1;
    utility::print_par(world, "E_X: ", eij.sum(), "\n");
    Eij_F12 += eij;
  }

  // compute B term
  TArray B_ijij_ijji =
      compute_B_ijij_ijji_db_df(lcao_factory(), ijij_ijji_shape);
  {
    Matrix eij = B_ijij_ijji("i1,j1,i2,j2")
                     .reduce(F12PairEnergyReductor<TA::TensorD>(
                         CC_ijij_bar, CC_ijji_bar, n_active_occ));
    utility::print_par(world, "E_B: ", eij.sum(), "\n");
    Eij_F12 += eij;
  }

  {
    utility::print_par(world, "Compute CC Term With DF \n");
    auto C_bar_ijab =
        f12::convert_C_ijab(C_ijab, n_occ, n_frozen, *(this->orbital_energy()));
    B_ijij_ijji("i1,j1,i2,j2") = (C_ijab("i1,j1,a,b") * C_bar_ijab("i2,j2,a,b"))
                                     .set_shape(ijij_ijji_shape);

    Matrix eij = B_ijij_ijji("i1,j1,i2,j2")
                     .reduce(F12PairEnergyReductor<TA::TensorD>(
                         CC_ijij_bar, CC_ijji_bar, n_active_occ));
    utility::print_par(world, "E_CC: ", eij.sum(), "\n");
    Eij_F12 += eij;
  }

  return std::make_tuple(Eij_MP2, Eij_F12);
}

void RIDBRMP2F12::compute(qc::PropertyBase* pb) {}

}  // end of namespace f12
}  // end of namespace mpqc
