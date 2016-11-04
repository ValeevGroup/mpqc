//
// Created by Chong Peng on 10/10/16.
//

#include "dbmp2f12.h"
#include "mpqc/util/keyval/forcelink.h"
MPQC_CLASS_EXPORT2("RI-DBRMP2F12", mpqc::f12::RIDBRMP2F12);

namespace mpqc {
namespace f12 {

using TArray = RIRMP2F12::TArray;

RIDBRMP2F12::RIDBRMP2F12(const KeyVal& kv) : RIRMP2F12(kv), kv_(kv) {
  redo_mp2_ = kv.value<bool>("redo_mp2", false);
}

double RIDBRMP2F12::value() {
  if (this->energy_ == 0.0) {

    auto& world = this->wfn_world()->world();

    double time;
    auto time0 = mpqc::fenced_now(world);

    double ref_energy = ref_wfn_->value();

    auto time1 = mpqc::fenced_now(world);
    time = mpqc::duration_in_s(time0, time1);
    utility::print_par(world,"Total Ref Time: ", time, " S \n");

    // initialize
    auto mol = this->lcao_factory().ao_factory().molecule();
    Eigen::VectorXd orbital_energy;
    this->trange1_engine_ = closed_shell_dualbasis_mo_build_eigen_solve_svd(
        this->lcao_factory(), orbital_energy, mol, is_frozen_core(),
        occ_block(), unocc_block());
    this->orbital_energy_ = std::make_shared<Eigen::VectorXd>(orbital_energy);

    // create shape
    auto occ_tr1 = this->trange1_engine()->get_occ_tr1();
    TiledArray::TiledRange occ4_trange({occ_tr1, occ_tr1, occ_tr1, occ_tr1});
    ijij_ijji_shape_ = f12::make_ijij_ijji_shape(occ4_trange);

    closed_shell_dualbasis_cabs_mo_build_svd(
        this->lcao_factory(), this->trange1_engine(), "VBS", unocc_block());

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

    auto time2 = mpqc::fenced_now(world);
    time = mpqc::duration_in_s(time1, time2);
    utility::print_par(world,"Total F12 Time: ", time, " S \n");

    if (!redo_mp2_) {
      this->energy_ = ref_energy + emp2 + ef12 + e_s;
    } else {
      // clear registry
      obsolete();

      // change basis
      auto vbs = lcao_factory().ao_factory().orbital_basis_registry().retrieve(
          OrbitalIndex(L"Α"));
      lcao_factory().ao_factory().orbital_basis_registry().remove(
          OrbitalIndex(L"κ"));
      lcao_factory().ao_factory().orbital_basis_registry().insert(
          OrbitalIndex(L"κ"), vbs);
      lcao_factory().ao_factory().orbital_basis_registry().remove(
          OrbitalIndex(L"Α"));

      auto mp2 = mbpt::RIRMP2(kv_);
      double new_mp2 = mp2.value();

      this->energy_ = new_mp2 + ef12 + e_s;

      auto time3 = mpqc::fenced_now(world);
      time = mpqc::duration_in_s(time2, time3);
      utility::print_par(world,"Total New MP2 Time: ", time, " S \n");

      time = mpqc::duration_in_s(time0, time3);
      utility::print_par(world,"Total DBMP2F12 Time: ", time, " S \n");

    }

  }
  return this->energy_;
}

TArray RIDBRMP2F12::compute_B() {
  TArray B;
  if (approximation_ == 'C') {
    B = f12::compute_B_ijij_ijji_db_df(this->lcao_factory(), ijij_ijji_shape_);
  } else if (approximation_ == 'D') {
    throw std::invalid_argument(
        "Dual Basis RIRMP2F12 D Approximation Not Implemented!!");
  }
  return B;
}

TArray RIDBRMP2F12::compute_V() {
  return f12::compute_V_ijij_ijji_db_df(this->lcao_factory(), ijij_ijji_shape_);
}

TArray RIDBRMP2F12::compute_X() {
  TArray X =
      f12::compute_X_ijij_ijji_db_df(this->lcao_factory(), ijij_ijji_shape_);
  auto Fij = this->lcao_factory().compute(L"(i|F|j)[df]");
  auto Fij_eigen = array_ops::array_to_eigen(Fij);
  f12::convert_X_ijkl(X, Fij_eigen);
  return X;
}

double RIDBRMP2F12::compute_cabs_singles() {
  double es;

  if (approximation_ == 'D') {
    throw std::invalid_argument(
        "Dual Basis RIRMP2F12 D Approximation Not Implemented!!");
  } else {
    CABSSingles<TA::TensorD> cabs_singles(lcao_factory());
    es = cabs_singles.compute(true, false, true);
    // if redo mp2, this should not include the non-canonical virtual part
    if (redo_mp2_) {
      auto F_ma = this->lcao_factory().compute(L"<m|F|a>[df]");
      es -= 2 *
            F_ma("m,a").reduce(mbpt::detail::ScfCorrection<TA::TensorD>(
                this->orbital_energy(), this->trange1_engine()->get_occ()));
    }
  }

  return es;
}

}  // end of namespace f12
}  // end of namespace mpqc
