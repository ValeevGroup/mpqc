//
// Created by Chong Peng on 10/10/16.
//

#include "dbmp2f12.h"
#include "mpqc/util/keyval/forcelink.h"
#include <mpqc/chemistry/qc/scf/rhf.h>
MPQC_CLASS_EXPORT2("RI-DBRMP2F12", mpqc::f12::RIDBRMP2F12);

namespace mpqc {
namespace f12 {

using TArray = RIRMP2F12::TArray;

RIDBRMP2F12::RIDBRMP2F12(const KeyVal& kv) : RIRMP2F12(kv), kv_(kv) {
  mp2_method_ = kv.value<std::string>("mp2", "none");
  if (this->approximation_ == 'D') {
    throw std::invalid_argument(
        "Dual Basis RIRMP2F12 D Approximation Not Implemented!!");
  }
}

double RIDBRMP2F12::value() {
  if (this->energy_ == 0.0) {
    auto& world = this->wfn_world()->world();

    double time;
    auto time0 = mpqc::fenced_now(world);

    double ref_energy = ref_wfn_->value();

    auto time1 = mpqc::fenced_now(world);
    time = mpqc::duration_in_s(time0, time1);
    utility::print_par(world, "Total Ref Time: ", time, " S \n");

    // initialize
    auto mol = this->lcao_factory().ao_factory().molecule();
    Eigen::VectorXd orbital_energy;
    this->trange1_engine_ = closed_shell_dualbasis_mo_build_eigen_solve_svd(
        this->lcao_factory(), orbital_energy, mol, is_frozen_core(),
        occ_block(), unocc_block());
    this->orbital_energy_ = std::make_shared<Eigen::VectorXd>(orbital_energy);

    // create shape
    auto occ_tr1 = this->trange1_engine()->get_active_occ_tr1();
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
    utility::print_par(world, "Total F12 Time: ", time, " S \n");

    if (mp2_method_ == "none") {
      this->energy_ = ref_energy + emp2 + ef12 + e_s;
    }
    else if(mp2_method_ == "redo") {
      double new_mp2 = compute_new_mp2();
      this->energy_ = new_mp2 + ef12 + e_s;
    }
    else if(mp2_method_ == "update"){
      double new_mp2 = compute_new_mp2();
      this->energy_ = ref_energy + new_mp2 + ef12 + e_s;
    }

    auto time3 = mpqc::fenced_now(world);
    time = mpqc::duration_in_s(time0, time3);
    utility::print_par(world, "Total DBMP2F12 Time: ", time, " S \n");
  }
  return this->energy_;
}

double RIDBRMP2F12::compute_new_mp2() {
  auto& world = this->wfn_world()->world();
  auto time1 = mpqc::fenced_now(world);

  auto& lcao_factory = this->lcao_factory();

  double new_mp2 = 0.0;
  /// recompute the mp2 energy in VBS
  if (mp2_method_ == "redo") {
    auto fock = lcao_factory.ao_factory().compute(L"<Α|F|Α>[df]");

    // clear registry
    obsolete();

    // change basis
    auto vbs = lcao_factory.ao_factory().orbital_basis_registry().retrieve(
        OrbitalIndex(L"Α"));
    lcao_factory.ao_factory().orbital_basis_registry().remove(
        OrbitalIndex(L"κ"));
    lcao_factory.ao_factory().orbital_basis_registry().insert(
        OrbitalIndex(L"κ"), vbs);
    lcao_factory.ao_factory().orbital_basis_registry().remove(
        OrbitalIndex(L"Α"));

    auto mp2 = mbpt::RIRMP2(kv_);

    auto ref = mp2.refwfn();

    std::shared_ptr<scf::RHF<TA::TensorD,TA::SparsePolicy>> rhf = std::dynamic_pointer_cast<scf::RHF<TA::TensorD,TA::SparsePolicy>>(ref);

    rhf->set_fock(fock);

    new_mp2 = mp2.value();
  } else if (mp2_method_ == "update") {
    std::size_t n = this->trange1_engine()->get_all();
    std::size_t v = this->trange1_engine()->get_vir();
    std::size_t o = this->trange1_engine()->get_occ();
    std::size_t f = this->trange1_engine()->get_nfrozen();
    // get fock matrix
    RowMatrixXd fock;
    {
      auto f_ii = lcao_factory.compute(L"<m|F|m>[df]");
      auto f_aa = lcao_factory.compute(L"<a|F|a>[df]");
      auto f_ia = lcao_factory.compute(L"<m|F|a>[df]");

      RowMatrixXd f_ii_eigen = array_ops::array_to_eigen(f_ii);
      RowMatrixXd f_aa_eigen = array_ops::array_to_eigen(f_aa);
      RowMatrixXd f_ia_eigen = array_ops::array_to_eigen(f_ia);

      // form the new fock matrix
      fock = RowMatrixXd::Zero(n, n);
      fock.block(0, 0, o, o) << f_ii_eigen;
      fock.block(o, o, v, v) << f_aa_eigen;
      fock.block(0, o, o, v) << f_ia_eigen;
      fock.block(o, 0, v, o) << f_ia_eigen.transpose();
    }

    // solve new orbitals
    TArray C_occ_ta, C_corr_occ_ta, C_vir_ta;
    Eigen::VectorXd ens_all;
    {
      // get the basis map
      Eigen::RowVectorXi basis_map;
      {
        auto vbs =
            this->wfn_world()->basis_registry()->retrieve(OrbitalIndex(L"Α"));
        auto obs =
            this->wfn_world()->basis_registry()->retrieve(OrbitalIndex(L"κ"));
        basis_map = basis::sub_basis_map(vbs, obs);

//        std::cout << "Basis Map" << basis_map << std::endl;
      }

      // update old C_occ to vbs
      RowMatrixXd old_coeffs = RowMatrixXd::Zero(n, n);
      {
        auto old_occ =
            lcao_factory.orbital_space().retrieve(OrbitalIndex(L"m")).array();
        RowMatrixXd old_occ_eigen = array_ops::array_to_eigen(old_occ);

        // convert occ to vbs
        RowMatrixXd occ_eigen = RowMatrixXd::Zero(n, o);
        for (std::size_t i = 0; i < n; i++) {
          if (basis_map[i] != 0) {
            occ_eigen.row(i) = old_occ_eigen.row(basis_map[i] - 1);
          }
        }

        auto old_vir =
            lcao_factory.orbital_space().retrieve(OrbitalIndex(L"a")).array();
        RowMatrixXd old_vir_eigen = array_ops::array_to_eigen(old_vir);

        old_coeffs.block(0, 0, n, o) << occ_eigen;
        old_coeffs.block(0, o, n, v) << old_vir_eigen;
      }

      // eigen solve
      Eigen::SelfAdjointEigenSolver<RowMatrixXd> es(fock);
      ens_all = es.eigenvalues();
      RowMatrixXd C_all = es.eigenvectors();
      utility::print_par(world, "New Orbitals\n", ens_all, "\n");

      // update the coefficient

      C_all = old_coeffs * C_all;

      RowMatrixXd C_occ = C_all.block(0, 0, n, o);
      RowMatrixXd C_corr_occ = C_all.block(0, f, n, o - f);
      RowMatrixXd C_vir = C_all.rightCols(n - o);

      auto tr_vbs = this->wfn_world()
                        ->basis_registry()
                        ->retrieve(OrbitalIndex(L"Α"))
                        .create_trange1();
      auto tr_corr_occ = this->trange1_engine()->get_active_occ_tr1();
      auto tr_occ = this->trange1_engine()->get_occ_tr1();
      auto tr_vir = this->trange1_engine()->get_vir_tr1();

      // convert to TA
      C_occ_ta =
          array_ops::eigen_to_array<TA::TensorD>(world, C_occ, tr_vbs, tr_occ);
      C_corr_occ_ta = array_ops::eigen_to_array<TA::TensorD>(
          world, C_corr_occ, tr_vbs, tr_corr_occ);
      C_vir_ta =
          array_ops::eigen_to_array<TA::TensorD>(world, C_vir, tr_vbs, tr_vir);
    }

    lcao_factory.orbital_space().clear();
    lcao_factory.registry().purge(world);
    lcao_factory.ao_factory().registry().purge_operator(world,L"F");
    lcao_factory.ao_factory().registry().purge_operator(world,L"J");
    lcao_factory.ao_factory().registry().purge_operator(world,L"K");

    auto& orbital_registry = lcao_factory.orbital_space();

    // insert to registry
    using OrbitalSpaceTArray =
        OrbitalSpace<TA::DistArray<TA::TensorD, TA::SparsePolicy>>;
    auto occ_space =
        OrbitalSpaceTArray(OrbitalIndex(L"m"), OrbitalIndex(L"Α"), C_occ_ta);
    orbital_registry.add(occ_space);

    auto corr_occ_space = OrbitalSpaceTArray(OrbitalIndex(L"i"),
                                             OrbitalIndex(L"Α"), C_corr_occ_ta);
    orbital_registry.add(corr_occ_space);

    auto vir_space =
        OrbitalSpaceTArray(OrbitalIndex(L"a"), OrbitalIndex(L"Α"), C_vir_ta);
    orbital_registry.add(vir_space);

    auto g_ijab = lcao_factory.compute(L"<i j|G|a b>[df]");

    new_mp2 =
        (g_ijab("i,j,a,b") * (2 * g_ijab("i,j,a,b") - g_ijab("i,j,b,a")))
            .reduce(mbpt::detail::Mp2Energy<TA::TensorD>(
                std::make_shared<Eigen::VectorXd>(ens_all),
//                this->orbital_energy_,
                trange1_engine_->get_occ(), trange1_engine_->get_nfrozen()));

    utility::print_par(world, "New MP2: ", new_mp2, "\n");

//    double d_scf = ens_all.segment(0,o).sum() - this->orbital_energy_->segment(0,o).sum();
//    utility::print_par(world, "SCF Correction: ", d_scf, "\n");

//    auto f_ia = lcao_factory.compute(L"<m|F|a>[df]");
//    d_scf += 2*f_ia("i,a").reduce(mbpt::detail::ScfCorrection<TA::TensorD>(std::make_shared<Eigen::VectorXd>(ens_all), o));

//    utility::print_par(world, "SCF Correction: ", d_scf, "\n");

//    new_mp2 += d_scf;

  } else {
    utility::print_par(
        world, "\n Warning! Not valid method to update MP2 energy! Skip! \n");
  }

  auto time2 = mpqc::fenced_now(world);
  double time = mpqc::duration_in_s(time1, time2);
  utility::print_par(world, "Total New MP2 Time: ", time, " S \n");

  return new_mp2;
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
//     if redo mp2, this should not include the non-canonical virtual part
//    if (mp2_method_ == "redo" || mp2_method_ == "update") {
    if (mp2_method_ == "redo") {
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
