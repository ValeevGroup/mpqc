/*
 *
 */

#include "mpqc/chemistry/qc/lcao/ci/cis_d.h"

namespace mpqc {
namespace lcao {

template <typename Tile, typename Policy>
cc::Intermediates<TA::DistArray<Tile, Policy>>
EOM_CCSD<Tile, Policy>::compute_FWintermediates() {
  ExEnv::out0() << indent << "\nInitialize Intermediates in EOM-CCSD\n";

  auto& world = this->wfn_world()->world();
  bool accurate_time = this->lcao_factory().accurate_time();
  auto time0 = mpqc::now(world, accurate_time);

  bool df = this->is_df();
  TArray Tau;
  auto t2 = this->t2();
  auto t1 = this->t1();
  Tau("a,b,i,j") = t2("a,b,i,j") + (t1("a,i") * t1("b,j"));

  cc::Intermediates<TArray> imds;

  if (df) {
    imds.Xia = this->lcao_factory().compute(L"( Λ |G|i a)[inv_sqr]");
    imds.Xab = this->lcao_factory().compute(L"( Λ |G|a b)[inv_sqr]");
  }

  std::tie(imds.FAB, imds.FIA, imds.FIJ) = cc::compute_cs_ccsd_F(
      this->lcao_factory(), this->ao_factory(), t1, Tau, df);

  // \cal{W}abef
  if (this->method_ != "direct" && this->method_ != "direct_df") {
    imds.Wabcd = cc::compute_cs_ccsd_W_AbCd(this->lcao_factory(), t1, Tau, df);
  }

  // \cal{W}mnij
  imds.Wijkl = cc::compute_cs_ccsd_W_KlIj(this->lcao_factory(), t1, Tau, df);

  // \cal{W}mbij
  imds.Wiajk = cc::compute_cs_ccsd_W_KaIj(this->lcao_factory(), t1, t2, Tau,
                                          imds.FIA, imds.Wijkl, df);
  // \cal{W}abei
  imds.Wabci =
      cc::compute_cs_ccsd_W_AbCi(this->lcao_factory(), this->ao_factory(), t1,
                                 t2, Tau, imds.FIA, imds.Wabcd, df);
  if (!df) {
    // \cal{W}amef
    imds.Waibc = cc::compute_cs_ccsd_W_AkCd(this->lcao_factory(), t1, df);
  }

  // \cal{W}mbej
  imds.Wiabj = cc::compute_cs_ccsd_W_IbAj(this->lcao_factory(), t1, t2, df);

  imds.Wiajb = cc::compute_cs_ccsd_W_IaJb(this->lcao_factory(), t1, t2, df);

  // \cal{W}mnie
  imds.Wijka = cc::compute_cs_ccsd_W_KlIc(this->lcao_factory(), t1, df);

  std::wstring postfix = df ? L"[df]" : L"";
  imds.Wijab = this->lcao_factory().compute(L"<i j|G|a b>" + postfix);

  auto time1 = mpqc::now(world, accurate_time);
  auto time = mpqc::duration_in_s(time0, time1);

  ExEnv::out0() << indent
                << "\nTime to Initialize Intermediates in EOM-CCSD: " << time
                << " S\n";
  return imds;
}

// compute [HSS_HSD C]^A_I
template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> EOM_CCSD<Tile, Policy>::compute_HSS_HSD_C(
    const TArray& Cai, const TArray& Cabij,
    const cc::Intermediates<TArray>& imds) {
  TArray HSS_HSD_C;

  {
    // HSS * C part
    HSS_HSD_C("a,i") =  //   Fac C^c_i
        imds.FAB("a,c") * Cai("c,i")
        // - Fki C^a_k
        - imds.FIJ("k,i") * Cai("a,k")
        // + Wakic C^c_k
        + (2.0 * imds.Wiabj("k,a,c,i") + imds.Wiajb("k,c,i,a")) * Cai("c,k");
  }

  // HSD * C part
  {
    TArray C;
    C("a,c,i,k") = 2.0 * Cabij("a,c,i,k") - Cabij("c,a,i,k");
    HSS_HSD_C("a,i") +=  //   Fkc C^ac_ik
        imds.FIA("k,c") * C("a,c,i,k")
        // - 1/2 Wklic C^ac_kl
        - imds.Wijka("k,l,i,c") * C("a,c,k,l");

    if (this->df_) {
      auto t1 = this->t1();
      // + 1/2 Wakcd C^cd_ik

      HSS_HSD_C("a,i") += imds.Xia("K,k,d") * C("c,d,i,k") *
                          (imds.Xab("K,a,c") - t1("a,l") * imds.Xia("K,l,c"));

    } else {
      // + 1/2 Wakcd C^cd_ik
      HSS_HSD_C("a,i") += imds.Waibc("a,k,c,d") * C("c,d,i,k");
    }
  }

  return HSS_HSD_C;
}

// compute [HDS_HDD C]^Ab_Ij
template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> EOM_CCSD<Tile, Policy>::compute_HDS_HDD_C(
    const TArray& Cai, const TArray& Cabij,
    const cc::Intermediates<TArray>& imds) {
  TArray HDS_HDD_C;
  auto t2 = this->t2();
  auto t1 = this->t1();

  // HDS * C part
  {
    HDS_HDD_C("a,b,i,j") =
        // P(ij) Wabcj C^c_i
        // WAbCj C^C_I + WbAcI
        imds.Wabci("a,b,c,j") * Cai("c,i")
        // P(ab) Wkaij C^b_k C^c_j
        // - WkAjI C^b_k - WKbIj C^A_K
        - imds.Wiajk("k,b,i,j") * Cai("a,k")

        // - P(ij) Wlkjc C^c_k t^ab_il
        // - WlKjC C^C_K T^Ab_Il - Wlkjc C^c_k T^Ad_Ij
        // + WLKIC C^C_K T^Ab_jL + WLkIc C^c_k T^Ab_jL
        - (2.0 * imds.Wijka("l,k,j,c") - imds.Wijka("k,l,j,c")) * Cai("c,k") *
              t2("a,b,i,l");

    if (this->df_) {
      TArray tmp;
      tmp("K,b,d") = imds.Xab("K,b,d") - t1("b,l") * imds.Xia("K,l,d");

      HDS_HDD_C("a,b,i,j") +=
          (2.0 * imds.Xia("K,k,c") * Cai("c,k") * tmp("K,b,d") -
           tmp("K,b,c") * Cai("c,k") * imds.Xia("K,k,d")) *
          t2("a,d,i,j");

    } else {
      // + P(ab) Wbkdc C^c_k T^ad_ij
      // + WbKdC C^C_K T^Ad_Ij + Wbkdc C^c_k T^Ad_Ij
      // - WAKDC C^C_K T^bD_Ij - WAkDc C^c_k T^bD_Ij
      HDS_HDD_C("a,b,i,j") +=
          (2.0 * imds.Waibc("b,k,d,c") - imds.Waibc("b,k,c,d")) * Cai("c,k") *
          t2("a,d,i,j");
    }
    //    HDS_HDD_C("a,b,i,j") += HDS_HDD_C("b,a,j,i");
  }

  // HDD * C part
  {
    TArray C, GC_ab, GC_ij;
    C("a,c,i,k") = 2.0 * Cabij("a,c,i,k") - Cabij("c,a,i,k");
    GC_ab("a,b") = imds.Wijab("k,l,a,c") * C("b,c,k,l");
    GC_ij("i,j") = imds.Wijab("i,k,c,d") * C("c,d,j,k");

    TArray tmp;
    HDS_HDD_C("a,b,i,j") +=

        //   P(ab) Fbc C^ac_ij
        //   Fbc C^ac_ij + Fac C^cb_ij
        imds.FAB("a,c") * Cabij("c,b,i,j")

        // - P(ij) Fkj C^ab_ik
        // - Fkj C^ab_ik - Fki C^ab_kj
        - imds.FIJ("k,j") * Cabij("a,b,i,k")

        // + P(ab) P(ij) Wbkjc C^ac_ik
        // + Wbkjc C^ac_ik - Wbkic C^ac_jk
        // - Wakjc C^bc_ik + Wakic C^bc_jk
        + imds.Wiabj("k,b,c,j") * C("a,c,i,k") +
        imds.Wiajb("k,c,j,b") * Cabij("a,c,i,k") +
        imds.Wiajb("k,c,i,b") * Cabij("a,c,k,j")

        // - 1/2 P(ab) g^lk_dc C^ca_kl t^db_ij
        // - 1/2 g^kl_dc C^ac_kl t^db_ij
        // - 1/2 g^kl_cd C^cb_kl t^ad_ij
        - GC_ab("d,a") * t2("d,b,i,j")

        // + 1/2 P(ij) Wlkdc C^dc_ik t^ab_jl
        // - 1/2 Wlkcd C^cd_ik t^ab_lj
        // - 1/2 Wlkcd C^cd_jk t^ab_il
        - GC_ij("l,j") * t2("a,b,i,l");

    HDS_HDD_C("a,b,i,j") += HDS_HDD_C("b,a,j,i")
                            // + 1/2 Wabcd C^cd_ij
                            //        + WAbCd_("a,b,c,d") * Cabij("c,d,i,j")
                            // + 1/2 Wklij C^ab_kl
                            + imds.Wijkl("k,l,i,j") * Cabij("a,b,k,l");

    if (imds.Wabcd.is_initialized()) {
      // + 1/2 Wabcd C^cd_ij
      //        + WAbCd_("a,b,c,d") * Cabij("c,d,i,j")
      HDS_HDD_C("a,b,i,j") += imds.Wabcd("a,b,c,d") * Cabij("c,d,i,j");
    } else {
      TArray tau;
      tau("a,b,i,j") = t2("a,b,i,j") + t1("a,i") * t1("b,j");

      // integral direct term
      auto direct_integral = this->get_direct_ao_integral();

      auto Ca =
          this->lcao_factory().orbital_registry().retrieve(OrbitalIndex(L"a"));
      auto Ci =
          this->lcao_factory().orbital_registry().retrieve(OrbitalIndex(L"i"));

      TArray U;
      U("p,r,i,j") =
          Cabij("c,d,i,j") * Ca("q,c") * Ca("s,d") * direct_integral("p,q,r,s");
      //      U("p,r,i,j") = 0.5 * (U("p,r,i,j") + U("r,p,j,i"));
      HDS_HDD_C("a,b,i,j") +=
          U("p,r,i,j") * Ca("p,a") * Ca("r,b") -
          U("r,p,i,j") * Ci("p,k") * Ca("r,a") * t1("b,k") -
          U("p,r,i,j") * Ci("p,k") * Ca("r,b") * t1("a,k") +
          U("p,r,i,j") * Ci("p,k") * Ci("r,l") * tau("a,b,k,l");
    }
  }

  return HDS_HDD_C;
}

template <typename Tile, typename Policy>
EigenVector<typename Tile::numeric_type>
EOM_CCSD<Tile, Policy>::eom_ccsd_davidson_solver(
    std::size_t n_roots, const std::vector<TArray>& cis_vector,
    const std::vector<numeric_type>& cis_eigs, std::size_t max_iter,
    double convergence) {
  cc::Intermediates<TArray> imds = init();

  std::size_t n_guess = cis_vector.size();
  std::vector<GuessVector> C(n_guess);
  auto t2 = this->t2();
  for (std::size_t i = 0; i < n_guess; i++) {
    C[i] = ::mpqc::cc::TPack<TArray>(2);
    C[i][0]("a,i") = cis_vector[i]("i,a");
    C[i][1] = TArray(t2.world(), t2.trange(), t2.shape());
    C[i][1].fill(0.0);
  }

  /// make preconditioner
  std::shared_ptr<DavidsonDiagPred<GuessVector>> pred;
  {
    EigenVector<numeric_type> eps_o =
        array_ops::array_to_eigen(imds.FIJ).diagonal();

    RowMatrix<numeric_type> FAB_eigen = array_ops::array_to_eigen(imds.FAB);
    EigenVector<numeric_type> eps_v = FAB_eigen.diagonal();

    pred = std::make_shared<cc::EEPred<TArray>>(eps_o, eps_v);

    // if simulate PNO, need to compute guess and initialize PNO, OSV
    if (!eom_pno_.empty()) {
      //       make initial guess, by run 2 iteration of
      //      DavidsonDiag<GuessVector> dvd(n_roots, false, 2, max_vector_,
      //                                    vector_threshold_);
      //
      //      for (std::size_t i = 0; i < 2; i++) {
      //        std::size_t dim = C.size();
      //        std::vector<GuessVector> HC(dim);
      //        for (std::size_t i = 0; i < dim; ++i) {
      //          if (C[i].t1.is_initialized() && C[i].t2.is_initialized()) {
      //            HC[i].t1 = compute_HSS_HSD_C(C[i].t1, C[i].t2);
      //            HC[i].t2 = compute_HDS_HDD_C(C[i].t1, C[i].t2);
      //
      //          } else {
      //            throw ProgrammingError("Guess Vector not initialized",
      //            __FILE__,
      //                                   __LINE__);
      //          }
      //        }
      //        EigenVector<numeric_type> eig_new, norms;
      //        std::tie(eig_new, norms) = dvd.extrapolate(HC, C, pred.get());
      //      }
      //
      //      C = dvd.eigen_vector();

      //      std::vector<TArray> guess(n_roots);
      //      for (std::size_t i = 0; i < n_roots; i++) {
      //        guess[i] = C[i].t2;
      //                        std::cout << guess[i] << std::endl;
      //      }

      // make CIS(D) doubles amplitudes

      if (eom_pno_ == "default") {
        auto guess = compute_cis_d_double_amplitude(
            this->lcao_factory(), cis_vector[0], cis_eigs[0], this->is_df());

        pred = std::make_shared<cc::PNOEEPred<TArray>>(
            guess, n_roots, eps_o, FAB_eigen, eom_tpno_, eom_tosv_,
            eom_pno_canonical_);
      } else if (eom_pno_ == "state-average") {
        std::vector<TArray> guess(n_roots);
        for (std::size_t i = 0; i < n_roots; i++) {
          guess[i] = compute_cis_d_double_amplitude(
              this->lcao_factory(), cis_vector[i], cis_eigs[i], this->is_df());
        }

        pred = std::make_shared<cc::StateAveragePNOEEPred<TArray>>(
            guess, n_roots, eps_o, FAB_eigen, eom_tpno_, eom_tosv_,
            eom_pno_canonical_);
      } else if(eom_pno_ == "state-merged"){
        std::vector<TArray> guess(n_roots);
        for (std::size_t i = 0; i < n_roots; i++) {
          guess[i] = compute_cis_d_double_amplitude(
              this->lcao_factory(), cis_vector[i], cis_eigs[i], this->is_df());
        }

        pred = std::make_shared<cc::StateMergedPNOEEPred<TArray>>(
            guess, n_roots, eps_o, FAB_eigen, eom_tpno_, eom_tosv_,
            eom_pno_canonical_);
      }
      else {
        throw InputError("Invalid option!", __FILE__, __LINE__, "eom_pno",
                         eom_pno_.c_str());
      }
    }
  }

  auto op = [this, &imds](const std::vector<GuessVector>& vec) {
    std::size_t dim = vec.size();

    // compute product of H with guess vector
    std::vector<GuessVector> HC(dim, GuessVector(2));
    for (std::size_t i = 0; i < dim; ++i) {
      if (vec[i][0].is_initialized() && vec[i][1].is_initialized()) {
        HC[i][0] = compute_HSS_HSD_C(vec[i][0], vec[i][1], imds);
        HC[i][1] = compute_HDS_HDD_C(vec[i][0], vec[i][1], imds);

      } else {
        throw ProgrammingError("Guess Vector not initialized", __FILE__,
                               __LINE__);
      }
    }

    return HC;
  };

  EigenVector<numeric_type> eig;
  std::vector<GuessVector> eig_vector;
  if (davidson_solver_ == "multi-state") {
    /// make davidson object
    DavidsonDiag<GuessVector> dvd(n_roots, false, 2, max_vector_,
                                  vector_threshold_);
    eig = dvd.solve(C, op, pred.get(), convergence, max_iter);
    eig_vector = dvd.eigen_vector();
  }
  // single-state
  else {
    // choose the shift according to the following paper

    /// Helmich, B. and Hättig, C. (2013) The Journal of Chemical Physics.
    /// 139(8), p. 84114. doi: 10.1063/1.4819071.
    double shift = 1.5 * cis_eigs[n_roots - 1];

    /// make davidson object
    SingleStateDavidsonDiag<GuessVector> dvd(n_roots, shift, false, 2,
                                             max_vector_, vector_threshold_);
    eig = dvd.solve(C, op, pred.get(), convergence, max_iter);
    eig_vector = dvd.eigen_vector();
  }

  ExEnv::out0() << "\n";
  util::print_excitation_energy(eig, false);

  for (std::size_t i = 0; i < n_roots; i++) {
    auto t1_dominants = array_abs_max_n_index(eig_vector[i][0], 5);

    /// double the size due to permutation symmetry in t2
    auto t2_dominants = array_abs_max_n_index(eig_vector[i][1], 5);

    ExEnv::out0() << "Dominant determinants of excited wave function " << i + 1
                  << "\n";
    util::print_t1_dominant_elements(t1_dominants);
    ExEnv::out0() << "\n";
    util::print_t2_dominant_elements(t2_dominants);
    ExEnv::out0() << "\n";
  }

  return eig;
}

template <typename Tile, typename Policy>
void EOM_CCSD<Tile, Policy>::evaluate(ExcitationEnergy* ex_energy) {
  auto target_precision = ex_energy->target_precision(0);
  auto ccsd_precision = 0.1 * target_precision;
  if (vector_threshold_ == 0) {
    vector_threshold_ = 10 * target_precision;
  }
  if (!this->computed()) {
    auto& world = this->wfn_world()->world();

    bool triplets = ex_energy->triplets();

    if (triplets) {
      throw InputError("EOM-CCSD only supports Singlets at this moment!\n",
                       __FILE__, __LINE__, "triplet", "true");
    }

    auto ccsd_energy =
        std::make_shared<Energy>(this->shared_from_this(), ccsd_precision);
    // do CCSD energy
    CCSD<Tile, Policy>::evaluate(ccsd_energy.get());

    auto time0 = mpqc::fenced_now(world);

    ExEnv::out0() << indent << "\nEOM-CCSD Excitation Energy \n";
    auto n_roots = ex_energy->n_roots();
    ExEnv::out0() << indent << "Number of roots: " << n_roots << "\n";
    ExEnv::out0() << indent << "Target Precision: " << target_precision << "\n";
    ExEnv::out0() << indent << "Max number of vector per root: " << max_vector_
                  << "\n";
    ExEnv::out0() << indent
                  << "Threshold for norm of new vector: " << vector_threshold_
                  << "\n";
    ExEnv::out0() << indent << "Davidson Solver: " << davidson_solver_ << "\n";
    ExEnv::out0() << indent << "PNO Simulation: "
                  << (eom_pno_.empty() ? "none" : eom_pno_) << "\n";

    if (!eom_pno_.empty()) {
      ExEnv::out0() << indent << "PNO Canonical: "
                    << (eom_pno_canonical_ ? "true" : "false") << "\n";
      ExEnv::out0() << indent << "TcutPNO : " << eom_tpno_ << "\n";
      ExEnv::out0() << indent << "TcutOSV : " << eom_tosv_ << "\n";
    }

    // initialize guest
    ExEnv::out0() << indent << "\nInitialize Guess Vector From CIS \n";

    std::vector<TArray> cis_vector;
    std::vector<numeric_type> cis_eig;
    {
      auto n_guess = ex_energy->n_guess();

      // compute number of guess CIS roots
      if (davidson_solver_ == "multi-state") {
        ex_energy->set_n_roots(n_guess);
      }

      ::mpqc::evaluate(*ex_energy, cis_guess_wfn_);

      ex_energy->set_n_roots(n_roots);

      cis_vector = cis_guess_wfn_->eigen_vector();
      cis_eig = cis_guess_wfn_->eigen_value();
    }

    auto max_iter = this->max_iter_;
    auto result = eom_ccsd_davidson_solver(n_roots, cis_vector, cis_eig,
                                           max_iter, target_precision);

    this->computed_ = true;
    ExcitationEnergy::Provider::set_value(
        ex_energy, std::vector<numeric_type>(result.data(),
                                             result.data() + result.size()));

    auto time1 = mpqc::fenced_now(world);
    ExEnv::out0() << "EOM-CCSD Total Time: "
                  << mpqc::duration_in_s(time0, time1) << " S\n";
  }
}
}  // namespace lcao

}  // namespace mpqc
