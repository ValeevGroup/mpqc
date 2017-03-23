/*
 * eom_cc.cpp
 *
 *  Created on: Feb 26, 2016
 *      Author: jinmei
 */

namespace mpqc {
namespace lcao {

template <typename Tile, typename Policy>
void EOM_CCSD<Tile, Policy>::compute_FWintermediates() {
  bool df = this->is_df();
  TArray Tau;
  auto T2_ = this->t2();
  auto T1_ = this->t1();
  Tau("a,b,i,j") = T2_("a,b,i,j") + (T1_("a,i") * T1_("b,j"));

  std::tie(FAB_, FIA_, FIJ_) = cc::compute_cs_ccsd_F(
      this->lcao_factory(), this->ao_factory(), T1_, Tau, df);

  // \cal{F}
  //   fia + t^b_j g^ij_ab
  //  FIA_("i,a") = Fai_("a,i") + gabij_temp("a,b,i,j") * T1_("b,j");
  //  FAB_("a,b") =  //   fab (1 - delta_ab) - fkb t^a_k
  //      Fab_("a,b") - T1_("a,m")*Fai_("b,m") +
  //      // + t^d_k g^ak_bd
  //      T1_("d,k") * (2.0 * Giabc_("k,a,d,b") - Giabc_("k,a,b,d"))
  //      // - 1/2 tau^ad_kl g^kl_bd
  //      - tau_ab("a,d,k,l") * gabij_temp("b,d,k,l");
  //  FIJ_("i,j") =  //   fij (1 - delta_ij) + 1/2 fic t^c_j
  //      Fij_("i,j") +
  //      // + t^c_l g^il_jc
  //      T1_("c,l") * (2.0 * Gijka_("i,l,j,c") - Gijka_("l,i,j,c"))
  //      // + 1/2 tau^cd_jl g^il_cd
  //      + tau_ab("c,d,j,l") * gabij_temp("c,d,i,l");

  // \cal{W}mnij
  WKlIj_ = cc::compute_cs_ccsd_W_KlIj(this->lcao_factory(), T1_, Tau, df);
  //    WKlIj_("k,l,i,j") =  // g^kl_ij
  //        Gijkl_("k,l,i,j")
  //
  //        // + P(ij) t^c_j g^kl_ic:
  //        // + t^c_j g^kl_ic
  //        + T1_("c,j") * Gijka_("k,l,i,c")
  //        // + t^c_i g^kl_cj
  //        + T1_("c,i") * Gijka_("l,k,j,c")
  //
  //        // + 1/2 * tau^cd_ij * g^kl_cd
  //        + Tau("c,d,i,j") * Gabij_("c,d,k,l");

  // \cal{W}abef
  WAbCd_ = cc::compute_cs_ccsd_W_AbCd(this->lcao_factory(), T1_, Tau, df);
  //  WAbCd_("a,b,c,d") =  //  g^ab_cd
  //      Gabcd_("a,b,c,d")
  //
  //      // - P(ab) t^b_k g^ak_cd:
  //      // - t^b_k g^ak_cd
  //      - T1_("b,k") * Giabc_("k,a,d,c")
  //      // - t^a_k g^kb_cd
  //      - T1_("a,k") * Giabc_("k,b,c,d")
  //
  //      // + 1/2 tau^ab_kl g^kl_cd
  //      + Tau("a,b,k,l") * Gabij_("c,d,k,l");

  // \cal{W}mbej
  WIbAj_ = cc::compute_cs_ccsd_W_IbAj(this->lcao_factory(), T1_, T2_, df);
  //  WIbAj_("i,b,a,j") =  // g^ib_aj
  //      Gabij_("a,b,i,j")
  //      // + t^d_j g^ib_ad
  //      + T1_("d,j") * Giabc_("i,b,a,d")
  //      // - t^b_l g^il_aj
  //      - T1_("b,l") * Gijka_("l,i,j,a")
  //
  //      // - (t^db_jl + t^d_j t^b_l) g^il_ad:
  //      // + t^bd_jl g^il_ad
  //      + T2_("b,d,j,l") * gabij_temp("a,d,i,l") -
  //      T2_("d,b,j,l") * Gabij_("a,d,i,l")
  //      // - t^d_j t^b_l g^il_ad
  //      - T1_("d,j") * (T1_("b,l") * Gabij_("a,d,i,l"));

  WIbaJ_ = cc::compute_cs_ccsd_W_IbaJ(this->lcao_factory(), T1_, T2_, df);
  //  WIbaJ_("i,b,a,j") =     // g^ib_aj
  //      -Giajb_("j,a,i,b")  // g_aibj("c,j,b,k")
  //      // + t^d_j g^ib_ad
  //      - T1_("d,j") * Giabc_("i,b,d,a")  // g_abci("c,d,b,k")
  //      // - t^b_l g^il_aj
  //      + T1_("b,l") * Gijka_("i,l,j,a")  // g_aikl("c,j,l,k")
  //
  //      // - (t^db_jl + t^d_j t^b_l) g^il_ad:
  //      // + t^bd_jl g^il_ad
  //      + T2_("d,b,j,l") * Gabij_("d,a,i,l")
  //      // - t^d_j t^b_l g^il_ad
  //      + T1_("d,j") * (T1_("b,l") * Gabij_("d,a,i,l"));

  // \cal{W}abei
  WAbCi_ = cc::compute_cs_ccsd_W_AbCi(this->lcao_factory(), T1_, T2_, Tau, FIA_,
                                      WAbCd_, df);
  //  if (WAbCd_.is_initialized()) {
  //    WAbCi_("a,b,c,i") =  //   g^ab_ci
  //        Giabc_("i,c,b,a")
  //
  //        // - \cal{F}kc t^ab_ki
  //        - FIA_("k,c") * T2_("a,b,k,i")
  //
  //        // + t^d_i \cal{W}abcd
  //        + T1_("d,i") * WAbCd_("a,b,c,d")
  //
  //        // + 1/2 g^kl_ci tau^ab_kl
  //        + Gijka_("l,k,i,c") * Tau("a,b,k,l")
  //
  //        // - P(ab) g^kb_cd t^ad_ki:
  //        // - g^kb_cd t^ad_ki
  //        - Giabc_("k,b,c,d") * T2_("a,d,k,i")
  //        // + g^ka_cd t^bd_ki = + g^ak_cd t^db_ki
  //        + (Giabc_("k,a,d,c") * t2_temp("b,d,i,k") -
  //           Giabc_("k,a,c,d") * T2_("b,d,i,k"))
  //
  //        // - P(ab) t^a_k (g^kb_ci - t^bd_li g^kl_cd):
  //        // - t^a_k (g^kb_ci - t^bd_li g^kl_cd) = - t^a_k (g^kb_ci + t^db_li
  //        // g^kl_cd)
  //        -
  //        T1_("a,k") * (Gabij_("b,c,i,k")  // g_aijb("b,k,i,c") ***
  //                      + (t2_temp("d,b,l,i") * Gabij_("c,d,k,l") -
  //                         T2_("d,b,l,i") * Gabij_("d,c,k,l")))
  //        // + t^b_k (g^ka_ci - t^ad_li g^kl_cd) = + t^b_k (- g^ak_ci +
  //        t^ad_li
  //        // g^lk_cd)
  //        +
  //        T1_("b,k") * (-Giajb_("k,a,i,c")  // g_aibj("a,k,c,i") ***
  //                      + T2_("a,d,l,i") * Gabij_("c,d,l,k"));
  //  } else {
  //    throw std::runtime_error("WAbCd_ has not been computed");
  //  }

  // \cal{W}mbij
  WKaIj_ = cc::compute_cs_ccsd_W_KaIj(this->lcao_factory(), T1_, T2_, Tau, FIA_,
                                      WKlIj_, df);
  //  if (WKlIj_.is_initialized()) {
  //    WKaIj_("k,a,i,j") =  //   g^ka_ij
  //        Gijka_("i,j,k,a")
  //
  //        // - \cal{F}kc t^bc_ij = + CFkc t^ca_ij
  //        + FIA_("k,c") * T2_("c,a,i,j")
  //
  //        // - t^a_l \cal{W}klij
  //        - T1_("a,l") * WKlIj_("k,l,i,j")
  //
  //        // + 0.5 g^ka_cd tau^cd_ij
  //        + Giabc_("k,a,c,d") * Tau("c,d,i,j")
  //
  //        // + P(ij) g^kl_ic t^ac_jl
  //        // + g^kl_ic t^ac_jl
  //        + Gijka_("k,l,i,c") * t2_temp("a,c,j,l") -
  //        Gijka_("l,k,i,c") * T2_("a,c,j,l")
  //        // - g^kl_jc t^ac_il = - g^lk_jc t^ca_il
  //        - Gijka_("l,k,j,c") * T2_("c,a,i,l")
  //
  //        // + P(ij) t^c_i (g^ka_cj - t^ad_lj g^kl_cd)
  //        // + t^c_i (g^ka_cj - t^ad_lj g^kl_cd) = + t^c_i (g^ka_cj + t^ad_jl
  //        // g^kl_cd)
  //        +
  //        T1_("c,i") * (Gabij_("a,c,j,k")  // g_aijb("c,j,k,b") ***
  //                      + t2_temp("a,d,j,l") * Gabij_("c,d,k,l") -
  //                      T2_("a,d,j,l") * Gabij_("d,c,k,l"))
  //        // - t^c_j (g^ka_ci - t^ad_li g^kl_cd) = - t^c_j (- g^ka_ic +
  //        t^ad_li
  //        // g^kl_dc)
  //        -
  //        T1_("c,j") * (-Giajb_("i,c,k,a")  // g_aibj("c,i,b,k") ***
  //                      + T2_("a,d,l,i") * Gabij_("d,c,k,l"));
  //  } else {
  //    throw std::runtime_error("WKlIj_ has not been computed");
  //  }

  // \cal{W}amef
  WAkCd_ = cc::compute_cs_ccsd_W_AkCd(this->lcao_factory(), T1_, df);
  //  WAkCd_("a,k,c,d") =  //   g^ak_cd
  //      Giabc_("k,a,d,c")
  //      // - t^a_l g^lk_cd
  //      - T1_("a,l") * Gabij_("c,d,l,k");

  // \cal{W}mnie
  WKlIc_ = cc::compute_cs_ccsd_W_KlIc(this->lcao_factory(), T1_, df);
  //  WKlIc_("k,l,i,c") =  //   g^kl_ic
  //      Gijka_("k,l,i,c")
  //      // + t^d_i g^kl_dc
  //      + T1_("d,i") * Gabij_("d,c,k,l");
  //    WKliC_("k,l,i,c") =  //   g^kl_ic
  //                       - Gijka_("l,k,i,c")
  //                         // + t^d_i g^kl_dc
  //                       - T1_("d,i") * Gabij_("d,c,l,k")
  //                       ;
  // WKliC_("k,l,i,c") = - WKliC_("l,k,i,c")
}

// compute [HSS C]^A_I
template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> EOM_CCSD<Tile, Policy>::compute_HSSC(TArray Cai) {
  TArray HSSC;
  HSSC("a,i") =  //   Fac C^c_i
      FAB_("a,c") * Cai("c,i")
      // - Fki C^a_k
      - FIJ_("k,i") * Cai("a,k")
      // + Wakic C^c_k
      + (2.0 * WIbAj_("k,a,c,i") + WIbaJ_("k,a,c,i")) * Cai("c,k");
  return HSSC;
}

// compute [HSD C]^A_I
template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> EOM_CCSD<Tile, Policy>::compute_HSDC(TArray Cabij) {
  TArray HSDC;
  HSDC("a,i") =  //   Fkc C^ac_ik
      FIA_("k,c") * (2.0 * Cabij("a,c,i,k") - Cabij("c,a,i,k"))
      // + 1/2 Wakcd C^cd_ik
      + WAkCd_("a,k,c,d") * (2.0 * Cabij("c,d,i,k") - Cabij("d,c,i,k"))
      // - 1/2 Wklic C^ac_kl
      - WKlIc_("k,l,i,c") * (2.0 * Cabij("a,c,k,l") - Cabij("a,c,l,k"));
  return HSDC;
}

// compute [HDS C]^Ab_Ij
template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> EOM_CCSD<Tile, Policy>::compute_HDSC(TArray Cai) {
  TArray HDSC;
  auto T2_ = this->t2();
  auto T1_ = this->t1();
  HDSC("a,b,i,j") =  //   P(ab) Wkaij C^b_k
      // - WkAjI C^b_k - WKbIj C^A_K
      -WKaIj_("k,a,j,i") * Cai("b,k") - WKaIj_("k,b,i,j") * Cai("a,k")
      // + P(ij) Wabcj C^c_i
      // + WAbCj C^C_I + WbAcI C^c_j
      + WAbCi_("a,b,c,j") * Cai("c,i") + WAbCi_("b,a,c,i") * Cai("c,j")
      // + P(ab) Wbkdc C^c_k T^ad_ij
      // + WbKdC C^C_K T^Ad_Ij + Wbkdc C^c_k T^Ad_Ij
      // - WAKDC C^C_K T^bD_Ij - WAkDc C^c_k T^bD_Ij
      +
      (2.0 * WAkCd_("b,k,d,c") - WAkCd_("b,k,c,d")) * Cai("c,k") *
          T2_("a,d,i,j") +
      (2.0 * WAkCd_("a,k,d,c") - WAkCd_("a,k,c,d")) * Cai("c,k") *
          T2_("d,b,i,j")
      // - P(ij) Wlkjc C^c_k t^ab_il
      // - WlKjC C^C_K T^Ab_Il - Wlkjc C^c_k T^Ad_Ij
      // + WLKIC C^C_K T^Ab_jL + WLkIc C^c_k T^Ab_jL
      -
      (2.0 * WKlIc_("l,k,j,c") - WKlIc_("k,l,j,c")) * Cai("c,k") *
          T2_("a,b,i,l") -
      (2.0 * WKlIc_("l,k,i,c") - WKlIc_("k,l,i,c")) * Cai("c,k") *
          T2_("a,b,l,j");
  return HDSC;
}

// compute [HDD C]^Ab_Ij
template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> EOM_CCSD<Tile, Policy>::compute_HDDC(TArray Cabij) {
  TArray GC_ab, GC_ij, HDDC;
  auto T2_ = this->t2();
  auto T1_ = this->t1();
  GC_ab("a,b") =
      g_ijab_("k,l,a,c") * (2.0 * Cabij("b,c,k,l") - Cabij("c,b,k,l"));
  GC_ij("i,j") =
      g_ijab_("i,k,c,d") * (2.0 * Cabij("c,d,j,k") - Cabij("d,c,j,k"));

  HDDC("a,b,i,j") =  //   P(ab) Fbc C^ac_ij
      //   Fbc C^ac_ij + Fac C^cb_ij
      FAB_("b,c") * Cabij("a,c,i,j") + FAB_("a,c") * Cabij("c,b,i,j")
      // - P(ij) Fkj C^ab_ik
      // - Fkj C^ab_ik - Fki C^ab_kj
      - FIJ_("k,j") * Cabij("a,b,i,k") - FIJ_("k,i") * Cabij("a,b,k,j")
      // + 1/2 Wabcd C^cd_ij
      + WAbCd_("a,b,c,d") * Cabij("c,d,i,j")
      // + 1/2 Wklij C^ab_kl
      + WKlIj_("k,l,i,j") * Cabij("a,b,k,l")
      // + P(ab) P(ij) Wbkjc C^ac_ik
      // + Wbkjc C^ac_ik - Wbkic C^ac_jk
      // - Wakjc C^bc_ik + Wakic C^bc_jk
      + WIbAj_("k,b,c,j") * (2.0 * Cabij("a,c,i,k") - Cabij("c,a,i,k")) +
      WIbaJ_("k,b,c,j") * Cabij("a,c,i,k")
      //
      + WIbaJ_("k,b,c,i") * Cabij("a,c,k,j")
      //
      + WIbaJ_("k,a,c,j") * Cabij("b,c,k,i")
      //
      + WIbAj_("k,a,c,i") * (2.0 * Cabij("b,c,j,k") - Cabij("c,b,j,k")) +
      WIbaJ_("k,a,c,i") * Cabij("b,c,j,k")
      // - 1/2 P(ab) g^lk_dc C^ca_kl t^db_ij
      // - 1/2 g^kl_dc C^ac_kl t^db_ij
      // - 1/2 g^kl_cd C^cb_kl t^ad_ij
      -
      GC_ab(
          "d,a")  // Gabij_("d,c,k,l")*(2.0*Cabij_("a,c,k,l")-Cabij_("c,a,k,l"))
          * T2_("d,b,i,j") -
      GC_ab(
          "d,b")  // Gabij_("c,d,k,l")*(2.0*Cabij_("c,b,k,l")-Cabij_("b,c,k,l"))
          * T2_("a,d,i,j")
      // + 1/2 P(ij) Wlkdc C^dc_ik t^ab_jl
      // - 1/2 Wlkcd C^cd_ik t^ab_lj
      // - 1/2 Wlkcd C^cd_jk t^ab_il
      -
      GC_ij(
          "l,i")  // Gabij_("c,d,l,k")*(2.0*Cabij_("c,d,i,k")-Cabij_("d,c,i,k"))
          * T2_("a,b,l,j") -
      GC_ij(
          "l,j")  // Gabij_("c,d,l,k")*(2.0*Cabij_("c,d,j,k")-Cabij_("d,c,j,k"))
          * T2_("a,b,i,l");
  return HDDC;
}

template <typename Tile, typename Policy>
EigenVector<typename Tile::numeric_type>
EOM_CCSD<Tile, Policy>::davidson_solver(std::size_t max_iter,
                                        double convergence) {
  madness::World& world =
      C_[0].t1.is_initialized() ? C_[0].t1.world() : C_[0].t2.world();
  std::size_t iter = 0;
  std::size_t n_roots = C_.size();
  double norm_r = 1.0;

  /// make preconditioner
  Preconditioner pred;
  {
    EigenVector<numeric_type> eps_o =
        array_ops::array_to_eigen(FIJ_).diagonal();
    EigenVector<numeric_type> eps_v =
        array_ops::array_to_eigen(FAB_).diagonal();

    pred = Preconditioner(eps_o, eps_v);
  }

  /// make davidson object
  DavidsonDiag<GuessVector> dvd(n_roots, false);

  EigenVector<double> eig = EigenVector<double>::Zero(n_roots);

  while (iter < max_iter && norm_r > convergence) {
    auto time0 = mpqc::fenced_now(world);
    std::size_t dim = C_.size();
    //    ExEnv::out0() << "vector dimension: " << dim << std::endl;

    // compute product of H with guess vector
    std::vector<GuessVector> HC(dim);
    for (std::size_t i = 0; i < dim; ++i) {
      if (C_[i].t1.is_initialized() && C_[i].t2.is_initialized()) {
        TArray HSSC_ai = compute_HSSC(C_[i].t1);
        TArray HDSC_abij = compute_HDSC(C_[i].t1);
        TArray HSDC_ai = compute_HSDC(C_[i].t2);
        TArray HDDC_abij = compute_HDDC(C_[i].t2);

        HC[i].t1("a,i") = HSSC_ai("a,i") + HSDC_ai("a,i");
        HC[i].t2("a,b,i,j") = HDSC_abij("a,b,i,j") + HDDC_abij("a,b,i,j");

      } else {
        throw ProgrammingError("Guess Vector not initialized", __FILE__,
                               __LINE__);
      }
    }

    auto time1 = mpqc::fenced_now(world);
    EigenVector<double> eig_new = dvd.extrapolate(HC, C_, pred);
    auto time2 = mpqc::fenced_now(world);

    norm_r = (eig - eig_new).norm();

    detail::print_excitation_energy_iteration(
        iter, norm_r, eig_new, mpqc::duration_in_s(time0, time1),
        mpqc::duration_in_s(time1, time2));

    eig = eig_new;
    iter++;

  }  // end of while loop

  if (iter == max_iter) {
    throw MaxIterExceeded("Davidson Diagonalization Exceeded Max Iteration",
                          __FILE__, __LINE__, max_iter, "EOM-CCSD");
  }

  ExEnv::out0() << "\n";
  detail::print_excitation_energy(eig, false);

  return eig;
}

template <typename Tile, typename Policy>
void EOM_CCSD<Tile, Policy>::evaluate(ExcitationEnergy* ex_energy) {
  auto target_precision = ex_energy->target_precision(0);
  if (!this->computed()) {
    auto& world = this->wfn_world()->world();

    bool triplets = ex_energy->triplets();

    if (triplets) {
      throw InputError("EOM-CCSD only supports Singlets at this moment!\n",
                       __FILE__, __LINE__, "triplet","true");
    }

    auto ccsd_energy =
        std::make_shared<Energy>(this->shared_from_this(), target_precision);
    // do CCSD energy
    CCSD<Tile, Policy>::evaluate(ccsd_energy.get());

    auto time0 = mpqc::fenced_now(world);

    ExEnv::out0() << indent << "\nEOM-CCSD Excitation Energy \n";
    auto n_roots = ex_energy->n_roots();

    // initialize guest
    ExEnv::out0() << indent << "\nInitialize Guess Vector From CIS \n";

    std::vector<TArray> guess;
    {
      auto cis = std::make_shared<CIS<Tile, Policy>>(this->kv_);
      ::mpqc::evaluate(*ex_energy, cis);
      guess = cis->eigen_vector();
    }

    ExEnv::out0() << indent << "\nInitialize Intermediates in EOM-CCSD\n";
    this->init();

    C_ = std::vector<GuessVector>(n_roots);
    auto t2 = this->t2();
    for (std::size_t i = 0; i < n_roots; i++) {
      C_[i].t1("a,i") = guess[i]("i,a");
      C_[i].t2 = TArray(t2.world(), t2.trange(), t2.shape());
      C_[i].t2.fill(0.0);
    }

    auto result = davidson_solver(30, target_precision);

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
