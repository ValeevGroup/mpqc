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

  ExEnv::out0() << indent << "\nInitialize Intermediates in EOM-CCSD\n";

  auto& world = this->wfn_world()->world();
  bool accurate_time = this->lcao_factory().accurate_time();
  auto time0 = mpqc::now(world, accurate_time);

  bool df = this->is_df();
  TArray Tau;
  auto t2 = this->t2();
  auto t1 = this->t1();
  Tau("a,b,i,j") = t2("a,b,i,j") + (t1("a,i") * t1("b,j"));

  std::tie(FAB_, FIA_, FIJ_) = cc::compute_cs_ccsd_F(
      this->lcao_factory(), this->ao_factory(), t1, Tau, df);

  // \cal{W}mnij
  WKlIj_ = cc::compute_cs_ccsd_W_KlIj(this->lcao_factory(), t1, Tau, df);

  // \cal{W}abef
  if (this->method_ != "direct" && this->method_ != "direct_df") {
    WAbCd_ = cc::compute_cs_ccsd_W_AbCd(this->lcao_factory(), t1, Tau, df);
  }

  // \cal{W}mbej
  WIbAj_ = cc::compute_cs_ccsd_W_IbAj(this->lcao_factory(), t1, t2, df);

  WIbaJ_ = cc::compute_cs_ccsd_W_IbaJ(this->lcao_factory(), t1, t2, df);

  // \cal{W}abei
  WAbCi_ = cc::compute_cs_ccsd_W_AbCi(this->lcao_factory(), this->ao_factory(),
                                      t1, t2, Tau, FIA_, WAbCd_, df);

  // \cal{W}mbij
  WKaIj_ = cc::compute_cs_ccsd_W_KaIj(this->lcao_factory(), t1, t2, Tau, FIA_,
                                      WKlIj_, df);

  // \cal{W}amef
  WAkCd_ = cc::compute_cs_ccsd_W_AkCd(this->lcao_factory(), t1, df);

  // \cal{W}mnie
  WKlIc_ = cc::compute_cs_ccsd_W_KlIc(this->lcao_factory(), t1, df);

  auto time1 = mpqc::now(world, accurate_time);
  auto time = mpqc::duration_in_s(time0, time1);

  ExEnv::out0() << indent << "\nTime to Initialize Intermediates in EOM-CCSD: " << time << " S\n";

}

// compute [HSS_HSD C]^A_I
template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> EOM_CCSD<Tile, Policy>::compute_HSS_HSD_C(
    const TArray& Cai, const TArray& Cabij) {
  TArray HSS_HSD_C;

  // HSS * C part
  HSS_HSD_C("a,i") =  //   Fac C^c_i
      FAB_("a,c") * Cai("c,i")
      // - Fki C^a_k
      - FIJ_("k,i") * Cai("a,k")
      // + Wakic C^c_k
      + (2.0 * WIbAj_("k,a,c,i") + WIbaJ_("k,a,c,i")) * Cai("c,k");

  // HSD * C part
  {
    TArray C;
    C("a,c,i,k") = 2.0 * Cabij("a,c,i,k") - Cabij("c,a,i,k");
    HSS_HSD_C("a,i") +=  //   Fkc C^ac_ik
        FIA_("k,c") * C("a,c,i,k")
        // + 1/2 Wakcd C^cd_ik

        + WAkCd_("a,k,c,d") * C("c,d,i,k")
        // - 1/2 Wklic C^ac_kl
        - WKlIc_("k,l,i,c") * C("a,c,k,l");
  }

  return HSS_HSD_C;
}

// compute [HDS_HDD C]^Ab_Ij
template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> EOM_CCSD<Tile, Policy>::compute_HDS_HDD_C(
    const TArray& Cai, const TArray& Cabij) {
  TArray HDS_HDD_C;
  auto t2 = this->t2();
  auto t1 = this->t1();

  // HDS * C part
  {
    HDS_HDD_C("a,b,i,j") =  //   P(ab) Wkaij C^b_k
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
            t2("a,d,i,j") +
        (2.0 * WAkCd_("a,k,d,c") - WAkCd_("a,k,c,d")) * Cai("c,k") *
            t2("d,b,i,j")
        // - P(ij) Wlkjc C^c_k t^ab_il
        // - WlKjC C^C_K T^Ab_Il - Wlkjc C^c_k T^Ad_Ij
        // + WLKIC C^C_K T^Ab_jL + WLkIc C^c_k T^Ab_jL
        -
        (2.0 * WKlIc_("l,k,j,c") - WKlIc_("k,l,j,c")) * Cai("c,k") *
            t2("a,b,i,l") -
        (2.0 * WKlIc_("l,k,i,c") - WKlIc_("k,l,i,c")) * Cai("c,k") *
            t2("a,b,l,j");
  }

  // HDD * C part
  {
    TArray GC_ab, GC_ij;
    GC_ab("a,b") =
        g_ijab_("k,l,a,c") * (2.0 * Cabij("b,c,k,l") - Cabij("c,b,k,l"));
    GC_ij("i,j") =
        g_ijab_("i,k,c,d") * (2.0 * Cabij("c,d,j,k") - Cabij("d,c,j,k"));

    HDS_HDD_C("a,b,i,j") +=  //   P(ab) Fbc C^ac_ij
        //   Fbc C^ac_ij + Fac C^cb_ij
        FAB_("b,c") * Cabij("a,c,i,j") + FAB_("a,c") * Cabij("c,b,i,j")
        // - P(ij) Fkj C^ab_ik
        // - Fkj C^ab_ik - Fki C^ab_kj
        - FIJ_("k,j") * Cabij("a,b,i,k") - FIJ_("k,i") * Cabij("a,b,k,j")
        // + 1/2 Wabcd C^cd_ij
        //        + WAbCd_("a,b,c,d") * Cabij("c,d,i,j")
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

        // Gabij_("d,c,k,l")*(2.0*Cabij_("a,c,k,l")-Cabij_("c,a,k,l"))
        // Gabij_("c,d,k,l")*(2.0*Cabij_("c,b,k,l")-Cabij_("b,c,k,l"))
        - GC_ab("d,a") * t2("d,b,i,j") - GC_ab("d,b") * t2("a,d,i,j")

        // + 1/2 P(ij) Wlkdc C^dc_ik t^ab_jl
        // - 1/2 Wlkcd C^cd_ik t^ab_lj
        // - 1/2 Wlkcd C^cd_jk t^ab_il
        // Gabij_("c,d,l,k")*(2.0*Cabij_("c,d,i,k")-Cabij_("d,c,i,k"))
        // Gabij_("c,d,l,k")*(2.0*Cabij_("c,d,j,k")-Cabij_("d,c,j,k"))
        - GC_ij("l,i") * t2("a,b,l,j") - GC_ij("l,j") * t2("a,b,i,l");

    if (WAbCd_.is_initialized()) {
      HDS_HDD_C("a,b,i,j") += WAbCd_("a,b,c,d") * Cabij("c,d,i,j");
    } else {
      //      auto direct_integral = this->ao_factory().compute_direct(L"(μ ν|
      //      G|κ λ)");
      //      auto Ca =
      //          this->lcao_factory().orbital_registry().retrieve(OrbitalIndex(L"a"));
      //      auto Ci =
      //          this->lcao_factory().orbital_registry().retrieve(OrbitalIndex(L"i"));
      //
      //      TArray tau;
      //      tau("a,b,i,j") = t2("a,b,i,j") + t1("a,i") * t1("b,j");
      //
      //      TArray U;
      //      U("p,q,i,j") = Cabij("r,s,i,j") * direct_integral("p,r,q,s");
      //
      //      TArray tmp;
      //      tmp("a,b,i,j") = -U("p,q,i,j") * Ci("p,k") * Ca("q,a") * t1("b,k")
      //      -
      //          U("p,q,i,j") * Ci("p,k") * Ca("q,b") * t1("a,k") +
      //          U("p,q,i,j") * Ci("p,k") * Ci("q,l") * tau("a,b,k,l") +
      //          U("p,q,i,j") * Ca("p,a") * Ca("q,b");
      //
      //      HDS_HDD_C("a,b,i,j") += tmp("a,b,i,j");

      TArray tau;
      tau("a,b,i,j") = t2("a,b,i,j") + t1("a,i") * t1("b,j");

      // integral direct term
      auto direct_integral = this->get_direct_ao_integral();

      auto Ca =
          this->lcao_factory().orbital_registry().retrieve(OrbitalIndex(L"a"));
      auto Ci =
          this->lcao_factory().orbital_registry().retrieve(OrbitalIndex(L"i"));

      TArray U;
      U("p,r,i,j") = Cabij("c,d,i,j") * Ca("q,c") * Ca("s,d") *
          direct_integral("p,q,r,s");
      U("p,r,i,j") =  0.5*( U("p,r,i,j") + U("r,p,j,i") );
      HDS_HDD_C("a,b,i,j") += U("p,r,i,j") * Ca("p,a") * Ca("r,b")
          - U("r,p,i,j") * Ci("p,k") * Ca("r,a") * t1("b,k")
          - U("p,r,i,j") * Ci("p,k") * Ca("r,b") * t1("a,k")
          + U("p,r,i,j") * Ci("p,k") * Ci("r,l") * tau("a,b,k,l");

    }
  }

  return HDS_HDD_C;
}

template <typename Tile, typename Policy>
EigenVector<typename Tile::numeric_type>
EOM_CCSD<Tile, Policy>::eom_ccsd_davidson_solver(std::size_t max_iter,
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
        HC[i].t1 = compute_HSS_HSD_C(C_[i].t1, C_[i].t2);
        HC[i].t2 = compute_HDS_HDD_C(C_[i].t1, C_[i].t2);

      } else {
        throw ProgrammingError("Guess Vector not initialized", __FILE__,
                               __LINE__);
      }
    }

    auto time1 = mpqc::fenced_now(world);
    EigenVector<double> eig_new = dvd.extrapolate(HC, C_, pred);
    auto time2 = mpqc::fenced_now(world);

    EigenVector<numeric_type> delta_e = eig - eig_new;
    norm_r = delta_e.norm();

    util::print_excitation_energy_iteration(iter, delta_e, eig_new,
                                            mpqc::duration_in_s(time0, time1),
                                            mpqc::duration_in_s(time1, time2));

    eig = eig_new;
    iter++;

  }  // end of while loop

  if (iter == max_iter) {
    throw MaxIterExceeded("Davidson Diagonalization Exceeded Max Iteration",
                          __FILE__, __LINE__, max_iter, "EOM-CCSD");
  }

  ExEnv::out0() << "\n";
  util::print_excitation_energy(eig, false);

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
                       __FILE__, __LINE__, "triplet", "true");
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
      // do not use cis direct method, not efficient
      KeyVal& kv_nonconst = const_cast<KeyVal&>(this->kv_);
      std::string cis_method = (this->df_ ? "df" : "standard");
      kv_nonconst.assign("method", cis_method);
      auto cis = std::make_shared<CIS<Tile, Policy>>(this->kv_);
      ::mpqc::evaluate(*ex_energy, cis);
      guess = cis->eigen_vector();
      kv_nonconst.assign("method", this->method_);
    }

    this->init();

    C_ = std::vector<GuessVector>(n_roots);
    auto t2 = this->t2();
    for (std::size_t i = 0; i < n_roots; i++) {
      C_[i].t1("a,i") = guess[i]("i,a");
      C_[i].t2 = TArray(t2.world(), t2.trange(), t2.shape());
      C_[i].t2.fill(0.0);
    }

    std::size_ max_iter = this->max_iter_;
    auto result = eom_ccsd_davidson_solver(max_iter, target_precision);

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
