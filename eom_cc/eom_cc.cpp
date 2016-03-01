/*
 * eom_cc.cpp
 *
 *  Created on: Feb 26, 2016
 *      Author: jinmei
 */
#include "eom_cc.h"

namespace mpqc{
  void EOM_CCSD::compute_intermediates() {

    TArray tau_ab;
    tau_ab("a,b,i,j") = T2_("a,b,i,j") + (T1_("a,i") * T1_("b,j"));

    TArray gabij_temp, t2_temp;
    gabij_temp("a,b,i,j") = 2.0 * Gabij_("a,b,i,j") - Gabij_("a,b,j,i");
    t2_temp("a,b,i,j") = 2.0 * T2_("a,b,i,j") - T2_("a,b,j,i");

    // \cal{F}
    FIA_("i,a") = //   fia
                  // + t^b_j g^ij_ab
                  gabij_temp("a,b,i,j") * T1_("b,j")
                  ;
    FAB_("a,b") =  //   fab (1 - delta_ab) - fkb t^a_k
                   // + t^d_k g^ak_bd
                   T1_("d,k") * (2.0 * Giabc_("k,a,d,b") - Giabc_("k,a,b,d"))
                   // - 1/2 tau^ad_kl g^kl_bd
                 - tau_ab("a,d,k,l") * gabij_temp("b,d,k,l")
                 ;
    FIJ_("i,j") =  //   fij (1 - delta_ij) + 1/2 fic t^c_j
                   // + t^c_l g^il_jc
                   T1_("c,l") * (2.0 * Gijka_("i,l,j,c") - Gijka_("l,i,j,c"))
                   // + 1/2 tau^cd_jl g^il_cd
                 + tau_ab("c,d,j,l") * gabij_temp("c,d,i,l");

    // \cal{W}mnij
    WKlIj_("k,l,i,j") =  // g^kl_ij
                         Gijkl_("k,l,i,j")

                         // + P(ij) t^c_j g^kl_ic:
                         // + t^c_j g^kl_ic
                       + T1_("c,j") * Gijka_("k,l,i,c")
                         // + t^c_i g^kl_cj
                       + T1_("c,i") * Gijka_("l,k,j,c")

                         // + 1/2 * tau^cd_ij * g^kl_cd
                       + tau_ab("c,d,i,j") * Gabij_("c,d,k,l")
                       ;

    // \cal{W}abef
    WAbCd_("a,b,c,d") =  //  g^ab_cd
                         Gabcd_("a,b,c,d")

                         // - P(ab) t^b_k g^ak_cd:
                         // - t^b_k g^ak_cd
                       - T1_("b,k") * Giabc_("k,a,d,c")
                         // - t^a_k g^kb_cd
                       - T1_("a,k") * Giabc_("k,b,c,d")

                         // + 1/2 tau^ab_kl g^kl_cd
                       + tau_ab("a,b,k,l") * Gabij_("c,d,k,l")
                       ;

    // \cal{W}mbej
    WIbAj_("i,b,a,j") =  // g^ib_aj
                         Gabij_("a,b,i,j")
                         // + t^d_j g^ib_ad
                       + T1_("d,j") * Giabc_("i,b,a,d")
                         // - t^b_l g^il_aj
                       - T1_("b,l") * Gijka_("l,i,j,a")

                         // - (t^db_jl + t^d_j t^b_l) g^il_ad:
                         // + t^bd_jl g^il_ad
                       + T2_("b,d,j,l") * gabij_temp("a,d,i,l")
                                - T2_("d,b,j,l") * Gabij_("a,d,i,l")
                         // - t^d_j t^b_l g^il_ad
                       - T1_("d,j") * (T1_("b,l") * Gabij_("a,d,i,l"))
                       ;

    WIbaJ_("i,b,a,j") =  // g^ib_aj
                       - Giajb_("j,a,i,b") //g_aibj("c,j,b,k")
                         // + t^d_j g^ib_ad
                       - T1_("d,j") * Giabc_("i,b,d,a") // g_abci("c,d,b,k")
                         // - t^b_l g^il_aj
                       + T1_("b,l") * Gijka_("i,l,j,a") // g_aikl("c,j,l,k")

                         // - (t^db_jl + t^d_j t^b_l) g^il_ad:
                         // + t^bd_jl g^il_ad
                       + T2_("d,b,j,l") * Gabij_("d,a,i,l")
                         // - t^d_j t^b_l g^il_ad
                       + T1_("d,j") * (T1_("b,l") * Gabij_("d,a,i,l"))
                       ;

    // \cal{W}abei
    if (WAbCd_.is_initialized()) {
      WAbCi_("a,b,c,i") =  //   g^ab_ci
                           Giabc_("i,c,b,a")

                           // - \cal{F}kc t^ab_ki
                         - FIA_("k,c") * T2_("a,b,k,i")

                           // + t^d_i \cal{W}abcd
                         + T1_("d,i") * WAbCd_("a,b,c,d")

                           // + 1/2 g^kl_ci tau^ab_kl
                         + Gijka_("l,k,i,c") * tau_ab("a,b,k,l")

                           // - P(ab) g^kb_cd t^ad_ki:
                           // - g^kb_cd t^ad_ki
                         - Giabc_("k,b,c,d") * T2_("a,d,k,i")
                           // + g^ka_cd t^bd_ki = + g^ak_cd t^db_ki
                         + (  Giabc_("k,a,d,c") * t2_temp("b,d,i,k")
                            - Giabc_("k,a,c,d") * T2_("b,d,i,k"))

                           // - P(ab) t^a_k (g^kb_ci - t^bd_li g^kl_cd):
                           // - t^a_k (g^kb_ci - t^bd_li g^kl_cd) = - t^a_k (g^kb_ci + t^db_li g^kl_cd)
                         - T1_("a,k")
                           * (  Gabij_("b,c,i,k") // g_aijb("b,k,i,c") ***
                               + (  t2_temp("d,b,l,i") * Gabij_("c,d,k,l")
                                  - T2_("d,b,l,i") * Gabij_("d,c,k,l"))
                              )
                            // + t^b_k (g^ka_ci - t^ad_li g^kl_cd) = + t^b_k (- g^ak_ci + t^ad_li g^lk_cd)
                         + T1_("b,k")
                           * (- Giajb_("k,a,i,c") // g_aibj("a,k,c,i") ***
                              + T2_("a,d,l,i") * Gabij_("c,d,l,k")
                             )
                         ;
    } else {
        throw std::runtime_error(
              "WAbCd_ has not been computed");
    }


    // \cal{W}mbij
    if (WKlIj_.is_initialized()) {
      WKbIj_("k,b,i,j") =  //   g^kb_ij
                           Gijka_("i,j,k,b")

                           // - \cal{F}kc t^bc_ij = + CFkc t^cb_ij
                         + FIA_("k,c") * T2_("c,b,i,j")

                           // - t^b_l \cal{W}klij
                         - T1_("b,l") * WKlIj_("k,l,i,j")

                           // + 0.5 g^kb_cd tau^cd_ij
                         + Giabc_("k,b,c,d") * tau_ab("c,d,i,j")

                           // + P(ij) g^kl_ic t^bc_jl
                           // + g^kl_ic t^bc_jl
                         + Gijka_("k,l,i,c") * t2_temp("b,c,j,l")
                         - Gijka_("l,k,i,c") * T2_("b,c,j,l")
                           // - g^kl_jc t^bc_il = - g^lk_jc t^cb_il
                         - Gijka_("l,k,j,c") * T2_("c,b,i,l")

                           // + P(IJ) t^c_i (g^kb_cj - t^bd_lj g^kl_cd)
                           // + t^c_i (g^kb_cj - t^bd_lj g^kl_cd) = + t^c_i (g^kb_cj + t^bd_jl g^kl_cd)
                         + T1_("c,i")
                           * (  Gabij_("b,c,j,k") // g_aijb("c,j,k,b") ***
                              + t2_temp("b,d,j,l") * Gabij_("c,d,k,l")
                              - T2_("b,d,j,l") * Gabij_("d,c,k,l")
                             )
                            // - t^c_j (g^kb_ci - t^bd_li g^kl_cd) = - t^c_j (- g^kb_ic + t^bd_li g^kl_dc)
                         - T1_("c,j")
                           * (- Giajb_("i,c,k,b") // g_aibj("c,i,b,k") ***
                              + T2_("b,d,l,i") * Gabij_("d,c,k,l")
                             )
                         ;
    } else {
        throw std::runtime_error(
              "WKlIj_ has not been computed");
    }

    // \cal{W}amef
    WAkCd_("a,k,c,d") =  //   g^ak_cd
                         Giabc_("k,a,d,c")
                         // - t^a_l g^lk_cd
                       - T1_("a,l") * Gabij_("c,d,l,k")
                       ;

    // \cal{W}mnie
    WKlIc_("k,l,i,c") =  //   g^kl_ic
                         Gijka_("k,l,i,c")
                         // + t^d_i g^kl_dc
                       + T1_("d,i") * Gabij_("d,c,k,l")
                       ;
    WKliC_("k,l,i,c") =  //   g^kl_ic
                       - Gijka_("l,k,i,c")
                         // + t^d_i g^kl_dc
                       - T1_("d,i") * Gabij_("d,c,l,k")
                       ;
  }

  double EOM_CCSD::compute_energy() {
    // check if intermediates are computed
    compute_intermediates();

    // test intermediates
    if (T1_.get_world().rank() == 0) {
      std::cout << "T1_ " << T1_ << std::endl;
    }

    TArray L1_test, CGac, CGki;
    TArray t2_temp;
    t2_temp("a,b,i,j") = 2.0 * T2_("a,b,i,j") - T2_("b,a,i,j");

    CGac("a,c") =  // - 1/2 t^cd_kl lambda^kl_ad
                 - T2_("c,d,k,l") * t2_temp("a,d,k,l");
    CGki("k,i") =  // 1/2 t^cd_kl lambda^il_cd
                   T2_("c,d,k,l") * t2_temp("c,d,i,l");

    L1_test("a,i") =  FIA_("i,a")
                    + T1_("c,i") * FAB_("c,a")
                    - T1_("a,k") * FIJ_("i,k")
                    + T1_("c,k") * (2.0 * WIbAj_("i,c,a,k") + WIbaJ_("i,c,a,k"))
                    + t2_temp("c,d,i,k") * WAbCi_("c,d,a,k")
                    - t2_temp("a,c,k,l") * WKbIj_("i,c,k,l")
                    - CGac("c,d") * (2.0 * WAkCd_("c,i,d,a") - WAkCd_("c,i,a,d"))
                    - CGki("k,l") * (2.0 * WKlIc_("k,i,l,a") + WKliC_("k,i,l,a"))
                    ;
    if (T1_.get_world().rank() == 0) {
      std::cout << "L1_test " << L1_test << std::endl;
    }

    return 1.0;
  }
}


