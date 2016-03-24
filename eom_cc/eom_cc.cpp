/*
 * eom_cc.cpp
 *
 *  Created on: Feb 26, 2016
 *      Author: jinmei
 */
#include "eom_cc.h"
#include <rapidjson/document.h>
#include <tiledarray.h>

namespace mpqc{

  using TArray = TiledArray::TSpArrayD;

  void EOM_CCSD::compute_FWintermediates() {

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
      WKaIj_("k,a,i,j") =  //   g^ka_ij
                           Gijka_("i,j,k,a")

                           // - \cal{F}kc t^bc_ij = + CFkc t^ca_ij
                         + FIA_("k,c") * T2_("c,a,i,j")

                           // - t^a_l \cal{W}klij
                         - T1_("a,l") * WKlIj_("k,l,i,j")

                           // + 0.5 g^ka_cd tau^cd_ij
                         + Giabc_("k,a,c,d") * tau_ab("c,d,i,j")

                           // + P(ij) g^kl_ic t^ac_jl
                           // + g^kl_ic t^ac_jl
                         + Gijka_("k,l,i,c") * t2_temp("a,c,j,l")
                         - Gijka_("l,k,i,c") * T2_("a,c,j,l")
                           // - g^kl_jc t^ac_il = - g^lk_jc t^ca_il
                         - Gijka_("l,k,j,c") * T2_("c,a,i,l")

                           // + P(ij) t^c_i (g^ka_cj - t^ad_lj g^kl_cd)
                           // + t^c_i (g^ka_cj - t^ad_lj g^kl_cd) = + t^c_i (g^ka_cj + t^ad_jl g^kl_cd)
                         + T1_("c,i")
                           * (  Gabij_("a,c,j,k") // g_aijb("c,j,k,b") ***
                              + t2_temp("a,d,j,l") * Gabij_("c,d,k,l")
                              - T2_("a,d,j,l") * Gabij_("d,c,k,l")
                             )
                            // - t^c_j (g^ka_ci - t^ad_li g^kl_cd) = - t^c_j (- g^ka_ic + t^ad_li g^kl_dc)
                         - T1_("c,j")
                           * (- Giajb_("i,c,k,a") // g_aibj("c,i,b,k") ***
                              + T2_("a,d,l,i") * Gabij_("d,c,k,l")
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
//    WKliC_("k,l,i,c") =  //   g^kl_ic
//                       - Gijka_("l,k,i,c")
//                         // + t^d_i g^kl_dc
//                       - T1_("d,i") * Gabij_("d,c,l,k")
//                       ;
    // WKliC_("k,l,i,c") = - WKliC_("l,k,i,c")
  }

  // compute [HSS C]^A_I
  TArray EOM_CCSD::compute_HSSC(TArray Cai) {
    TArray HSSC;
    HSSC("a,i") =  //   Fae C^e_i
                   FAB_("a,c") * Cai("c,i")
                   // - Fmi C^a_m
                 - FIJ_("k,i") * Cai("a,k")
                   // + Wamie C^e_m
                 + (2.0 * WIbAj_("i,c,a,k") + WIbaJ_("i,c,a,k"))
                   * Cai("c,k")
                 ;
    return HSSC;
  }
  // compute [HSD C]^A_I
  TArray EOM_CCSD::compute_HSDC(TArray Cabij) {
    TArray HSDC;
    HSDC("a,i") =  //   Fme C^ae_im
                   FIA_("k,c")
                   * (2.0 * Cabij("a,c,i,k") - Cabij("c,a,i,k"))
                   // + 1/2 Wamef C^ef_im
                 + WAkCd_("a,k,c,d")
                   * (2.0 * Cabij("c,d,i,k") - Cabij("d,c,i,k"))
                   // - 1/2 Wmnie C^ae_mn
                 - WKlIc_("k,l,i,c")
                   * (2.0 * Cabij("a,c,k,l") - Cabij("d,c,i,k"))
                 ;
    return HSDC;
  }

  // compute [HDS C]^Ab_Ij
  TArray EOM_CCSD::compute_HDSC(TArray Cai) {
    TArray HDSC;
    HDSC("a,b,i,j") =  //   P(ab) Wmaij C^b_m
                       // - WmAjI C^b_m - WMbIj C^A_M
                     - WKaIj_("k,a,j,i") * Cai("b,k")
                     - WKaIj_("k,b,i,j") * Cai("a,k")
                       // + P(ij) Wabej C^e_i
                       // + WAbEj C^E_I + WbAeI C^e_j
                     + WAbCi_("a,b,c,j") * Cai("c,i")
                     + WAbCi_("b,a,c,i") * Cai("c,j")
                       // + P(ab) Wbmfe C^e_m T^af_ij
                       // + WbMfE C^E_M T^Af_Ij + Wbmfe C^e_m T^Af_Ij
                       // - WAMFE C^E_M T^bF_Ij - WAmFe C^e_m T^bF_Ij
                     + (2.0 * WAkCd_("b,k,d,c") - WAkCd_("b,k,c,d") )
                       * Cai("c,k") * T2_("a,d,i,j")
                     + (2.0 * WAkCd_("a,k,d,c") - WAkCd_("a,k,c,d") )
                       * Cai("c,k") * T2_("d,b,i,j")
                       // - P(ij) Wnmje C^e_m t^ab_in
                       // - WnMjE C^E_M T^Ab_In - Wnmje C^e_m T^Af_Ij
                       // + WNMIE C^E_M T^Ab_jN + WNmIe C^e_m T^Ab_jN
                     - (2.0 * WKlIc_("l,k,j,c") - WKlIc_("k,l,j,c") )
                       * Cai("c,k") * T2_("a,b,i,l")
                     - (2.0 * WKlIc_("l,k,i,c") - WKlIc_("k,l,i,c") )
                       * Cai("c,k") * T2_("a,b,l,j")
                     ;
    return HDSC;
  }

  // compute [HDD C]^Ab_Ij
  TArray EOM_CCSD::compute_HDDC(TArray Cabij) {
    TArray GC_ab, GC_ij, HDDC;
    GC_ab("a,b") = Gabij_("a,c,k,l")
                * (2.0 * Cabij("b,c,k,l") - Cabij("c,b,k,l"));
    GC_ij("i,j") = Gabij_("c,d,i,k")
                   * (2.0 * Cabij("c,d,j,k") - Cabij("d,c,j,k"));

    HDDC("a,b,i,j") =  //   P(ab) Fbe C^ae_ij
                       //   Fbe C^ae_ij + Fae C^eb_ij
                       FAB_("b,c") * Cabij("a,c,i,j")
                     + FAB_("a,c") * Cabij("c,b,i,j")
                       // - P(ij) Fmj C^ab_im
                       // - Fmj C^ab_im - Fmi C^ab_mj
                     - FIJ_("k,j") * Cabij("a,b,i,k")
                     + FIJ_("k,i") * Cabij("a,b,k,j")
                       // + 1/2 Wabef C^ef_ij
                     + WAbCd_("a,b,c,d") * Cabij("c,d,i,j")
                       // + 1/2 Wmnij C^ab_mn
                     + WKlIj_("k,l,i,j") * Cabij("a,b,k,l")
                       // + P(ab) P(ij) Wbmje C^ae_im
                       // + Wbmje C^ae_im - Wbmie C^ae_jm
                       // - Wamje C^be_im + Wamie C^be_jm
                     + WIbAj_("k,b,c,j")
                       * (2.0 * Cabij("a,c,i,k") - Cabij("c,a,i,k"))
                     + WIbaJ_("k,b,c,j") * Cabij("a,c,i,k") //
                     + WIbaJ_("k,b,c,i") * Cabij("a,c,k,j") //
                     + WIbaJ_("k,a,c,i") * Cabij("b,c,k,i") //
                     + WIbAj_("k,a,c,i")
                       * (2.0 * Cabij("b,c,j,k") - Cabij("c,b,j,k"))
                     + WIbaJ_("k,a,c,i") * Cabij("b,c,j,k")
                       // - 1/2 P(ab) Wnmfe C^ea_mn t^fb_ij
                       // - 1/2 Wmnfe C^ae_mn t^fb_ij
                       // - 1/2 Wmnef C^eb_mn t^af_ij
                     - GC_ab("d,a") //Gabij_("d,c,k,l")*(2.0*Cabij_("a,c,k,l")-Cabij_("c,a,k,l"))
                       * T2_("d,b,i,j")
                     - GC_ab("d,b") //Gabij_("c,d,k,l")*(2.0*Cabij_("c,b,k,l")-Cabij_("b,c,k,l"))
                       * T2_("a,d,i,j")
                       // + 1/2 P(ij) Wnmfe C^fe_im t^ab_jn
                       // - 1/2 Wnmef C^ef_im t^ab_nj
                       // - 1/2 Wnmef C^ef_jm t^ab_in
                     - GC_ij("l,i") //Gabij_("c,d,l,k")*(2.0*Cabij_("c,d,i,k")-Cabij_("d,c,i,k"))
                       * T2_("a,b,l,j")
                     - GC_ij("l,j") //Gabij_("c,d,l,k")*(2.0*Cabij_("c,d,j,k")-Cabij_("d,c,j,k"))
                       * T2_("a,n,i,l")
                     ;
    return HDDC;

  }

  // read guess vectors from input, eg:
  // "RootNumber": 3 (specify the number of guess vector)
  // "Root1": [1, 0, 2, 0, 1.0] (1st vector, electron excited from 1st occ. to 2nd vir.)
  // "Root2": [1, 1, 1, 2, 1.0] (electrons excited from 1st occ. to 1st & 2nd vir.)
  // "Root3": [1, 0, 2, 0, 2.0, 1, 1, 1, 2, 1.0]
  void EOM_CCSD::read_guess_vectors(rapidjson::Document& in) {

    std::size_t Nroots = in["RootNumber"].GetInt();
    C_.resize(Nroots);

    madness::World& world = T1_.get_world();
    for (std::size_t i = 0; i < Nroots; ++i) {

      const std::string rootn_str = "Root" + std::to_string(i+1);
      const char* RootN = rootn_str.c_str();
      const rapidjson::Value& RN = in[RootN];

      const int num_roots = RN.Size()/5;
      std::vector<int> indices_RN(num_roots);

      if (world.rank() == 0) {
        std::cout << std::endl << "Eigenstate " << i+1
                  << "  Number of Roots: " << num_roots << std::endl;
      }
      for (std::size_t j = 0; j < num_roots; ++j) {
        if (world.rank() == 0) {
          std::cout << "Root " << j+1 << ": ";
        }

        if (RN[j*5] != 0 && RN[j*5+1] != 0) { // double excitations

            const int idx_i = RN[j*5].GetInt() - 1;
            const int idx_j = RN[j*5+1].GetInt() - 1;
            const int idx_a = RN[j*5+2].GetInt() - 1;
            const int idx_b = RN[j*5+3].GetInt() - 1;
            const double coefficient = RN[j*5+4].GetDouble();
            if (world.rank() == 0) {
              std::cout << "Double excitation form i " << idx_i
                        << " j " << idx_i
                        << " to a " << idx_a
                        << " b " << idx_b
                        << " with coefficient " << coefficient << std::endl;
            }

            const std::array<int,4> idx = {{idx_a,idx_b,idx_i,idx_j}};
            const auto tile_idx = T2_.trange().element_to_tile(idx);
            const auto tile_range = T2_.trange().make_tile_range(tile_idx);

            // Construct sparse shape for array
            TiledArray::TensorF tile_norm(T2_.trange().tiles(), 0);
            tile_norm[tile_idx] = tile_range.volume();
            typename TArray::shape_type shape(tile_norm, T2_.trange());

            // Construct the array
            C_[i].Cabij = TArray(world, T2_.trange(), shape);

            // Initialize the tile
            if(C_[i].Cabij.is_local(tile_idx)) {
              // Construct tile
              TiledArray::TensorD tile(tile_range, 0);
              tile[idx] = coefficient;
              C_[i].Cabij.set(tile_idx, tile);
            }
            if (world.rank() == 0) {
              std::cout << "Cabij: " << C_[i].Cabij << std::endl;
            }

        } else if (RN[j*5+1] == 0) { // single excitation

            const int idx_i = RN[j*5].GetInt() - 1;
            const int idx_a = RN[j*5+2].GetInt() - 1;
            const double coefficient = RN[j*5+4].GetDouble();
            if (world.rank() == 0) {
              std::cout << "Single excitation form i " << idx_i
                        << " to a " << idx_a
                        << " with coefficient " << coefficient << std::endl;
            }

            const std::array<int,2> idx = {{idx_a, idx_i}};
            const auto tile_idx = T1_.trange().element_to_tile(idx);
            const auto tile_range = T1_.trange().make_tile_range(tile_idx);

            // Construct sparse shape for array
            TiledArray::TensorF tile_norm(T1_.trange().tiles(), 0);
            tile_norm[tile_idx] = tile_range.volume();
            typename TArray::shape_type shape(tile_norm, T1_.trange());

            // Construct the array
            C_[i].Cai = TArray(world, T1_.trange(), shape);

            // Initialize the tile
            if(C_[i].Cai.is_local(tile_idx)) {
              // Construct tile
              TiledArray::TensorD tile(tile_range, 0);
              tile[idx] = coefficient;
              C_[i].Cai.set(tile_idx, tile);
            }
            if (world.rank() == 0) {
              std::cout << "Cai: " << C_[i].Cai << std::endl;
            }
        }

      } // end of looping over roots in each guess vector

    } // end of looping over guess vectors

  }

  double EOM_CCSD::compute_energy() {
    // check if intermediates are computed
    compute_FWintermediates();

    // test intermediates
    if (T1_.get_world().rank() == 0) {
      std::cout << "T1_ " << T1_ << std::endl;
    }

    TArray L1_test, L2_test, CGac, CGki;
    TArray gabij_temp, t2_temp;
    gabij_temp("a,b,i,j") = 2.0 * Gabij_("a,b,i,j") - Gabij_("a,b,j,i");
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
                    - t2_temp("a,c,k,l") * WKaIj_("i,c,k,l")
                    - CGac("c,d") * (2.0 * WAkCd_("c,i,d,a") - WAkCd_("c,i,a,d"))
                    - CGki("k,l") * (2.0 * WKlIc_("k,i,l,a") - WKlIc_("i,k,l,a"))//+ WKliC_("k,i,l,a"))
                    ;
    if (T1_.get_world().rank() == 0) {
      std::cout << "L1_test " << L1_test << std::endl;
    }

    L2_test("a,b,i,j") =  Gabij_("a,b,i,j")

                        + L1_test("c,i") * WAkCd_("c,j,a,b")
                        + L1_test("c,j") * WAkCd_("c,i,b,a")

                        - L1_test("a,k") * WKlIc_("i,j,k,b")
                        + L1_test("b,k") * - WKlIc_("j,i,k,a")//WKliC_("i,j,k,a")

                        + L1_test("a,i") * FIA_("j,b")
                        + L1_test("b,j") * FIA_("i,a")

                        + T2_("a,c,i,j") * FAB_("c,b")
                        + T2_("c,b,i,j") * FAB_("c,a")

                        - T2_("a,b,i,k") * FIJ_("j,k")
                        - T2_("a,b,k,j") * FIJ_("i,k")

                        + T2_("a,b,k,l") * WKlIj_("i,j,k,l")

                        + T2_("c,d,i,j") * WAbCd_("c,d,a,b")

                        + (  t2_temp("a,c,i,k") * WIbAj_("j,c,b,k")
                           + T2_("a,c,i,k") * WIbaJ_("j,c,b,k"))
                        + T2_("c,b,i,k") * WIbaJ_("j,c,a,k")
                        + T2_("a,c,k,j") * WIbaJ_("i,c,b,k")
                        + (  t2_temp("b,c,j,k") * WIbAj_("i,c,a,k")
                           + T2_("b,c,j,k") * WIbaJ_("i,c,a,k"))

                        + Gabij_("a,c,i,j") * CGac("b,c")
                        + Gabij_("c,b,i,j") * CGac("a,c")

                        - Gabij_("a,b,i,k") * CGki("k,j")
                        - Gabij_("a,b,k,j") * CGki("k,i")
                        ;
    double E_L = dot(gabij_temp("a,b,i,j"), L2_test("a,b,i,j") );
    if (T1_.get_world().rank() == 0) {
      std::cout << "CCSD test E_L: " << E_L << std::endl;

      //std::cout << "Gijkl_: " << Gijkl_ << std::endl;
    }
    return 1.0;
  }
}


