/*
 * eom_cc.cpp
 *
 *  Created on: Feb 26, 2016
 *      Author: jinmei
 */
#include "mpqc/chemistry/qc/lcao/cc/eom/eom_ccsd.h"

#include <math.h>
#include <rapidjson/document.h>
#include <tiledarray.h>
#include <Eigen/Eigenvalues>

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
    HSSC("a,i") =  //   Fac C^c_i
                   FAB_("a,c") * Cai("c,i")
                   // - Fki C^a_k
                 - FIJ_("k,i") * Cai("a,k")
                   // + Wakic C^c_k
                 + (2.0 * WIbAj_("k,a,c,i") + WIbaJ_("k,a,c,i"))
                   * Cai("c,k")
                 ;
    return HSSC;
  }

  // compute [HSD C]^A_I
  TArray EOM_CCSD::compute_HSDC(TArray Cabij) {
    TArray HSDC;
    HSDC("a,i") =  //   Fkc C^ac_ik
                   FIA_("k,c")
                   * (2.0 * Cabij("a,c,i,k") - Cabij("c,a,i,k"))
                   // + 1/2 Wakcd C^cd_ik
                 + WAkCd_("a,k,c,d")
                   * (2.0 * Cabij("c,d,i,k") - Cabij("d,c,i,k"))
                   // - 1/2 Wklic C^ac_kl
                 - WKlIc_("k,l,i,c")
                   * (2.0 * Cabij("a,c,k,l") - Cabij("a,c,l,k"))
                 ;
    return HSDC;
  }

  // compute [HDS C]^Ab_Ij
  TArray EOM_CCSD::compute_HDSC(TArray Cai) {
    TArray HDSC;
    HDSC("a,b,i,j") =  //   P(ab) Wkaij C^b_k
                       // - WkAjI C^b_k - WKbIj C^A_K
                     - WKaIj_("k,a,j,i") * Cai("b,k")
                     - WKaIj_("k,b,i,j") * Cai("a,k")
                       // + P(ij) Wabcj C^c_i
                       // + WAbCj C^C_I + WbAcI C^c_j
                     + WAbCi_("a,b,c,j") * Cai("c,i")
                     + WAbCi_("b,a,c,i") * Cai("c,j")
                       // + P(ab) Wbkdc C^c_k T^ad_ij
                       // + WbKdC C^C_K T^Ad_Ij + Wbkdc C^c_k T^Ad_Ij
                       // - WAKDC C^C_K T^bD_Ij - WAkDc C^c_k T^bD_Ij
                     + (2.0 * WAkCd_("b,k,d,c") - WAkCd_("b,k,c,d") )
                       * Cai("c,k") * T2_("a,d,i,j")
                     + (2.0 * WAkCd_("a,k,d,c") - WAkCd_("a,k,c,d") )
                       * Cai("c,k") * T2_("d,b,i,j")
                       // - P(ij) Wlkjc C^c_k t^ab_il
                       // - WlKjC C^C_K T^Ab_Il - Wlkjc C^c_k T^Ad_Ij
                       // + WLKIC C^C_K T^Ab_jL + WLkIc C^c_k T^Ab_jL
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

    HDDC("a,b,i,j") =  //   P(ab) Fbc C^ac_ij
                       //   Fbc C^ac_ij + Fac C^cb_ij
                       FAB_("b,c") * Cabij("a,c,i,j")
                     + FAB_("a,c") * Cabij("c,b,i,j")
                       // - P(ij) Fkj C^ab_ik
                       // - Fkj C^ab_ik - Fki C^ab_kj
                     - FIJ_("k,j") * Cabij("a,b,i,k")
                     - FIJ_("k,i") * Cabij("a,b,k,j")
                       // + 1/2 Wabcd C^cd_ij
                     + WAbCd_("a,b,c,d") * Cabij("c,d,i,j")
                       // + 1/2 Wklij C^ab_kl
                     + WKlIj_("k,l,i,j") * Cabij("a,b,k,l")
                       // + P(ab) P(ij) Wbkjc C^ac_ik
                       // + Wbkjc C^ac_ik - Wbkic C^ac_jk
                       // - Wakjc C^bc_ik + Wakic C^bc_jk
                     + WIbAj_("k,b,c,j")
                       * (2.0 * Cabij("a,c,i,k") - Cabij("c,a,i,k"))
                     + WIbaJ_("k,b,c,j") * Cabij("a,c,i,k")
                       //
                     + WIbaJ_("k,b,c,i") * Cabij("a,c,k,j")
                       //
                     + WIbaJ_("k,a,c,j") * Cabij("b,c,k,i")
                       //
                     + WIbAj_("k,a,c,i")
                       * (2.0 * Cabij("b,c,j,k") - Cabij("c,b,j,k"))
                     + WIbaJ_("k,a,c,i") * Cabij("b,c,j,k")
                       // - 1/2 P(ab) g^lk_dc C^ca_kl t^db_ij
                       // - 1/2 g^kl_dc C^ac_kl t^db_ij
                       // - 1/2 g^kl_cd C^cb_kl t^ad_ij
                     - GC_ab("d,a") //Gabij_("d,c,k,l")*(2.0*Cabij_("a,c,k,l")-Cabij_("c,a,k,l"))
                       * T2_("d,b,i,j")
                     - GC_ab("d,b") //Gabij_("c,d,k,l")*(2.0*Cabij_("c,b,k,l")-Cabij_("b,c,k,l"))
                       * T2_("a,d,i,j")
                       // + 1/2 P(ij) Wlkdc C^dc_ik t^ab_jl
                       // - 1/2 Wlkcd C^cd_ik t^ab_lj
                       // - 1/2 Wlkcd C^cd_jk t^ab_il
                     - GC_ij("l,i") //Gabij_("c,d,l,k")*(2.0*Cabij_("c,d,i,k")-Cabij_("d,c,i,k"))
                       * T2_("a,b,l,j")
                     - GC_ij("l,j") //Gabij_("c,d,l,k")*(2.0*Cabij_("c,d,j,k")-Cabij_("d,c,j,k"))
                       * T2_("a,b,i,l")
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

    madness::World& world = T1_.world();
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
            TiledArray::TensorF tile_norm(T2_.trange().tiles_range(), 0);
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
            TiledArray::TensorF tile_norm(T1_.trange().tiles_range(), 0);
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

  void EOM_CCSD::davidson_solver(std::size_t max_iter, double convergence) {

    madness::World& world = C_[0].Cai.is_initialized()? C_[0].Cai.world()
                                                      : C_[0].Cabij.world();
    std::size_t iter = 0;
    double norm_r = 1.0;

    std::size_t dim;
    while (iter < max_iter && norm_r > convergence) {
      dim = C_.size();
      Eigen::MatrixXd CHC(dim, dim);
      if (world.rank() == 0)
         std::cout << "CHC dimension: " << dim << std::endl;

      std::vector<guess_vector> HC(dim);
      for (std::size_t i = 0; i < dim; ++i) {

        if (C_[i].Cai.is_initialized()
            && C_[i].Cabij.is_initialized()) {
          TArray HSSC_ai = compute_HSSC(C_[i].Cai);
          TArray HDSC_abij = compute_HDSC(C_[i].Cai);
          TArray HSDC_ai = compute_HSDC(C_[i].Cabij);
          TArray HDDC_abij = compute_HDDC(C_[i].Cabij);

          HC[i].Cai("a,i") = HSSC_ai("a,i") + HSDC_ai("a,i");
          HC[i].Cabij("a,b,i,j") = HDSC_abij("a,b,i,j") + HDDC_abij("a,b,i,j");

        } else if (C_[i].Cai.is_initialized()) {
            HC[i].Cai = compute_HSSC(C_[i].Cai);
            HC[i].Cabij = compute_HDSC(C_[i].Cai);

        } else if (C_[i].Cabij.is_initialized()) {
            HC[i].Cai = compute_HSDC(C_[i].Cabij);
            HC[i].Cabij = compute_HDDC(C_[i].Cabij);
        }

        double CHC_ai = 0.0;
        double CHC_abij = 0.0;
        for (std::size_t j = 0; j < dim; ++j) {
          if (C_[j].Cai.is_initialized())
            CHC_ai =  C_[j].Cai("a,i") * HC[i].Cai("a,i");

          if (C_[j].Cabij.is_initialized())
            CHC_abij = C_[j].Cabij("a,b,i,j") * HC[i].Cabij("a,b,i,j");

          CHC(i,j) = CHC_ai + CHC_abij;
        }
      }

      // compute CHC a^k = \lambda^k a^k
      Eigen::EigenSolver<Eigen::MatrixXd> es(CHC);
      Eigen::VectorXd es_values = es.eigenvalues().real();
      Eigen::MatrixXd es_vectors = es.eigenvectors().real();

      if (world.rank() == 0) {
        std::cout << "CHC:" << std::endl
                  << CHC << std::endl;
        std::cout << "The eigenvalues of CHC are:" << std::endl
                  << es.eigenvalues() << std::endl;
        std::cout << "The matrix of eigenvectors, V, is:" << std::endl
                  << es.eigenvectors() << std::endl;
        std::cout << "real eigenvalues is:" << std::endl
                  << es_values << std::endl;
        std::cout << "real eigenvectors is:" << std::endl
                   << es_vectors << std::endl;
      }

      // compute residual r^k = a^k_i (HC_i - \lambda^k C_i)
      std::vector<guess_vector> r(dim);
      for (std::size_t k = 0; k < dim; ++k) {
        const double e_k = es_values(k);

        for (std::size_t i = 0; i < dim; ++i) {
          if (C_[i].Cai.is_initialized()) {
              if (!r[k].Cai.is_initialized()) {
                TiledArray::TensorF tile_norms(C_[i].Cai.trange().tiles_range(), 0.0f);
                TArray::shape_type shape(tile_norms, C_[i].Cai.trange());
                r[k].Cai = TArray(C_[i].Cai.world(), C_[i].Cai.trange(),
                                  shape, C_[i].Cai.pmap());
              }
              r[k].Cai("a,i") += es_vectors(i,k)
                                * (HC[i].Cai("a,i") - C_[i].Cai("a,i") * e_k);
          }
          if (C_[i].Cabij.is_initialized()) {
              if (!r[k].Cabij.is_initialized()) {
                TiledArray::TensorF tile_norms(C_[i].Cabij.trange().tiles_range(), 0.0f);
                TArray::shape_type shape(tile_norms, C_[i].Cabij.trange());
                r[k].Cabij = TArray(C_[i].Cabij.world(), C_[i].Cabij.trange(),
                                  shape, C_[i].Cabij.pmap());
              }
              r[k].Cabij("a,b,i,j") += es_vectors(i,k)
                                      * (HC[i].Cabij("a,b,i,j") - C_[i].Cabij("a,b,i,j") * e_k);
          }
        }
      }

      // compute norm of r
      double norm_square_r = 0.0;
      for (std::size_t i = 0; i < dim; ++i) {
        if (r[i].Cai.is_initialized())
          norm_square_r += r[i].Cai("a,i").squared_norm();

        if (r[i].Cabij.is_initialized())
          norm_square_r += r[i].Cabij("a,b,i,j").squared_norm();
      }
      norm_r = sqrt(norm_square_r);

      // compute preconditioned residual vector:
      // \gamma^k = r^k / (D - \lambda^k I)
      std::vector<guess_vector> gamma(dim);

      // Schmidt orthonormalize \gamma^k against C^k

      // append \gamma^k to C
      C_.resize(dim*2);
      for (std::size_t i = 0; i < dim; ++i) {
        if (gamma[i].Cai.is_initialized())
          C_[i+dim].Cai = gamma[i].Cai;

        if (gamma[i].Cabij.is_initialized())
          C_[i+dim].Cabij = gamma[i].Cabij;

      }
   } // end of while loop

  }

  double EOM_CCSD::compute_energy(std::size_t max_iter, double convergence) {
    // check if intermediates are computed
    compute_FWintermediates();

    // test intermediates
//    if (T1_.world().rank() == 0) {
//      std::cout << "T1_ " << T1_ << std::endl;
//    }
//
//    TArray L1_test, L2_test, CGac, CGki;
//    TArray gabij_temp, t2_temp;
//    gabij_temp("a,b,i,j") = 2.0 * Gabij_("a,b,i,j") - Gabij_("a,b,j,i");
//    t2_temp("a,b,i,j") = 2.0 * T2_("a,b,i,j") - T2_("b,a,i,j");
//
//    CGac("a,c") =  // - 1/2 t^cd_kl lambda^kl_ad
//                 - T2_("c,d,k,l") * t2_temp("a,d,k,l");
//    CGki("k,i") =  // 1/2 t^cd_kl lambda^il_cd
//                   T2_("c,d,k,l") * t2_temp("c,d,i,l");
//
//    L1_test("a,i") =  FIA_("i,a")
//                    + T1_("c,i") * FAB_("c,a")
//                    - T1_("a,k") * FIJ_("i,k")
//                    + T1_("c,k") * (2.0 * WIbAj_("i,c,a,k") + WIbaJ_("i,c,a,k"))
//                    + t2_temp("c,d,i,k") * WAbCi_("c,d,a,k")
//                    - t2_temp("a,c,k,l") * WKaIj_("i,c,k,l")
//                    - CGac("c,d") * (2.0 * WAkCd_("c,i,d,a") - WAkCd_("c,i,a,d"))
//                    - CGki("k,l") * (2.0 * WKlIc_("k,i,l,a") - WKlIc_("i,k,l,a"))//+ WKliC_("k,i,l,a"))
//                    ;
//    if (T1_.world().rank() == 0) {
//      std::cout << "L1_test " << L1_test << std::endl;
//    }
//
//    L2_test("a,b,i,j") =  Gabij_("a,b,i,j")
//
//                        + L1_test("c,i") * WAkCd_("c,j,a,b")
//                        + L1_test("c,j") * WAkCd_("c,i,b,a")
//
//                        - L1_test("a,k") * WKlIc_("i,j,k,b")
//                        + L1_test("b,k") * - WKlIc_("j,i,k,a")//WKliC_("i,j,k,a")
//
//                        + L1_test("a,i") * FIA_("j,b")
//                        + L1_test("b,j") * FIA_("i,a")
//
//                        + T2_("a,c,i,j") * FAB_("c,b")
//                        + T2_("c,b,i,j") * FAB_("c,a")
//
//                        - T2_("a,b,i,k") * FIJ_("j,k")
//                        - T2_("a,b,k,j") * FIJ_("i,k")
//
//                        + T2_("a,b,k,l") * WKlIj_("i,j,k,l")
//
//                        + T2_("c,d,i,j") * WAbCd_("c,d,a,b")
//
//                        + (  t2_temp("a,c,i,k") * WIbAj_("j,c,b,k")
//                           + T2_("a,c,i,k") * WIbaJ_("j,c,b,k"))
//                        + T2_("c,b,i,k") * WIbaJ_("j,c,a,k")
//                        + T2_("a,c,k,j") * WIbaJ_("i,c,b,k")
//                        + (  t2_temp("b,c,j,k") * WIbAj_("i,c,a,k")
//                           + T2_("b,c,j,k") * WIbaJ_("i,c,a,k"))
//
//                        + Gabij_("a,c,i,j") * CGac("b,c")
//                        + Gabij_("c,b,i,j") * CGac("a,c")
//
//                        - Gabij_("a,b,i,k") * CGki("k,j")
//                        - Gabij_("a,b,k,j") * CGki("k,i")
//                        ;
//    double E_L = dot(gabij_temp("a,b,i,j"), L2_test("a,b,i,j") );
//    if (T1_.world().rank() == 0) {
//      std::cout << "CCSD test E_L: " << E_L << std::endl;
//
//      //std::cout << "Gijkl_: " << Gijkl_ << std::endl;
//    }

    davidson_solver(max_iter, convergence);
    return 1.0;
  }
}


