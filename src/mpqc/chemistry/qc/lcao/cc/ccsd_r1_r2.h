//
// Created by Chong Peng on 11/2/17.
//

#ifndef SRC_MPQC_CHEMISTRY_QC_LCAO_CC_CCSD_R1_R2_H_
#define SRC_MPQC_CHEMISTRY_QC_LCAO_CC_CCSD_R1_R2_H_

#include <tiledarray.h>
#include "mpqc/chemistry/qc/lcao/cc/ccsd_intermediates.h"
#include "mpqc/chemistry/qc/lcao/integrals/direct_task_integrals.h"

namespace mpqc {
namespace lcao {
namespace cc {

/**
 * @brief this file contains functions to compute spin-adapted closed-shell CCSD
 * amplitude equation
 */

/**
 * Structure to hold integrals needed to compute CCSD amplitudes
 */
template <typename Array>
struct Integrals {
  Integrals() = default;
  ~Integrals() = default;

  Array Fia;  // fock matrix <i|F|a>
  Array Fij;  // fock matrix <i|F|j>
  Array Fab;  // fock matrix <a|F|a>

  // effective one body Hamiltonian
  Array FIJ;  // <i|H_bar|j>
  Array FAB;  // <a|H_bar|b>

  // mo two electron integrals
  Array Gabcd;  // <a b|G|c d>
  Array Gijab;  // <i j|G|a b>
  Array Gijkl;  // <i j|G|k l>
  Array Giajb;  // <i a|G|j b>
  Array Giabc;  // <i a|G|b c>
  Array Gijka;  // <i j|G|k a>

  // mo three center integrals
  Array Xai;  // (X|a i)(X|K)^-1/2
  Array Xij;  // (X|i j)(X|K)^-1/2
  Array Xab;  // (X|a b)(X|K)^-1/2

  // mo coefficients
  Array Ci;
  Array Ca;

  // CP decompositions of mo three center integrals
  std::vector<Array> Xab_factors;
};

/**
 * computes closed-shell CCSD R1 residual
 *
 * @param t1 CCSD T1 amplitudes in T1("a,i")
 * @param t2 CCSD T2 amplitudes in T2("a,b,i,j")
 * @param tau T2("a,b,i,j") + T1("a,i")*T1("b,j")
 * @param ints cc::Integrals, requires Fia, Fij, Fab, Gijab, Giajb, Gijka and
 * Giabc if u is not initialized
 * @param u half transformed intermediates U("p,r,i,j") =
 * (Tau("a,b,i,j")*Ca("q,a")*Ca("s,b"))*(p q|r s)
 * @return R1 residual, and update FIJ and FAB in ints
 */
template <typename Array>
Array compute_cs_ccsd_r1(const Array& t1, const Array& t2, const Array& tau,
                         cc::Integrals<Array>& ints, const Array& u = Array()) {
  Array r1;
  {
    // intermediates for t1
    // external index i and a
    // vir index a b c d
    // occ index i j k l
    Array h_kc, h_ki, h_ac;

    // compute residual r1(n) (at convergence r1 = 0)
    // external index i and a
    r1("a,i") =
        ints.Fia("i,a") - 2.0 * (ints.Fia("k,c") * t1("c,i")) * t1("a,k");

    Array g_ijab_bar;
    g_ijab_bar("i,j,a,b") = 2.0 * ints.Gijab("i,j,a,b") - ints.Gijab("j,i,a,b");
    {
      h_ac = compute_cs_ccsd_F_vv(ints.Fab, g_ijab_bar, tau);
      ints.FAB = h_ac;
      r1("a,i") += h_ac("a,c") * t1("c,i");
    }

    {
      h_ki = compute_cs_ccsd_F_oo(ints.Fij, g_ijab_bar, tau);
      ints.FIJ = h_ki;
      r1("a,i") -= t1("a,k") * h_ki("k,i");
    }

    {
      h_kc = compute_cs_ccsd_F_ov(ints.Fia, g_ijab_bar, t1);

      r1("a,i") += h_kc("k,c") * (2.0 * t2("c,a,k,i") - t2("c,a,i,k") +
                                  t1("c,i") * t1("a,k"));
    }

    r1("a,i") +=
        (2.0 * ints.Gijab("k,i,c,a") - ints.Giajb("k,a,i,c")) * t1("c,k");

    if (!u.is_initialized()) {
      Array g_iabc_bar;
      g_iabc_bar("k,a,c,d") =
          2.0 * ints.Giabc("k,a,c,d") - ints.Giabc("k,a,d,c");
      r1("a,i") += g_iabc_bar("k,a,c,d") * tau("c,d,k,i");
    } else {
      r1("a,i") +=
          (2.0 * u("p,r,k,i") - u("p,r,i,k")) * ints.Ci("p,k") * ints.Ca("r,a");
    }

    r1("a,i") -=
        (2.0 * ints.Gijka("l,k,i,c") - ints.Gijka("k,l,i,c")) * tau("c,a,k,l");
  }
  return r1;
};

/**
 * computes closed-shell CCSD R1 residual with density-fitting
 *
 * @param t1 CCSD T1 amplitudes in T1("a,i")
 * @param t2 CCSD T2 amplitudes in T2("a,b,i,j")
 * @param tau T2("a,b,i,j") + T1("a,i")*T1("b,j")
 * @param ints cc::Integrals, requires Fia, Fij, Fab, Gijab, Giajb, Xai, Xab,
 * Xij
 * @return R1 residual, and update FIJ and FAB in ints
 */
template <typename Array>
Array compute_cs_ccsd_r1_df(const Array& t1, const Array& t2, const Array& tau,
                            cc::Integrals<Array>& ints) {
  Array r1;
  {
    // intermediates for t1
    // external index i and a
    // vir index a b c d
    // occ index i j k l
    Array h_kc, h_ki, h_ac;

    Array X_ai_tau;

    X_ai_tau("X,a,i") = 2.0 * ints.Xai("X,b,j") * tau("a,b,i,j") -
                        ints.Xai("X,b,j") * tau("a,b,j,i");

    // compute residual r1(n) (at convergence r1 = 0)
    // external index i and a
    r1("a,i") = ints.Fia("i,a") - 2.0 * ints.Fia("k,c") * t1("c,i") * t1("a,k");

    {
      h_ac("a,c") = ints.Fab("a,c") - X_ai_tau("K,a,l") * ints.Xai("K,c,l");
      ints.FAB = h_ac;
      r1("a,i") += h_ac("a,c") * t1("c,i");
    }

    {
      h_ki("k,i") = ints.Fij("k,i") + X_ai_tau("K,c,i") * ints.Xai("K,c,k");
      ints.FIJ = h_ki;
      r1("a,i") -= t1("a,k") * h_ki("k,i");
    }

    {
      h_kc("k,c") =
          ints.Fia("k,c") +
          (2.0 * ints.Gijab("k,l,c,d") - ints.Gijab("k,l,d,c")) * t1("d,l");
      r1("a,i") += h_kc("k,c") * (2.0 * t2("c,a,k,i") - t2("c,a,i,k") +
                                  t1("c,i") * t1("a,k"));
    }

    r1("a,i") +=
        (2.0 * ints.Gijab("k,i,c,a") - ints.Giajb("k,a,i,c")) * t1("c,k");

    r1("a,i") += X_ai_tau("K,c,i") * ints.Xab("K,a,c");

    r1("a,i") -= X_ai_tau("K,a,l") * ints.Xij("K,l,i");
  }
  return r1;
}

/**
 * computes closed-shell CCSD R2 residual
 *
 * @param t1 CCSD T1 amplitudes in T1("a,i")
 * @param t2 CCSD T2 amplitudes in T2("a,b,i,j")
 * @param tau T2("a,b,i,j") + T1("a,i")*T1("b,j")
 * @param ints cc::Integrals, requires Fia, FIJ, FAB, Gijab, Giajb, Gijka,
 * Gijkl, Giabc and Gabcd if u is not initialized
 * @param u half transformed intermediates U("p,r,i,j") =
 * (Tau("a,b,i,j")*Ca("q,a")*Ca("s,b"))*(p q|r s)
 * @return R2 residual
 */
template <typename Array>
Array compute_cs_ccsd_r2(const Array& t1, const Array& t2, const Array& tau,
                         const cc::Integrals<Array>& ints,
                         const Array& u = Array()) {
  Array r2;

  // compute residual r2(n) (at convergence r2 = 0)
  // permutation part
  {
    r2("a,b,i,j") =
        (ints.Giabc("i,c,a,b") - ints.Giajb("k,b,i,c") * t1("a,k")) * t1("c,j");

    r2("a,b,i,j") -=
        (ints.Gijka("j,i,k,a") + ints.Gijab("i,k,a,c") * t1("c,j")) * t1("b,k");
  }

  {
    // compute g intermediate
    Array g_ki, g_ac;

    Array g_ijka_bar;
    g_ijka_bar("i,j,k,a") = 2.0 * ints.Gijka("i,j,k,a") - ints.Gijka("j,i,k,a");
    g_ki = cc::compute_cs_ccsd_F_OO(ints.FIJ, ints.Fia, g_ijka_bar, t1);

    Array g_iabc_bar;
    g_iabc_bar("k,a,d,b") = 2.0 * ints.Giabc("k,a,d,b") - ints.Giabc("k,a,b,d");
    g_ac = cc::compute_cs_ccsd_F_VV(ints.FAB, ints.Fia, g_iabc_bar, t1);

    r2("a,b,i,j") += g_ac("a,c") * t2("c,b,i,j") - g_ki("k,i") * t2("a,b,k,j");
  }

  {
    Array j_akic;
    Array k_kaic;
    // compute j and k intermediate
    {
      Array T;

      T("d,b,i,l") = 0.5 * t2("d,b,i,l") + t1("d,i") * t1("b,l");

      j_akic("a,k,i,c") =
          ints.Gijab("i,k,a,c") - ints.Gijka("l,k,i,c") * t1("a,l") +
          ints.Giabc("k,a,c,d") * t1("d,i") -
          ints.Gijab("k,l,c,d") * T("d,a,i,l") +
          0.5 * (2.0 * ints.Gijab("k,l,c,d") - ints.Gijab("k,l,d,c")) *
              t2("a,d,i,l");

      k_kaic("k,a,i,c") = ints.Giajb("k,a,i,c") -
                          ints.Gijka("k,l,i,c") * t1("a,l") +
                          ints.Giabc("k,a,d,c") * t1("d,i") -
                          ints.Gijab("k,l,d,c") * T("d,a,i,l");
    }

    r2("a,b,i,j") += 0.5 * (2.0 * j_akic("a,k,i,c") - k_kaic("k,a,i,c")) *
                     (2.0 * t2("c,b,k,j") - t2("b,c,k,j"));

    Array tmp;
    tmp("a,b,i,j") = k_kaic("k,a,i,c") * t2("c,b,j,k");
    r2("a,b,i,j") -= 0.5 * tmp("a,b,i,j") + tmp("b,a,i,j");
  }

  // perform the permutation
  r2("a,b,i,j") = r2("a,b,i,j") + r2("b,a,j,i");

  r2("a,b,i,j") += ints.Gijab("i,j,a,b");

  {
    // compute a intermediate
    Array a_klij =
        cc::compute_cs_ccsd_W_oooo(ints.Gijkl, ints.Gijka, ints.Gijab, tau, t1);

    r2("a,b,i,j") += a_klij("k,l,i,j") * tau("a,b,k,l");
  }

  {
    Array b_abij;
    if (!u.is_initialized()) {
      b_abij("a,b,i,j") = tau("c,d,i,j") * ints.Gabcd("a,b,c,d");
      Array tmp;
      tmp("k,a,i,j") = ints.Giabc("k,a,c,d") * tau("c,d,i,j");
      b_abij("a,b,i,j") -= tmp("k,a,j,i") * t1("b,k");

      b_abij("a,b,i,j") -= tmp("k,b,i,j") * t1("a,k");
    } else {
      b_abij("a,b,i,j") =
          (u("p,r,i,j") * ints.Ca("r,b") -
           ints.Ci("r,k") * t1("b,k") * u("p,r,i,j")) *
              ints.Ca("p,a") -
          u("p,r,i,j") * ints.Ci("p,k") * ints.Ca("r,b") * t1("a,k");
    }

    r2("a,b,i,j") += b_abij("a,b,i,j");
  }

  return r2;
}

/**
 * computes closed-shell CCSD R2 residual with density-fitting
 *
 * @param t1 CCSD T1 amplitudes in T1("a,i")
 * @param t2 CCSD T2 amplitudes in T2("a,b,i,j")
 * @param tau T2("a,b,i,j") + T1("a,i")*T1("b,j")
 * @param ints cc::Integrals, requires Fia, FIJ, FAB, Xai, Xab, Xij, Gijab,
 * Giajb, Gijka, Gijkl  and  Giabc with Gabcd
 * if Giabc and Gabcd is not initialized, it will evaluated with lazy density-fitting
 * if U is initialized, Giabc and Gabcd won't be used
 * @param u half transformed intermediates U("p,r,i,j") =
 * (Tau("a,b,i,j")*Ca("q,a")*Ca("s,b"))*(p q|r s)
 * @return R2 residual
 */
template <typename Array>
Array compute_cs_ccsd_r2_df(const Array& t1, const Array& t2, const Array& tau,
                            const cc::Integrals<Array>& ints,
                            const Array& u = Array()) {
  // compute residual r2(n) (at convergence r2 = 0)

  Array r2;
  // permutation part
  Array X_ab_t1;
  X_ab_t1("K,a,i") = ints.Xab("K,a,b") * t1("b,i");

  {
    r2("a,b,i,j") =
        X_ab_t1("K,b,j") * (ints.Xai("K,a,i") - ints.Xij("K,k,i") * t1("a,k")) -
        (ints.Gijka("j,i,k,a") + ints.Gijab("i,k,a,c") * t1("c,j")) * t1("b,k");
  }
  {
    // compute g intermediate
    Array g_ki, g_ac;

    Array g_ijka_bar;
    g_ijka_bar("i,j,k,a") = 2.0 * ints.Gijka("i,j,k,a") - ints.Gijka("j,i,k,a");
    g_ki = cc::compute_cs_ccsd_F_OO(ints.FIJ, ints.Fia, g_ijka_bar, t1);

    g_ac("a,c") = ints.FAB("a,c") - ints.Fia("k,c") * t1("a,k") +
                  2.0 * ints.Xai("K,d,k") * t1("d,k") * ints.Xab("K,a,c") -
                  X_ab_t1("K,a,k") * ints.Xai("K,c,k");

    r2("a,b,i,j") += g_ac("a,c") * t2("c,b,i,j") - g_ki("k,i") * t2("a,b,k,j");
  }
  {
    Array j_akic;
    Array k_kaic;
    // compute j and k intermediate
    {
      Array T;

      T("d,b,i,l") = 0.5 * t2("d,b,i,l") + t1("d,i") * t1("b,l");

      j_akic("a,k,i,c") =
          ints.Gijab("i,k,a,c") - ints.Gijka("l,k,i,c") * t1("a,l") +
          X_ab_t1("K,a,i") * ints.Xai("K,c,k") -
          ints.Xai("x,d,l") * T("d,a,i,l") * ints.Xai("x,c,k") +
          0.5 * (2.0 * ints.Gijab("k,l,c,d") - ints.Gijab("k,l,d,c")) *
              t2("a,d,i,l");

      k_kaic("k,a,i,c") = ints.Giajb("k,a,i,c") -
                          ints.Gijka("k,l,i,c") * t1("a,l") +
                          ints.Xai("K,d,k") * t1("d,i") * ints.Xab("K,a,c") -
                          ints.Gijab("k,l,d,c") * T("d,a,i,l");
    }

    r2("a,b,i,j") += 0.5 * (2.0 * j_akic("a,k,i,c") - k_kaic("k,a,i,c")) *
                     (2.0 * t2("c,b,k,j") - t2("b,c,k,j"));

    Array tmp;
    tmp("a,b,i,j") = k_kaic("k,a,i,c") * t2("c,b,j,k");
    r2("a,b,i,j") -= 0.5 * tmp("a,b,i,j") + tmp("b,a,i,j");
  }

  // perform the permutation
  r2("a,b,i,j") = r2("a,b,i,j") + r2("b,a,j,i");
  r2("a,b,i,j") += ints.Gijab("i,j,a,b");

  {
    // compute a intermediate
    Array a_klij =
        cc::compute_cs_ccsd_W_oooo(ints.Gijkl, ints.Gijka, ints.Gijab, tau, t1);

    r2("a,b,i,j") += a_klij("k,l,i,j") * tau("a,b,k,l");
  }

  {
    // compute b intermediate
    // avoid store b_abcd
    Array b_abij;

    if (!u.is_initialized()) {
      auto time1 = std::chrono::high_resolution_clock::now();
      auto time2 = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> time_span = time2 - time1;

      if (ints.Gabcd.is_initialized()) {
        time1 = std::chrono::high_resolution_clock::now();
        b_abij("a,b,i,j") = tau("c,d,i,j") * ints.Gabcd("a,b,c,d");
        time2 = std::chrono::high_resolution_clock::now();
        time_span = time2 - time1;

        double old_time = time_span.count();
        //std::cout << "Time to compute G with contraction is " << old_time << std::endl;

        Array tmp;
        tmp("k,a,i,j") = ints.Giabc("k,a,c,d") * tau("c,d,i,j");
        b_abij("a,b,i,j") -= tmp("k,a,j,i") * t1("b,k");
        b_abij("a,b,i,j") -= tmp("k,b,i,j") * t1("a,k");

      }
      // if cp_ccsd_ compute the contraction of t_ijcd g_cdab with factor matrices
      else if(!ints.Xab_factors.empty()) {
        time1 = std::chrono::high_resolution_clock::now();
        Array temp, Zrrp;
        Zrrp("r, rp") = ints.Xab_factors[0]("X, r") * ints.Xab_factors[0]("X, rp");

        auto &world = Zrrp.world();
        Zrrp.make_replicated();
        world.gop.fence();

        temp("i, j, r, rp") = (tau("c,d,i,j") * ints.Xab_factors[1]("c, r")) * ints.Xab_factors[1]("d, rp");

        auto text = temp.trange().tiles_range().extent_data();
        auto ni = text[0]; // tile range for ij dimension
        auto nr = text[2]; // tile range for the r dimension

        // for each tile in temp scale said tile by a value in Z(r,r')
        auto hadamard = [&Zrrp](int r, int rp, TA::Tensor<double> tile) {
          // if the scaling factor is 0 fill the tile with 0s
          if(Zrrp.is_zero({r,rp})){
            const auto volume = tile.range().volume();
            std::fill(tile.data(), tile.data() + volume, 0.0);
            return 0.0;
          }
          // use future to get the correct tile of Z(r r')
          auto Zt = Zrrp.find({r, rp}).get();

          auto lo = tile.range().lobound_data();
          auto up = tile.range().upbound_data();

          for (auto k = lo[0]; k < up[0]; ++k) {
            for (auto l = lo[1]; l < up[1]; ++l) {
              for(auto R = lo[2]; R < up[2]; ++R) {
                for (auto Rp = lo[3]; Rp < up[3]; ++Rp) {
                  tile(k,l,R,Rp) *= Zt(R,Rp);
                }
              }
            }
          }
          return tile.norm();
        };

        //TODO Foreach in place
        // These loops go over the tiles of temp and compute the scaling
        for (int i = 0; i < ni; ++i) {
          for (int j = 0; j < ni; ++j) {
            for (int r = 0; r < nr; ++r) {
              for (int rp = 0; rp < nr; ++rp) {

                // If the processor doesn't have the tile of temp or the tile of temp is 0 it is skipped
                if (!temp.is_local({i, j, r, rp}) || temp.is_zero({i, j, r, rp})) {
                  continue;
                }

                // Future grabs the tile
                auto tile = temp.find({i, j, r, rp});

                // tasks of hadamard on tiles are distributed in parallel processing.
                world.taskq.add(hadamard, r, rp, tile);

              }
            }
          }
        }

        world.gop.fence();
        temp.truncate();

        b_abij("a,b,i,j") = (temp("i, j, r, rp") * ints.Xab_factors[2]("b, rp")) * ints.Xab_factors[2]("a, r");

        time2 = std::chrono::high_resolution_clock::now();
        time_span = time2 - time1;

        double old_time = time_span.count();
        //std::cout << "Time to compute G with CP is " << old_time << std::endl;

        Array tmp;

        if(!ints.Giabc.is_initialized()){
          auto giabc = gaussian::df_direct_integrals(
                  ints.Xai, ints.Xab);
          tmp("k, a, i, j") = giabc("c,k,a,d") * tau("c,d,i,j");
          b_abij("a,b,i,j") -= tmp("k,a,j,i") * t1("b,k");
          b_abij("a,b,i,j") -= tmp("k,b,i,j") * t1("a,k");

        }else {
          tmp("k,a,i,j") = ints.Giabc("k,a,c,d") * tau("c,d,i,j");
          b_abij("a,b,i,j") -= tmp("k,a,j,i") * t1("b,k");
          b_abij("a,b,i,j") -= tmp("k,b,i,j") * t1("a,k");
        }
      } else {
        Array X_ab_t1;
        X_ab_t1("K,a,b") =
            0.5 * ints.Xab("K,a,b") - ints.Xai("K,b,i") * t1("a,i");

        auto g_abcd_iabc_direct = gaussian::df_direct_integrals(
            ints.Xab, X_ab_t1, Formula::Notation::Physical);

        b_abij("a,b,i,j") = tau("c,d,i,j") * g_abcd_iabc_direct("a,b,c,d");

        b_abij("a,b,i,j") += b_abij("b,a,j,i");
      }
    } else {
      b_abij("a,b,i,j") =
          (u("p,r,i,j") * ints.Ca("r,b") -
           ints.Ci("r,k") * t1("b,k") * u("p,r,i,j")) *
              ints.Ca("p,a") -
          u("p,r,i,j") * ints.Ci("p,k") * ints.Ca("r,b") * t1("a,k");
    }

    r2("a,b,i,j") += b_abij("a,b,i,j");
  }

  return r2;
}

}  // namespace cc
}  // namespace lcao
}  // namespace mpqc

#endif  // SRC_MPQC_CHEMISTRY_QC_LCAO_CC_CCSD_R1_R2_H_
