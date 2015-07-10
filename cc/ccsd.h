//
// Created by Chong Peng on 7/1/15.
//

#ifndef TILECLUSTERCHEM_CCSD_H
#define TILECLUSTERCHEM_CCSD_H

#include "../include/tiledarray.h"
#include "../common/namespaces.h"

#include "../mp2/trange1_engine.h"
#include "../mp2/mo_block.h"
#include "../ta_routines/tarray_block.h"

#include "./two_electron_int_mo.h"


namespace tcc {

  template<typename Tile, typename Policy>
  class CCSD {

  public:
    typedef TA::Array <double, 2, Tile, Policy> TArray2;
    typedef TA::Array <double, 3, Tile, Policy> TArray3;
    typedef TA::Array <double, 4, Tile, Policy> TArray4;

    typedef tcc::TArrayBlock<double, 2, Tile, Policy, tcc::MOBlock> TArrayBlock2;
    typedef tcc::TArrayBlock<double, 4, Tile, Policy, tcc::MOBlock> TArrayBlock4;


    CCSD(const TArray2 &fock, const Eigen::VectorXd &ens,
         const std::shared_ptr<TRange1Engine> &tre,
         const std::shared_ptr<TwoElectronIntMO<Tile, Policy>> &g) :
            ens_(ens), tre_(tre), g_(g)
    {
      auto mo_block = std::make_shared<tcc::MOBlock>(*tre_);
      fock_ = TArrayBlock2(fock, mo_block);
    }

    void compute_cc2(){

      auto n_occ = tre_->get_occ();

      TArray2 f_ai;
      f_ai("a,i") = fock_("a,i");

      auto g_abij = g_->get_abij();
      std::cout << g_abij << std::endl;

      TArray2 d1 = guess_t_ai(f_ai, ens_, n_occ);
      TArray4 d2 = guess_t_abij(g_abij, ens_, n_occ);

      TArray2 t1;
      TArray4 t2;
      t1("a,i") = f_ai("a,i")*d1("a,i");
      t2("a,b,i,j") = g_abij("a,b,i,j")*d2("a,b,i,j");

//      std::cout << t1 << std::endl;
//      std::cout << t2 << std::endl;
      TArray4 tau;
      tau("a,b,i,j") = t2("a,b,i,j") + t1("a,i")*t1("b,j");

      double E0 = 0.0;
      double E1 = 2.0* TA::dot(f_ai("a,i"), t1("a,i")) +
              TA::dot(2.0*g_abij("a,b,i,j") - g_abij("b,a,i,j"), tau("a,b,i,j"));
      double dE = std::abs(E1 - E0);
      std::cout << E1 << std::endl;


      TArray2 f_ab, f_ij;
      f_ab("a,b") = fock_("a,b");
      f_ij("i,j") = fock_("i,j");
      // get all two electron integrals
      TArray4 g_ijkl = g_->get_ijkl();
      TArray4 g_abcd = g_->get_abcd();
      TArray4 g_iajb;
      g_iajb("i,a,j,b") = g_->get_aibj()("a,i,b,j");
      TArray4 g_iacd;
      g_iacd("i,a,c,d") = g_->get_abci()("a,c,d,i");
      TArray4 g_aicd;
      g_aicd("a,i,c,d") = g_->get_abic()("d,a,i,c");
      TArray4 g_klai;
      g_klai("k,l,a,i") = g_->get_aikl()("a,l,k,i");
      TArray4 g_klia;
      g_klia("k,l,i,a") = g_->get_iakl()("i,a,k,l");
      TArray4 g_iabj;
      g_iabj("i,a,b,j") = g_->get_aijb()("a,i,j,b");

      // print out g
      std::cout << g_ijkl << std::endl;
      std::cout << g_klai << std::endl;
      std::cout << g_klia << std::endl;
      std::cout << g_iabj << std::endl;
      std::cout << g_iajb << std::endl;
      std::cout << g_iacd << std::endl;
      std::cout << g_aicd << std::endl;
      std::cout << g_abcd << std::endl;

      //optimize t1 and t2
      std::size_t iter = 0ul;
      while (dE >= 1.0e-12){

        // intermediates for t1
        // external index i and a
        TArray2 h_ac, h_ki, h_ck;
        {
          h_ac("a,c") = //- f_ab("a,c")
               - (2.0 * g_abij("c,d,k,l") - g_abij("c,d,l,k")) * tau("a,d,k,l");

          h_ki("k,i") = //f_ij("k,i") +
                  (2.0 * g_abij("c,d,k,l") - g_abij("d,c,k,l")) * tau("c,d,i,l");

          h_ck("c,k") = f_ai("c,k") + (2.0 * g_abij("c,d,k,l") - g_abij("d,c,k,l")) * t1("d,l");
        }

        // update t1
        TArray2 r1;
        {
          t1("a,i") = d1("a,i")*(
                  //
                  f_ai("a,i") - 2.0 * f_ai("c,k") * t1("a,k") * t1("c,i")
                  //
                  + t1("c,i") * h_ac("a,c") - t1("a,k") * h_ki("k,i")
                  //
                  + h_ck("c,k")
                    * (2.0 * t2("c,a,k,i") - t2("c,a,i,k") + t1("c,i") * t1("a,k"))
                  //
                  + (2.0 * g_iabj("k,a,c,i") - g_iajb("k,a,i,c")) * t1("c,k")
                  //
                  + (2.0 * g_iacd("k,a,c,d") - g_iacd("k,a,d,c")) * tau("c,d,k,i")
                  //
                  - (2.0 * g_klai("k,l,c,i") - g_klai("l,k,c,i")) * tau("c,a,k,l")
          );
        }

        // intermediates for t2
        // external index i j a b

        TArray4 a_klij, b_abcd;
        TArray2 g_ki, g_ac;
        {
          a_klij("k,l,i,j") =  g_ijkl("k,l,i,j")
                               + g_klia("k,l,i,c") * t1("c,j") + g_klai("k,l,c,j") * t1("c,i")
                               + g_abij("c,d,k,l") * t1("c,i") * t1("d,j");

          b_abcd("a,b,c,d") =  g_abcd("a,b,c,d")
                               - g_aicd("a,k,c,d") * t1("b,k") - g_iacd("k,b,c,d") * t1("a,k");

//          g_ki("k,i") = h_ki("k,i") + f_ai("c,k")*t1("c,i") + (2.0*g_klia("k,l,i,a")-g_klia("l,k,i,a"))*t1("a,l");
//          g_ac("a,c") = h_ac("a,c") - f_ai("c,k")*t1("a,k") + (2.0*g_aicd("a,k,c,d")-g_aicd("a,k,d,c"))*t1("d,k");
        }

        TArray4 r2;
        t2("a,b,i,j") = d2("a,b,i,j")*(
                //
                g_abij("a,b,i,j")
                //
                + a_klij("k,l,i,j") * t1("a,k") * t1("b,l")
                //
                + b_abcd("a,b,c,d") * t1("c,i") * t1("d,j")
                //
                //+ (fab("a,c") * t2("c,b,i,j") + fab("b,c") * t2("a,c,i,j"))
                //  * (1.0 - Iab("a,c"))
                //+ (fij("k,i") * t2("a,b,k,j") + fij("k,j") * t2("a,c,i,k"))
                //  * (1 - Iij("k,i"))
                //
//                + (g_ac("a,c")*t2("c,b,i,j") - g_ki("k,i")*t2("a,b,k,j"))
//                + (g_ac("b,c")*t2("c,a,j,i") - g_ki("k,j")*t2("b,a,k,i"))

                + (g_iacd("i,c,a,b") - g_iajb("k,b,i,c") * t1("a,k")) * t1("c,j")
                + (g_iacd("j,c,b,a") - g_iajb("k,a,j,c") * t1("b,k")) * t1("c,i")
                //
                - (g_klai("i,j,a,k") + g_iabj("i,c,a,k") * t1("c,j")) * t1("b,k")
                - (g_klai("j,i,b,k") + g_iabj("j,c,b,k") * t1("c,i")) * t1("a,k")
        );

//        t1("a,i") = r1("a,i");
//        t2("a,b,i,j") = r2("a,b,i,j");
        tau("a,b,i,j") = t2("a,b,i,j") + t1("a,i") * t1("b,j");

        // recompute energy
        E0 = E1;
        E1 = 2.0 * TA::dot(f_ai("a,i"), t1("a,i"))
              + TA::dot((2.0 * g_abij("a,b,i,j") - g_abij("b,a,i,j")), tau("a,b,i,j") );
        dE = std::abs(E0 - E1);
        iter += 1ul;
        std::cout << iter << "  " << dE << std::endl;
//        std::cout << indent << scprintf("%-5.0f", iter) << scprintf("%-20.10f", Delta_E)
//        << scprintf("%-15.10f", E_1) << std::endl;


      }
      std::cout << E1 << std::endl;
    }


    void compute_ccsd(){

      auto n_occ = tre_->get_occ();

      TArray2 f_ai;
      f_ai("a,i") = fock_("a,i");

      auto g_abij = g_->get_abij();
//      std::cout << g_abij << std::endl;

      TArray2 d1 = guess_t_ai(f_ai, ens_, n_occ);
      TArray4 d2 = guess_t_abij(g_abij, ens_, n_occ);

      TArray2 t1;
      TArray4 t2;
      t1("a,i") = f_ai("a,i")*d1("a,i");
      t2("a,b,i,j") = g_abij("a,b,i,j")*d2("a,b,i,j");

//      std::cout << t1 << std::endl;
//      std::cout << t2 << std::endl;
      TArray4 tau;
      tau("a,b,i,j") = t2("a,b,i,j") + t1("a,i")*t1("b,j");

      double E0 = 0.0;
      double E1 = 2.0* TA::dot(f_ai("a,i"), t1("a,i")) +
                  TA::dot(2.0*g_abij("a,b,i,j") - g_abij("b,a,i,j"), tau("a,b,i,j"));
      double dE = std::abs(E1 - E0);
      std::cout << E1 << std::endl;


      TArray2 f_ab, f_ij;
      f_ab("a,b") = fock_("a,b");
      f_ij("i,j") = fock_("i,j");
      // get all two electron integrals
      TArray4 g_ijkl = g_->get_ijkl();
      TArray4 g_abcd = g_->get_abcd();
      TArray4 g_iajb;
      g_iajb("i,a,j,b") = g_->get_aibj()("a,i,b,j");
      TArray4 g_iacd;
      g_iacd("i,a,c,d") = g_->get_abci()("a,c,d,i");
      TArray4 g_aicd;
      g_aicd("a,i,c,d") = g_->get_abic()("d,a,i,c");
      TArray4 g_klai;
      g_klai("k,l,a,i") = g_->get_aikl()("a,l,k,i");
      TArray4 g_klia;
      g_klia("k,l,i,a") = g_->get_iakl()("i,a,k,l");
      TArray4 g_iabj;
      g_iabj("i,a,b,j") = g_->get_aijb()("a,i,j,b");

      // print out g
      std::cout << g_ijkl << std::endl;
//      std::cout << g_klai << std::endl;
//      std::cout << g_klia << std::endl;
//      std::cout << g_iabj << std::endl;
//      std::cout << g_iajb << std::endl;
//      std::cout << g_iacd << std::endl;
//      std::cout << g_aicd << std::endl;
//      std::cout << g_abcd << std::endl;

      //optimize t1 and t2
      std::size_t iter = 0ul;
      while (dE >= 1.0e-12){

        // intermediates for t1
        // external index i and a
        TArray2 h_ac, h_ki, h_kc;
        {
          h_ac("a,c") = //- f_ab("a,c")
                        - (2.0 * g_abij("c,d,k,l") - g_abij("c,d,l,k")) * tau("a,d,k,l");

          h_ki("k,i") = // f_ij("k,i") +
                  (2.0 * g_abij("c,d,k,l") - g_abij("d,c,k,l")) * tau("c,d,i,l");

          h_kc("k,c") = f_ai("c,k") + (2.0 * g_abij("c,d,k,l") - g_abij("d,c,k,l")) * t1("d,l");
        }

        // update t1
        TArray2 r1;
        {
          t1("a,i") = d1("a,i")*(
                  //
                  f_ai("a,i") - 2.0 * f_ai("c,k") * t1("a,k") * t1("c,i")
                  //
                  + h_ac("a,c")*t1("c,i") - t1("a,k")*h_ki("k,i")
                  //
                  + h_kc("k,c")
                    * (2.0 * t2("c,a,k,i") - t2("c,a,i,k") + t1("c,i") * t1("a,k"))
                  //
                  + (2.0 * g_iabj("k,a,c,i") - g_iajb("k,a,i,c")) * t1("c,k")
                  //
                  + (2.0 * g_iacd("k,a,c,d") - g_iacd("k,a,d,c")) * tau("c,d,k,i")
                  //
                  - (2.0 * g_klai("k,l,c,i") - g_klai("l,k,c,i")) * tau("c,a,k,l")
          );
        }

        // intermediates for t2
        // external index i j a b

        TArray4 a_klij, b_abcd, j_akic, k_kaic, T;
        TArray2 g_ki, g_ac;
        {
          T("d,b,i,l") = 0.5*t2("d,b,i,l") + t1("d,i")*t1("b,l");

          a_klij("k,l,i,j") =  g_ijkl("k,l,i,j")
                               + g_klia("k,l,i,c") * t1("c,j") + g_klai("k,l,c,j") * t1("c,i")
                               + g_abij("c,d,k,l") * tau("c,d,i,j");

          b_abcd("a,b,c,d") =  g_abcd("a,b,c,d")
                               - g_aicd("a,k,c,d") * t1("b,k") - g_iacd("k,b,c,d") * t1("a,k");

          g_ki("k,i") = h_ki("k,i") + f_ai("c,k")*t1("c,i")
                        + (2.0*g_klia("k,l,i,c")-g_klia("l,k,i,c"))*t1("c,l");

          g_ac("a,c") = h_ac("a,c") - f_ai("c,k")*t1("a,k")
                        + (2.0*g_aicd("a,k,c,d")-g_aicd("a,k,d,c"))*t1("d,k");

          j_akic("a,k,i,c") = g_abij("a,c,i,k") - g_klia("l,k,i,c")*t1("a,l")
                              + g_aicd("a,k,d,c")*t1("d,i") - g_abij("c,d,k,l")*T("d,a,i,l")
                              + 0.5*(2.0*g_abij("c,d,k,l")- g_abij("d,c,k,l"))*t2("a,d,i,l");

          k_kaic("k,a,i,c") = g_iajb("k,a,i,c") - g_klia("k,l,i,c")*t1("a,l") +
                  g_iacd("k,a,c,d")*t1("d,i") - g_abij("d,c,k,l")*T("d,a,i,l");

        }

        TArray4 t2_unsymm;
        t2_unsymm("a,b,i,j") = d2("a,b,i,j")*(
                //
                g_abij("a,b,i,j")
                //
                + a_klij("k,l,i,j") * tau("a,b,k,l")
                //
                + b_abcd("a,b,c,d") * tau("c,d,i,j")

                // permutation part
                //
                + (g_ac("a,c")*t2("c,b,i,j") - g_ki("k,i")*t2("a,b,k,j"))
                + (g_ac("b,c")*t2("c,a,j,i") - g_ki("k,j")*t2("b,a,k,i"))

                + (g_iacd("i,c,a,b") - g_iajb("k,b,i,c") * t1("a,k")) * t1("c,j")
                + (g_iacd("j,c,b,a") - g_iajb("k,a,j,c") * t1("b,k")) * t1("c,i")
                //
                - (g_klai("i,j,a,k") + g_abij("a,c,i,k") * t1("c,j")) * t1("b,k")
                - (g_klai("j,i,b,k") + g_abij("b,c,j,k") * t1("c,i")) * t1("a,k")

                + 0.5*(2.0*j_akic("a,k,i,c") - k_kaic("k,a,i,c")) * (2.0*t2("c,b,k,j") - t2("b,c,k,j"))
                + 0.5*(2.0*j_akic("b,k,j,c") - k_kaic("k,b,j,c")) * (2.0*t2("c,a,k,i") - t2("a,c,k,i"))

                - 0.5*k_kaic("k,a,i,c")*t2("b,c,k,j") - k_kaic("k,b,i,c")*t2("a,c,k,j")
                - 0.5*k_kaic("k,b,j,c")*t2("a,c,k,i") - k_kaic("k,a,j,c")*t2("b,c,k,i")
        );

        t2("a,b,i,j") = 0.5*(t2_unsymm("a,b,i,j") + t2_unsymm("b,a,i,j"));
//        t1("a,i") = r1("a,i");
//        t2("a,b,i,j") = r2("a,b,i,j");
        tau("a,b,i,j") = t2("a,b,i,j") + t1("a,i") * t1("b,j");

        // recompute energy
        E0 = E1;
        E1 = 2.0 * TA::dot(f_ai("a,i"), t1("a,i"))
             + TA::dot((2.0 * g_abij("a,b,i,j") - g_abij("b,a,i,j")), tau("a,b,i,j") );
        dE = std::abs(E0 - E1);
        iter += 1ul;
        std::cout << iter << "  " << dE << "  " << E1 << std::endl;
//        std::cout << indent << scprintf("%-5.0f", iter) << scprintf("%-20.10f", Delta_E)
//        << scprintf("%-15.10f", E_1) << std::endl;

      }
      std::cout << E1 << std::endl;
    }

  private:

    TArray4 guess_t_abij(const TArray4& abij,
                        const Eigen::VectorXd& ens, std::size_t n_occ)
    {
      auto convert = [&ens, n_occ](Tile &result_tile, const Tile &arg_tile) {
        result_tile = Tile(arg_tile.range());

        // compute index
        const auto a0 = arg_tile.range().lobound()[0];
        const auto an = arg_tile.range().upbound()[0];
        const auto b0 = arg_tile.range().lobound()[1];
        const auto bn = arg_tile.range().upbound()[1];
        const auto i0 = arg_tile.range().lobound()[2];
        const auto in = arg_tile.range().upbound()[2];
        const auto j0 = arg_tile.range().lobound()[3];
        const auto jn = arg_tile.range().upbound()[3];

        auto tile_idx = 0;
        typename Tile::value_type norm = 0.0;
        for (auto a = a0; a < an; ++a) {
          const auto e_a = ens[a + n_occ];
          for (auto b = b0; b < bn; ++b) {
            const auto e_b = ens[b + n_occ];
            for (auto i = i0; i < in; ++i) {
              const auto e_i = ens[i];
              for (auto j = j0; j < jn; ++j, ++tile_idx) {
                const auto e_j = ens[j];
                const auto e_iajb = e_i + e_j - e_a - e_b;
//                const auto result_abij = arg_tile[tile_idx]/(e_iajb);
                const auto result_abij = 1.0/(e_iajb);
                norm += result_abij*result_abij;
                result_tile[tile_idx] = result_abij;
              }
            }
          }
        }
        return std::sqrt(norm);
      };

      return TA::foreach(abij, convert);
    }

    TArray2 guess_t_ai(const TArray2& f_ai, const Eigen::VectorXd& ens, int n_occ)
    {
      auto convert = [&ens, n_occ] (Tile& result_tile, const Tile& arg_tile){
        result_tile =  Tile(arg_tile.range());
        const auto a0 = arg_tile.range().lobound()[0];
        const auto an = arg_tile.range().upbound()[0];
        const auto i0 = arg_tile.range().lobound()[1];
        const auto in = arg_tile.range().upbound()[1];

        auto ai = 0;
        typename Tile::value_type norm = 0.0;
        for (auto a = a0; a < an; ++a) {
          const auto e_a = ens[a + n_occ];
          for (auto i = i0; i < in; ++i, ++ai) {
            const auto e_i = ens[i];
            const auto e_ia = e_i - e_a;
//            const auto result_ai = arg_tile[ai] / (e_ia);
            const auto result_ai = 1.0 / (e_ia);
            norm += result_ai * result_ai;
            result_tile[ai] = result_ai;
          }
        }
        return std::sqrt(norm);
      };

      TArray2 t_ai = TA::foreach(f_ai, convert);
      return  t_ai;
    }

  private:
    Eigen::VectorXd ens_;
    std::shared_ptr<tcc::TRange1Engine> tre_;
    std::shared_ptr<tcc::TwoElectronIntMO<Tile, Policy>> g_;
    TArrayBlock2 fock_;
  };
}


#endif //TILECLUSTERCHEM_CCSD_H
