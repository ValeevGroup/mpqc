//
// Created by Chong Peng on 7/29/15.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_CC_CC2_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_CC_CC2_H_

#include <tiledarray.h>

#include "mpqc/chemistry/qc/lcao/wfn/trange1_engine.h"

#include "../../../../../utility/cc_utility.h"
#include "mpqc/chemistry/qc/cc/diis_ccsd.h"

namespace mpqc {
namespace lcao {

template <typename Tile, typename Policy>
class CC2 {
 public:
  typedef TA::Array<double, 2, Tile, Policy> TArray2;
  typedef TA::Array<double, 3, Tile, Policy> TArray3;
  typedef TA::Array<double, 4, Tile, Policy> TArray4;

  typedef mpqc::TArrayBlock<double, 2, Tile, Policy, mpqc::MOBlock>
      TArrayBlock2;
  typedef mpqc::TArrayBlock<double, 4, Tile, Policy, mpqc::MOBlock>
      TArrayBlock4;

  typedef TA::Array<double, 4, LazyTwoElectronTile, Policy>
      DirectTwoElectronArray;

  CC2(const TArray2 &fock, const Eigen::VectorXd &ens,
      const std::shared_ptr<TRange1Engine> &tre,
      const std::shared_ptr<CCSDIntermediate<Tile, Policy>> &g)
      : ens_(ens), tre_(tre), intermediate_(g) {
    auto mo_block = std::make_shared<mpqc::MOBlock>(*tre_);
    fock_ = TArrayBlock2(fock, mo_block);
  }

  void compute_cc2() {
    auto n_occ = tre_->get_active_occ();

    TArray2 f_ai;
    f_ai("a,i") = fock_("a,i");

    auto g_abij = intermediate_->get_abij();

    TArray2 d1(f_ai.world(), f_ai.trange(), f_ai.shape(), f_ai.pmap());
    // store d1 to local
    d_ai(d1, ens_, n_occ);

    TArray4 d2(g_abij.world(), g_abij.trange(), g_abij.shape(), g_abij.pmap());
    // store d2 distributed
    d_abij_inplace(d2, ens_, n_occ);

    TArray2 t1;
    TArray4 t2;
    t1("a,i") = f_ai("a,i") * d1("a,i");
    t2("a,b,i,j") = g_abij("a,b,i,j") * d2("a,b,i,j");

    TArray4 tau;
    tau("a,b,i,j") = t2("a,b,i,j") + t1("a,i") * t1("b,j");

    double E0 = 0.0;
    double E1 =
        2.0 * TA::dot(f_ai("a,i"), t1("a,i")) +
        TA::dot(2.0 * g_abij("a,b,i,j") - g_abij("b,a,i,j"), tau("a,b,i,j"));
    double dE = std::abs(E1 - E0);

    // get all two electron integrals
    TArray4 g_ijkl = intermediate_->get_ijkl();
    TArray4 g_abcd = intermediate_->get_abcd();
    TArray4 g_iajb = intermediate_->get_iajb();
    TArray4 g_iabc = intermediate_->get_iabc();
    TArray4 g_aibc = intermediate_->get_aibc();
    TArray4 g_ijak = intermediate_->get_ijak();
    TArray4 g_ijka = intermediate_->get_ijka();

    intermediate_->clean();

    // optimize t1 and t2
    std::size_t iter = 0ul;
    while (dE >= 1.0e-12) {
      // intermediates for t1
      // external index i and a
      TArray2 h_ac, h_ki, h_ck;
      {
        h_ac("a,c") =
            -(2.0 * g_abij("c,d,k,l") - g_abij("c,d,l,k")) * tau("a,d,k,l");

        h_ki("k,i") =
            (2.0 * g_abij("c,d,k,l") - g_abij("d,c,k,l")) * tau("c,d,i,l");

        h_ck("c,k") = f_ai("c,k") +
                      (2.0 * g_abij("c,d,k,l") - g_abij("d,c,k,l")) * t1("d,l");
      }

      // update t1
      TArray2 r1;
      {
        t1("a,i") =
            d1("a,i") *
            (
                //
                f_ai("a,i") - 2.0 * f_ai("c,k") * t1("a,k") * t1("c,i")
                //
                + t1("c,i") * h_ac("a,c") - t1("a,k") * h_ki("k,i")
                //
                +
                h_ck("c,k") * (2.0 * t2("c,a,k,i") - t2("c,a,i,k") +
                               t1("c,i") * t1("a,k"))
                //
                + (2.0 * g_abij("c,a,k,i") - g_iajb("k,a,i,c")) * t1("c,k")
                //
                + (2.0 * g_iacd("k,a,c,d") - g_iacd("k,a,d,c")) * tau("c,d,k,i")
                //
                -
                (2.0 * g_klai("k,l,c,i") - g_klai("l,k,c,i")) * tau("c,a,k,l"));
      }

      // intermediates for t2
      // external index i j a b

      TArray4 a_klij, b_abcd;
      TArray2 g_ki, g_ac;
      {
        a_klij("k,l,i,j") = g_ijkl("k,l,i,j") + g_klia("k,l,i,c") * t1("c,j") +
                            g_klai("k,l,c,j") * t1("c,i") +
                            g_abij("c,d,k,l") * t1("c,i") * t1("d,j");

        b_abcd("a,b,c,d") = g_abcd("a,b,c,d") - g_aicd("a,k,c,d") * t1("b,k") -
                            g_iacd("k,b,c,d") * t1("a,k");

        //          g_ki("k,i") = h_ki("k,i") + f_ai("c,k")*t1("c,i") +
        //          (2.0*g_klia("k,l,i,a")-g_klia("l,k,i,a"))*t1("a,l");
        //          g_ac("a,c") = h_ac("a,c") - f_ai("c,k")*t1("a,k") +
        //          (2.0*g_aicd("a,k,c,d")-g_aicd("a,k,d,c"))*t1("d,k");
      }

      TArray4 r2;
      t2("a,b,i,j") =
          d2("a,b,i,j") *
          (
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
              //                + (g_ac("a,c")*t2("c,b,i,j") -
              //                g_ki("k,i")*t2("a,b,k,j"))
              //                + (g_ac("b,c")*t2("c,a,j,i") -
              //                g_ki("k,j")*t2("b,a,k,i"))

              +
              (g_iacd("i,c,a,b") - g_iajb("k,b,i,c") * t1("a,k")) * t1("c,j") +
              (g_iacd("j,c,b,a") - g_iajb("k,a,j,c") * t1("b,k")) * t1("c,i")
              //
              -
              (g_klai("i,j,a,k") + g_abij("a,c,i,k") * t1("c,j")) * t1("b,k") -
              (g_klai("j,i,b,k") + g_abij("b,c,j,k") * t1("c,i")) * t1("a,k"));

      //        t1("a,i") = r1("a,i");
      //        t2("a,b,i,j") = r2("a,b,i,j");
      tau("a,b,i,j") = t2("a,b,i,j") + t1("a,i") * t1("b,j");

      // recompute energy
      E0 = E1;
      E1 = 2.0 * TA::dot(f_ai("a,i"), t1("a,i")) +
           TA::dot((2.0 * g_abij("a,b,i,j") - g_abij("b,a,i,j")),
                   tau("a,b,i,j"));
      dE = std::abs(E0 - E1);
      iter += 1ul;
      ExEnv::out0() << iter << "  " << dE << std::endl;
    }
    ExEnv::out0() << E1 << std::endl;
  }

 private:
  Eigen::VectorXd ens_;
  std::shared_ptr<mpqc::TRange1Engine> tre_;
  std::shared_ptr<cc::CCSDIntermediate<Tile, Policy>> intermediate_;
  TArrayBlock2 fock_;
};

}  // namespace lcao
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_CC_CC2_H_
