//
// Created by Chong Peng on 11/2/17.
//

#ifndef SRC_MPQC_CHEMISTRY_QC_LCAO_CC_CCSD_R1_R2_H_
#define SRC_MPQC_CHEMISTRY_QC_LCAO_CC_CCSD_R1_R2_H_

#include <tiledarray.h>
#include "mpqc/chemistry/qc/lcao/cc/ccsd_intermediates.h"

namespace mpqc {
namespace lcao {
namespace cc {

template <typename Array>
struct Integrals {
  Integrals() = default;
  ~Integrals() = default;

  Array Fia;
  Array Fij;
  Array Fab;

  Array FIJ;
  Array FAB;

  Array Gabcd;
  Array Gijab;
  Array Gijkl;
  Array Giajb;
  Array Giabc;
  Array Gijka;
};

template <typename Array>
Array compute_cs_ccsd_r1(const Array& t1, const Array& t2, const Array& tau,
                         cc::Integrals<Array>& ints) {
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

    Array g_iabc_bar;
    g_iabc_bar("k,a,c,d") = 2.0 * ints.Giabc("k,a,c,d") - ints.Giabc("k,a,d,c");
    r1("a,i") += g_iabc_bar("k,a,c,d") * tau("c,d,k,i");

    r1("a,i") -=
        (2.0 * ints.Gijka("l,k,i,c") - ints.Gijka("k,l,i,c")) * tau("c,a,k,l");
  }
  return r1;
};


template <typename Array>
Array compute_cs_ccsd_r2(const Array& t1, const Array& t2, const Array& tau, const cc::Integrals<Array>& ints)
{

}

}  // namespace cc
}  // namespace lcao
}  // namespace mpqc

#endif  // SRC_MPQC_CHEMISTRY_QC_LCAO_CC_CCSD_R1_R2_H_
