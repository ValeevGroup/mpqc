//
// Created by Chong Peng on 3/21/17.
//

#ifndef SRC_MPQC_CHEMISTRY_QC_LCAO_CC_CCSD_HBAR_H_
#define SRC_MPQC_CHEMISTRY_QC_LCAO_CC_CCSD_HBAR_H_

/**
 * @brief This file contains functions to compute effective Hamiltonian of CCSD
 *
 * The spin orbital equations for intermediates can be found in
 * Gauss, J., & Stanton, J. F. (1995). The Journal of Chemical Physics, 103(9),
 * 3561–3577.
 *
 */

#include <tuple>

#include <tiledarray.h>

#include "mpqc/chemistry/qc/lcao/factory/lcao_factory.h"

namespace mpqc {
namespace lcao {
namespace cc {

/**
 *
 */
template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> compute_cs_ccsd_W_KlIj(
    LCAOFactoryBase<Tile, Policy>& lcao_factory,
    const TA::DistArray<Tile, Policy>& t1,
    const TA::DistArray<Tile, Policy>& tau, bool df) {
  std::wstring postfix = df ? L"[df]" : L"";
  TA::DistArray<Tile, Policy> W_KlIj;
  auto g_ijka = lcao_factory.compute(L"<i j|G|k a>" + postfix);
  auto g_ijkl = lcao_factory.compute(L"<i j|G|k l>" + postfix);
  auto g_ijab = lcao_factory.compute(L"<i j|G|a b>" + postfix);

  W_KlIj("k,l,i,j") = t1("c,j") * g_ijka("k,l,i,c");
  W_KlIj("k,l,i,j") = W_KlIj("k,l,i,j") + W_KlIj("l,k,j,i") +
                      g_ijkl("k,l,i,j") + g_ijab("k,l,c,d") * tau("c,d,i,j");
  return W_KlIj;
};

/*
 * //TODO need to permute to IjAb
 */
template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> compute_cs_ccsd_W_IbAj(
    LCAOFactoryBase<Tile, Policy>& lcao_factory,
    const TA::DistArray<Tile, Policy>& t1,
    const TA::DistArray<Tile, Policy>& t2, bool df) {
  std::wstring postfix = df ? L"[df]" : L"";

  TA::DistArray<Tile, Policy> W_IbAj;

  auto g_ijab = lcao_factory.compute(L"<i j|G|a b>" + postfix);
  auto g_ijka = lcao_factory.compute(L"<i j|G|k a>" + postfix);
  auto g_iabc = lcao_factory.compute(L"<i a|G|b c>" + postfix);

  W_IbAj("i,b,a,j") =
      g_ijab("i,j,a,b") + t1("d,j") * g_iabc("i,b,a,d") -
      t1("b,l") * g_ijka("l,i,j,a") +
      t2("b,d,j,l") * (2.0 * g_ijab("i,l,a,d") - g_ijab("l,i,a,d")) -
      t2("d,b,j,l") * g_ijab("i,l,a,d") -
      t1("d,j") * (t1("b,l") * g_ijab("i,l,a,d"));

  return W_IbAj;
};

template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> compute_cs_ccsd_W_IbaJ(
    LCAOFactoryBase<Tile, Policy>& lcao_factory,
    const TA::DistArray<Tile, Policy>& t1,
    const TA::DistArray<Tile, Policy>& t2, bool df) {
  std::wstring postfix = df ? L"[df]" : L"";

  TA::DistArray<Tile, Policy> W_IbaJ;

  auto g_ijab = lcao_factory.compute(L"<i j|G|a b>" + postfix);
  auto g_iajb = lcao_factory.compute(L"<i a|G|j b>" + postfix);
  auto g_ijka = lcao_factory.compute(L"<i j|G|k a>" + postfix);
  auto g_iabc = lcao_factory.compute(L"<i a|G|b c>" + postfix);

  W_IbaJ("i,b,a,j") = -g_iajb("j,a,i,b") - t1("d,j") * g_iabc("i,b,d,a") +
                      t1("b,l") * g_ijka("i,l,j,a") +
                      t2("d,b,j,l") * g_ijab("i,l,d,a") +
                      t1("d,j") * (t1("b,l") * g_ijab("i,l,d,a"));

  return W_IbaJ;
};

template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> compute_cs_ccsd_W_AkCd(
    LCAOFactoryBase<Tile, Policy>& lcao_factory,
    const TA::DistArray<Tile, Policy>& t1, bool df) {
  std::wstring postfix = df ? L"[df]" : L"";

  TA::DistArray<Tile, Policy> W_AkCd;

  auto g_ijab = lcao_factory.compute(L"<i j|G|a b>" + postfix);
  auto g_iabc = lcao_factory.compute(L"<i a|G|b c>" + postfix);

  W_AkCd("a,k,c,d") = g_iabc("k,a,d,c") - t1("a,l") * g_ijab("l,k,c,d");

  return W_AkCd;
};

template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> compute_cs_ccsd_W_KlIc(
    LCAOFactoryBase<Tile, Policy>& lcao_factory,
    const TA::DistArray<Tile, Policy>& t1, bool df) {
  std::wstring postfix = df ? L"[df]" : L"";

  TA::DistArray<Tile, Policy> W_KlIc;

  auto g_ijab = lcao_factory.compute(L"<i j|G|a b>" + postfix);
  auto g_ijka = lcao_factory.compute(L"<i j|G|k a>" + postfix);

  W_KlIc("k,l,i,c") = g_ijka("k,l,i,c") + t1("d,i") * g_ijab("k,l,d,c");

  return W_KlIc;
};

template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> compute_cs_ccsd_W_KaIj(
    LCAOFactoryBase<Tile, Policy>& lcao_factory,
    const TA::DistArray<Tile, Policy>& t1,
    const TA::DistArray<Tile, Policy>& t2,
    const TA::DistArray<Tile, Policy>& tau,
    const TA::DistArray<Tile, Policy>& FIA,
    const TA::DistArray<Tile, Policy>& W_KlIj, bool df) {
  std::wstring postfix = df ? L"[df]" : L"";

  TA::DistArray<Tile, Policy> W_KaIj;

  auto g_ijab = lcao_factory.compute(L"<i j|G|a b>" + postfix);
  auto g_iajb = lcao_factory.compute(L"<i a|G|j b>" + postfix);
  auto g_iabc = lcao_factory.compute(L"<i a|G|b c>" + postfix);
  auto g_ijka = lcao_factory.compute(L"<i j|G|k a>" + postfix);

  TA::DistArray<Tile, Policy> t22;
  t22("a,b,i,j") = 2.0 * t2("a,b,i,j") - t2("a,b,j,i");

  W_KaIj("k,a,i,j") =

      g_ijka("i,j,k,a") + FIA("k,c") * t2("c,a,i,j") -

      t1("a,l") * W_KlIj("k,l,i,j") + g_iabc("k,a,c,d") * tau("c,d,i,j") +

      g_ijka("k,l,i,c") * t22("a,c,j,l") -

      g_ijka("l,k,i,c") * t2("a,c,j,l") - g_ijka("l,k,j,c") * t2("c,a,i,l") +

      t1("c,i") * (g_ijab("j,k,a,c") + t22("a,d,j,l") * g_ijab("k,l,c,d") -
                   t2("a,d,j,l") * g_ijab("k,l,d,c")) -

      t1("c,j") * (-g_iajb("i,c,k,a") + t2("a,d,l,i") * g_ijab("k,l,d,c"));

  return W_KaIj;
};

template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> compute_cs_ccsd_W_AbCi(
    LCAOFactoryBase<Tile, Policy>& lcao_factory,
    gaussian::AOFactoryBase<Tile, Policy>& ao_factory,
    const TA::DistArray<Tile, Policy>& t1,
    const TA::DistArray<Tile, Policy>& t2,
    const TA::DistArray<Tile, Policy>& tau,
    const TA::DistArray<Tile, Policy>& FIA,
    const TA::DistArray<Tile, Policy>& W_AbCd, bool df) {
  std::wstring postfix = df ? L"[df]" : L"";

  auto g_ijab = lcao_factory.compute(L"<i j|G|a b>" + postfix);
  auto g_iajb = lcao_factory.compute(L"<i a|G|j b>" + postfix);
  auto g_iabc = lcao_factory.compute(L"<i a|G|b c>" + postfix);
  auto g_ijka = lcao_factory.compute(L"<i j|G|k a>" + postfix);

  TA::DistArray<Tile, Policy> t22;
  t22("a,b,i,j") = 2.0 * t2("a,b,i,j") - t2("a,b,j,i");

  TA::DistArray<Tile, Policy> W_AbCi;

  W_AbCi("a,b,c,i") =

      g_iabc("i,a,b,c")

      - FIA("k,c") * t2("a,b,k,i")

      + g_ijka("l,k,i,c") * tau("a,b,k,l")

      - g_iabc("k,b,c,d") * t2("a,d,k,i") + g_iabc("k,a,d,c") * t22("b,d,i,k") -
      g_iabc("k,a,c,d") * t2("b,d,i,k")

      -
      t1("a,k") * (g_ijab("i,k,b,c") + (t22("d,b,l,i") * g_ijab("k,l,c,d") -
                                        t2("d,b,l,i") * g_ijab("k,l,d,c"))) +
      t1("b,k") * (-g_iajb("k,a,i,c") + t2("a,d,l,i") * g_ijab("l,k,c,d"));

  // if W_AbCd term is computed and stored
  if (W_AbCd.is_initialized()) {
    W_AbCi("a,b,c,i") += t1("d,i") * W_AbCd("a,b,c,d");
  }
  else{
    W_AbCi("a,b,c,i") += -t1("d,i") * g_iabc("k,a,d,c") * t1("b,k") -
        t1("d,i") * g_iabc("k,b,c,d") * t1("a,k") +
        t1("d,i") * g_ijab("k,l,c,d") * tau("a,b,k,l");

    // integral direct term
    auto direct_integral = ao_factory.compute_direct(L"(μ ν| G|κ λ)");
    auto Ca = lcao_factory.orbital_registry().retrieve(OrbitalIndex(L"a"));

    W_AbCi("a,b,c,i") += (t1("d,i")*Ca("s,d")*direct_integral("p,q,r,s"))*Ca("r,b")*Ca("q,c")*Ca("p,a");

  }

  return W_AbCi;
};

/*
 *
 */
template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> compute_cs_ccsd_W_AbCd(
    LCAOFactoryBase<Tile, Policy>& lcao_factory,
    const TA::DistArray<Tile, Policy>& t1,
    const TA::DistArray<Tile, Policy>& tau, bool df) {
  std::wstring postfix = df ? L"[df]" : L"";

  TA::DistArray<Tile, Policy> W_AbCd;

  auto g_abcd = lcao_factory.compute(L"<a b|G|c d>" + postfix);
  auto g_ijab = lcao_factory.compute(L"<i j|G|a b>" + postfix);
  auto g_iabc = lcao_factory.compute(L"<i a|G|b c>" + postfix);

  W_AbCd("a,b,c,d") = g_abcd("a,b,c,d") - t1("b,k") * g_iabc("k,a,d,c") -
                      t1("a,k") * g_iabc("k,b,c,d") +
                      tau("a,b,k,l") * g_ijab("k,l,c,d");

  return W_AbCd;
};

/**
 * compute closed shell one particle matrix of effective Hamiltonian of CCSD
 *
 * @return std::tuple of FAB, FIA, FIJ
 *
 */
template <typename Tile, typename Policy>
std::tuple<TA::DistArray<Tile, Policy>, TA::DistArray<Tile, Policy>,
           TA::DistArray<Tile, Policy>>
compute_cs_ccsd_F(LCAOFactoryBase<Tile, Policy>& lcao_factory,
                  gaussian::AOFactoryBase<Tile, Policy>& ao_factory,
                  const TA::DistArray<Tile, Policy>& t1,
                  const TA::DistArray<Tile, Policy>& tau, bool df) {
  std::wstring postfix = df ? L"[df]" : L"";

  auto g_ijab_bar = lcao_factory.compute(L"<i j|G|a b>" + postfix);
  g_ijab_bar("i,j,a,b") = 2.0 * g_ijab_bar("i,j,a,b") - g_ijab_bar("j,i,a,b");

  auto f_ia = lcao_factory.compute(L"<i|F|a>" + postfix);

  // compute FIA
  TA::DistArray<Tile, Policy> FIA;
  { FIA("i,a") = f_ia("i,a") + g_ijab_bar("i,j,a,b") * t1("b,j"); }

  TA::DistArray<Tile, Policy> FIJ;
  {
    FIJ = lcao_factory.compute(L"<i|F|j>" + postfix);
    auto g_ijka = lcao_factory.compute(L"<i j|G|k a>" + postfix);
    FIJ("i,j") += t1("c,l") * (2.0 * g_ijka("i,l,j,c") - g_ijka("l,i,j,c")) +
                  tau("c,d,j,l") * g_ijab_bar("i,l,c,d");
  }

  TA::DistArray<Tile, Policy> FAB;
  {
    FAB = lcao_factory.compute(L"<a|F|b>" + postfix);
    FAB("a,b") -=
        t1("a,m") * f_ia("m,b") + tau("a,d,k,l") * g_ijab_bar("k,l,b,d");

    if (df) {
      //     refactorize with density fitting
      auto Xia = lcao_factory.compute(L"( Λ |G|i a)");
      auto Xab = lcao_factory.compute(L"( Λ |G|a b)");
      auto X = ao_factory.compute(L"( Κ |G| Λ)[inv]");

      FAB("a,b") += 2.0 * t1("d,k") * Xia("X,k,d") * X("X,Y") * Xab("Y,a,b") -
                    t1("d,k") * Xab("X,a,d") * X("X,Y") * Xia("Y,k,b");

    } else {
      auto g_iabc = lcao_factory.compute(L"<i a|G|b c>" + postfix);
      FAB("a,b") += t1("d,k") * (2.0 * g_iabc("k,a,d,b") - g_iabc("k,a,b,d"));
    }
  }

  return std::make_tuple(FAB, FIA, FIJ);
}

}  // namespace lcao
}  // namespace cc
}  // namespace mpqc

#endif  // SRC_MPQC_CHEMISTRY_QC_LCAO_CC_CCSD_HBAR_H_
