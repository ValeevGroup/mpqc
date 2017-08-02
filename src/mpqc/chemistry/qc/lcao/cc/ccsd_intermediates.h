//
// Created by Chong Peng on 7/21/17.
//

#ifndef SRC_MPQC_CHEMISTRY_QC_LCAO_CC_CCSD_INTERMEDIATES_H_
#define SRC_MPQC_CHEMISTRY_QC_LCAO_CC_CCSD_INTERMEDIATES_H_

#include "mpqc/chemistry/qc/lcao/factory/lcao_factory.h"
#include <tiledarray.h>

namespace mpqc {
namespace lcao {
namespace cc {

/*
 * @brief This file contains intermediates used in spin-adapted closed-shell
 * CCSD
 *
 */

/**
 *  Structure to hold effective Hamiltonian Components
 */
template <typename Tile, typename Policy>
struct Intermediates {
  using TArray = TA::DistArray<Tile, Policy>;

  Intermediates() = default;
  ~Intermediates() = default;

  TArray FIA;
  TArray FIJ;
  TArray FAB;

  TArray Wijkl;
  TArray Wijka;
  TArray Wiajk;
  TArray Wiabj;
  TArray Wiajb;
  TArray Wijab;
  TArray Wabci;
  TArray Waibc;
  TArray Wabcd;
};

///
/// @section CCSD intermediates in solving T1 T2 amplitudes
///

/**
 * @param f_oo  <i|F|j>
 * @param g_ijab_bar   2<i j|G|a b> - <j i|G|a b>
 * @param tau   T2("a,b,i,j") + T1("a,i")*T1("b,j")
 */
template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> compute_cs_ccsd_F_oo(
    const TA::DistArray<Tile, Policy>& f_oo,
    const TA::DistArray<Tile, Policy>& g_ijab_bar,
    const TA::DistArray<Tile, Policy>& tau) {
  TA::DistArray<Tile, Policy> F_oo;
  F_oo("k,l") = f_oo("k,l") + g_ijab_bar("k,i,a,b") * tau("a,b,l,i");
  return F_oo;
};

/**
 *
 * @param f_vv   <a|F|b>
 * @param g_ijab_bar 2<i j|G|a b> - <j i|G|a b>
 * @param tau   T2("a,b,i,j") + T1("a,i")*T1("b,j")
 */
template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> compute_cs_ccsd_F_vv(
    const TA::DistArray<Tile, Policy>& f_vv,
    const TA::DistArray<Tile, Policy>& g_ijab_bar,
    const TA::DistArray<Tile, Policy>& tau) {
  TA::DistArray<Tile, Policy> F_vv;
  F_vv("a,c") = f_vv("a,c") + g_ijab_bar("k,l,c,d") * tau("a,d,k,l");
  return F_vv;
};

///
///  @section CCSD one-particle effective Hamiltonian
///

/**
 *
 * @param f_ov  <i|F|a>
 * @param g_ijab_bar  2<i j|G|a b> - <j i|G|a b>
 * @param t1  CCSD T1 amplitude T1("a,i")
 */
template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> compute_cs_ccsd_F_ov(
    const TA::DistArray<Tile, Policy>& f_ov,
    const TA::DistArray<Tile, Policy>& g_ijab_bar,
    const TA::DistArray<Tile, Policy>& t1) {
  TA::DistArray<Tile, Policy> F_ov;
  F_ov("k,c") = f_ov("k,c") + g_ijab_bar("k,l,c,d") * t1("d,l");
  return F_ov;
};

/**
 *
 * @param F_oo  CCSD intermediate F_oo
 * @param f_ov  <i|F|a>
 * @param g_ijka_bar  2 <i j|G|k a> - <j i|G|k a>
 * @param t1 CCSD T1 amplitude T1("a,i")
 * @return
 */
template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> compute_cs_ccsd_F_OO(
    const TA::DistArray<Tile, Policy>& F_oo,
    const TA::DistArray<Tile, Policy>& f_ov,
    const TA::DistArray<Tile, Policy>& g_ijka_bar,
    const TA::DistArray<Tile, Policy>& t1) {
  TA::DistArray<Tile, Policy> F_OO;
  F_OO("k,i") =
      F_oo("k,i") + f_ov("k,c") * t1("c,i") + g_ijka_bar("k,l,i,c") * t1("c,l");
  return F_OO;
};
/**
 *
 * @param f_oo  <i|F|j>
 * @param f_ov <i|F|a>
 * @param g_ijab_bar   2<i j|G|a b> - <j i|G|a b>
 * @param g_ijka_bar 2 <i j|G|k a> - <j i|G|k a>
 * @param tau  T2("a,b,i,j") + T1("a,i")*T1("b,j")
 * @param t1   CCSD T1 amplitude T1("a,i")
 * @return
 */
template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> compute_cs_ccsd_F_OO(
    const TA::DistArray<Tile, Policy>& f_oo,
    const TA::DistArray<Tile, Policy>& f_ov,
    const TA::DistArray<Tile, Policy>& g_ijab_bar,
    const TA::DistArray<Tile, Policy>& g_ijka_bar,
    const TA::DistArray<Tile, Policy>& tau,
    const TA::DistArray<Tile, Policy>& t1) {
  TA::DistArray<Tile, Policy> F_OO;
  F_OO("k,i") = f_oo("k,i") + g_ijab_bar("k,l,a,b") * tau("a,b,i,l") +
                f_ov("k,c") * t1("c,i") + g_ijka_bar("k,l,i,c") * t1("c,l");
  return F_OO;
};

/**
 *
 * @param F_vv  CCSD intermediate F_vv
 * @param f_ov  <i|F|a>
 * @param g_iabc_bar  2 <i a|G|b c> - <i a|G|c b>
 * @param t1 CCSD T1 amplitude T1("a,i")
 * @return
 */
template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> compute_cs_ccsd_F_VV(
    const TA::DistArray<Tile, Policy>& F_vv,
    const TA::DistArray<Tile, Policy>& f_ov,
    const TA::DistArray<Tile, Policy>& g_iabc_bar,
    const TA::DistArray<Tile, Policy>& t1) {
  TA::DistArray<Tile, Policy> F_VV;
  F_VV("a,c") =
      F_vv("a,c") + f_ov("k,c") * t1("a,k") + g_iabc_bar("k,a,d,c") * t1("d,k");
  return F_VV;
};

/**
 *
 * @param f_vv  <a|F|b>
 * @param f_ov <i|F|a>
 * @param g_ijab_bar   2<i j|G|a b> - <j i|G|a b>
 * @param g_iabc_bar  2 <i a|G|b c> - <i a|G|c b>
 * @param tau  T2("a,b,i,j") + T1("a,i")*T1("b,j")
 * @param t1   CCSD T1 amplitude T1("a,i")
 * @return
 */
template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> compute_cs_ccsd_F_VV(
    const TA::DistArray<Tile, Policy>& f_vv,
    const TA::DistArray<Tile, Policy>& f_ov,
    const TA::DistArray<Tile, Policy>& g_ijab_bar,
    const TA::DistArray<Tile, Policy>& g_iabc_bar,
    const TA::DistArray<Tile, Policy>& tau,
    const TA::DistArray<Tile, Policy>& t1) {
  TA::DistArray<Tile, Policy> F_VV;
  F_VV("a,c") = f_vv("a,c") - g_ijab_bar("k,l,c,d") * tau("a,d,k,l") -
                f_ov("k,c") * t1("a,k") + g_iabc_bar("k,a,d,c") * t1("d,k");
  return F_VV;
};

///
///  @section CCSD two-particle effective Hamiltonian
///

/**
 * @param g_ijkl  <i j|G|k l>
 * @param g_ijka  <i j|G|k a>
 * @param g_ijab  <i j|G|a b>
 * @param tau     T2("a,b,i,j") + T1("a,i")*T1("b,j")
 * @param t1 CCSD T1 amplitude T1("a,i")
 * @return W_oooo intermediate
 */
template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> compute_cs_ccsd_W_oooo(
    const TA::DistArray<Tile, Policy>& g_ijkl,
    const TA::DistArray<Tile, Policy>& g_ijka,
    const TA::DistArray<Tile, Policy>& g_ijab,
    const TA::DistArray<Tile, Policy>& tau,
    const TA::DistArray<Tile, Policy>& t1) {
  TA::DistArray<Tile, Policy> W_oooo;

  W_oooo("k,l,i,j") = g_ijkl("k,l,i,j") + g_ijka("k,l,i,c") * t1("c,j") +
                      g_ijka("l,k,j,c") * t1("c,i") +
                      g_ijab("k,l,c,d") * tau("c,d,i,j");

  return W_oooo;
};

/**
 *
 * @param g_ijka  <i j|G|k a>
 * @param g_ijab  <i j|G|a b>
 * @param t1  CCSD T1 amplitude T1("a,i")
 * @return  W_ooov intermediate
 */
template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> compute_cs_ccsd_W_ooov(
    const TA::DistArray<Tile, Policy>& g_ijka,
    const TA::DistArray<Tile, Policy>& g_ijab,
    const TA::DistArray<Tile, Policy>& t1) {
  TA::DistArray<Tile, Policy> W_ooov;

  W_ooov("i,j,k,a") = g_ijka("i,j,k,a") + g_ijab("i,j,b,a") * t1("b,k");

  return W_ooov;
};

/**
 *
 * @param g_iajb  <i a|G|j b>
 * @param g_ijab  <i j|G|a b>
 * @param t2    CCSD T2 amplitude T2("a,b,i,j")
 * @return  first term in W_ovov intermediate
 */
template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> compute_cs_ccsd_W_ovov1(
    const TA::DistArray<Tile, Policy>& g_iajb,
    const TA::DistArray<Tile, Policy>& g_ijab,
    const TA::DistArray<Tile, Policy>& t2) {
  TA::DistArray<Tile, Policy> W_ovov1;
  W_ovov1("i,a,j,b") = g_iajb("i,a,j,b") - g_ijab("i,k,c,b") * t2("c,a,j,k");
  return W_ovov1;
};

/**
 *
 * @param W_ooov ccsd intermediate W_ooov
 * @param g_iabc  <i a|G|b c>
 * @param t1  CCSD T1 amplitude T1("a,i")
 * @return
 */
template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> compute_cs_ccsd_W_ovov2(
    const TA::DistArray<Tile, Policy>& W_ooov,
    const TA::DistArray<Tile, Policy>& g_iabc,
    const TA::DistArray<Tile, Policy>& t1) {
  TA::DistArray<Tile, Policy> W_ovov2;
  W_ovov2("i,a,j,b") =
      g_iabc("i,a,c,b") * t1("c,j") - W_ooov("i,k,j,b") * t1("a,k");
  return W_ovov2;
};

/**
 *
 * @param W_ooov  ccsd intermediate W_ooov
 * @param g_iabc   <i a|G|b c>
 * @param g_ijab   <i j|G|a b>
 * @param g_iajb   <i a|G|j b>
 * @param t2    CCSD T2 amplitude T2("a,b,i,j")
 * @param t1    CCSD T1 amplitude T1("a,i")
 * @return     W_ovov intermediate
 */
template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> compute_cs_ccsd_W_ovov(
    const TA::DistArray<Tile, Policy>& W_ooov,
    const TA::DistArray<Tile, Policy>& g_iabc,
    const TA::DistArray<Tile, Policy>& g_ijab,
    const TA::DistArray<Tile, Policy>& g_iajb,
    const TA::DistArray<Tile, Policy>& t2,
    const TA::DistArray<Tile, Policy>& t1) {
  TA::DistArray<Tile, Policy> W_ovov1 =
      compute_cs_ccsd_W_ovov1(g_iajb, g_ijab, t2);
  TA::DistArray<Tile, Policy> W_ovov2 =
      compute_cs_ccsd_W_ovov2(W_ooov, g_iabc, t1);
  W_ovov1("i,a,j,b") += W_ovov2("i,a,j,b");
  return W_ovov1;
};

/**
 *
 * @param g_ijab_bar  2<i j|G|a b> - <j i|G|a b>
 * @param g_ijab  <i j|G|a b>
 * @param g_iajb  <i a|G|j b>
 * @param t2  CCSD T2 amplitude T2("a,b,i,j")
 * @return first term in W_ovvo intermediate
 */
template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> compute_cs_ccsd_W_ovvo1(
    const TA::DistArray<Tile, Policy>& g_ijab_bar,
    const TA::DistArray<Tile, Policy>& g_ijab,
    const TA::DistArray<Tile, Policy>& g_iajb,
    const TA::DistArray<Tile, Policy>& t2) {
  TA::DistArray<Tile, Policy> W_ovvo1;
  W_ovvo1("i,a,b,j") = g_ijab("j,i,a,b") +
                       g_ijab_bar("i,k,b,c") * t2("a,c,j,k") -
                       g_ijab("i,k,b,c") * t2("a,c,k,j");
  return W_ovvo1;
};

/**
 *
 * @param W_ooov ccsd intermediate W_ooov
 * @param g_iabc  <i a|G|b c>
 * @param t1  CCSD T1 amplitude T1("a,i")
 * @return second term in W_ovvo intermediate
 */
template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> compute_cs_ccsd_W_ovvo2(
    const TA::DistArray<Tile, Policy>& W_ooov,
    const TA::DistArray<Tile, Policy>& g_iabc,
    const TA::DistArray<Tile, Policy>& t1) {
  TA::DistArray<Tile, Policy> W_ovvo2;
  W_ovvo2("i,a,b,j") =
      g_iabc("i,a,b,c") * t1("c,j") - W_ooov("k,i,j,b") * t1("a,k");
  return W_ovvo2;
};

/**
 *
 * @param W_ooov  ccsd intermediate W_ooov
 * @param g_iabc   <i a|G|b c>
 * @param g_iajb   <i a|G|j b>
 * @param g_ijab   <i j|G|a b>
 * @param t2    CCSD T2 amplitude T2("a,b,i,j")
 * @param t1    CCSD T1 amplitude T1("a,i")
 * @return     W_ovvo intermediate
 */
template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> compute_cs_ccsd_W_ovvo(
    const TA::DistArray<Tile, Policy>& W_ooov,
    const TA::DistArray<Tile, Policy>& g_iabc,
    const TA::DistArray<Tile, Policy>& g_ijab_bar,
    const TA::DistArray<Tile, Policy>& g_ijab,
    const TA::DistArray<Tile, Policy>& g_iajb,
    const TA::DistArray<Tile, Policy>& t2,
    const TA::DistArray<Tile, Policy>& t1) {
  TA::DistArray<Tile, Policy> W_ovvo1 =
      compute_cs_ccsd_W_ovvo1(g_ijab_bar, g_ijab, g_iajb, t2);
  TA::DistArray<Tile, Policy> W_ovvo2 =
      compute_cs_ccsd_W_ovvo2(W_ooov, g_iabc, t1);
  W_ovvo1("i,a,b,j") += W_ovvo2("i,a,b,j");
  return W_ovvo1;
};

template <typename Tile, typename Policy>
cc::Intermediates<Tile, Policy> compute_intermediates(
    LCAOFactoryBase<Tile, Policy>& lcao_factory,
    gaussian::AOFactoryBase<Tile, Policy>& ao_factory,
    const TA::DistArray<Tile, Policy>& t2,
    const TA::DistArray<Tile, Policy>& t1, bool df, std::string option) {
  cc::Intermediates<Tile, Policy> imds;

  std::wstring postfix = df ? L"[df]" : L"";

  auto f_ia = lcao_factory.compute(L"<i|F|a>" + postfix);
  auto f_ij = lcao_factory.compute(L"<i|F|j>" + postfix);
  auto f_ab = lcao_factory.compute(L"<a|F|b>" + postfix);

  auto g_ijkl = lcao_factory.compute(L"<i j|G|k l>" + postfix);
  auto g_ijab = lcao_factory.compute(L"<i j|G|a b>" + postfix);
  auto g_iajb = lcao_factory.compute(L"<i a|G|j b>" + postfix);
  auto g_iabc = lcao_factory.compute(L"<i a|G|b c>" + postfix);
  auto g_ijka = lcao_factory.compute(L"<i j|G|k a>" + postfix);

  TA::DistArray<Tile, Policy> tau;
  tau("a,b,i,j") = t2("a,b,i,j") + (t1("a,i") * t1("b,j"));

  TA::DistArray<Tile, Policy> g_ijab_bar;
  g_ijab_bar("i,j,a,b") = 2.0 * g_ijab("i,j,a,b") - g_ijab("j,i,a,b");

  // one body
  {
    TA::DistArray<Tile, Policy> g_ijka_bar;
    g_ijka_bar("i,j,k,a") = 2.0 * g_ijka("i,j,k,a") - g_ijka("j,i,k,a");
//
    imds.FIJ =
        compute_cs_ccsd_F_OO(f_ij, f_ia, g_ijab_bar, g_ijka_bar, tau, t1);
//
    imds.FIA = compute_cs_ccsd_F_ov(f_ia, g_ijab_bar, t1);

    TA::DistArray<Tile, Policy> g_iabc_bar;
    g_iabc_bar("k,a,d,b") = 2.0 * g_iabc("k,a,d,b") - g_iabc("k,a,b,d");
//
    imds.FAB =
        compute_cs_ccsd_F_VV(f_ab, f_ia, g_ijab_bar, g_iabc_bar, tau, t1);

//    std::tie(imds.FAB, imds.FIA, imds.FIJ) = compute_cs_ccsd_F(lcao_factory,ao_factory,t1,tau,df);
  }

  // two body
  {
    imds.Wijab = g_ijab;

    imds.Wijka = compute_cs_ccsd_W_ooov(g_ijka, g_ijab, t1);

    if (option == "ip") {
      imds.Wiajb =
          compute_cs_ccsd_W_ovov(imds.Wijka, g_iabc, g_ijab, g_iajb, t2, t1);

      imds.Wiabj = compute_cs_ccsd_W_ovvo(imds.Wijka, g_iabc, g_ijab_bar,
                                          g_ijab, g_iajb, t2, t1);

      imds.Wijkl = compute_cs_ccsd_W_oooo(g_ijkl, g_ijka, g_ijab, tau, t1);

      imds.Wiajk = compute_cs_ccsd_W_KaIj(lcao_factory, t1, t2, tau, imds.FIA,
                                          imds.Wijkl, df);
    } else if (option == "ea") {
      imds.Wiajb =
          compute_cs_ccsd_W_ovov(imds.Wijka, g_iabc, g_ijab, g_iajb, t2, t1);

      imds.Wiabj = compute_cs_ccsd_W_ovvo(imds.Wijka, g_iabc, g_ijab_bar,
                                          g_ijab, g_iajb, t2, t1);

      imds.Wabcd = compute_cs_ccsd_W_AbCd(lcao_factory, t1, tau, df);

      imds.Waibc = compute_cs_ccsd_W_AkCd(lcao_factory, t1, df);

      imds.Wabci = compute_cs_ccsd_W_AbCi(lcao_factory, ao_factory, t1, t2, tau,
                                          imds.FIA, imds.Wabcd, df);

    } else {
      throw InputError("Wrong Option in compute_intermediates", __FILE__,
                       __LINE__);
    }
  }

  return imds;
};

}  // end of namespace cc
}  // end of namespace lcao
}  // end of namespace mpqc

#endif  // SRC_MPQC_CHEMISTRY_QC_LCAO_CC_CCSD_INTERMEDIATES_H_
