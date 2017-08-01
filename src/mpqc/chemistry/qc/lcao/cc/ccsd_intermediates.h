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
 * @brief This file contains intermediates used in CCSD
 *
 */


/**
 *  Structure to hold effective Hamiltonian Components
 */
template <typename Tile, typename Policy>
struct Intermediates {

  using TArray = TA::DistArray<Tile,Policy>;

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
 * @param g_ijab_bar   2<i j|G|a b> - <i j|G|b a>
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
 * @param g_ijab_bar  2<i j|G|a b> - <i j|G|b a>
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
 * @param g_ijab_bar  2<i j|G|a b> - <i j|G|b a>
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
    const TA::DistArray <Tile, Policy> &F_oo,
    const TA::DistArray <Tile, Policy> &f_ov,
    const TA::DistArray <Tile, Policy> &g_ijka_bar,
    const TA::DistArray <Tile, Policy> &t1) {
  TA::DistArray<Tile, Policy> F_OO;
  F_OO("k,i") =
      F_oo("k,i") + f_ov("k,c") * t1("c,i") + g_ijka_bar("k,l,i,c") * t1("c,l");
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
                      g_ijka("k,l,j,c") * t1("c,i") +
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
      g_iajb("i,a,c,b") * t1("c,j") - W_ooov("i,k,j,b") * t1("a,k");
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
 * @param g_ijab_bar  2<i j|G|a b> - <i j|G|b a>
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
  W_ovvo1("i,a,b,j") = g_iajb("j,a,i,b") +
                       g_ijab_bar("i,k,b,c") * t2("a,c,j,k") -
                       g_ijab("i,k,b,c") * t2("a,c,k,j");
  return W_ovov1;
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
  W_ovvo2("i,a,b,j") =  g_iabc("i,a,b,c")*t1("c,j") - W_ooov("k,i,j,b")*t1("a,k");
  return W_ovov2;
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
  return W_ovov1;
};

template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> compute_cs_ccsd_W_ovoo(

)
{

};

}  // end of namespace cc
}  // end of namespace lcao
}  // end of namespace mpqc

#endif  // SRC_MPQC_CHEMISTRY_QC_LCAO_CC_CCSD_INTERMEDIATES_H_
