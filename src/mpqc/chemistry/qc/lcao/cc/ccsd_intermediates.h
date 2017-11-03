//
// Created by Chong Peng on 7/21/17.
//

#ifndef SRC_MPQC_CHEMISTRY_QC_LCAO_CC_CCSD_INTERMEDIATES_H_
#define SRC_MPQC_CHEMISTRY_QC_LCAO_CC_CCSD_INTERMEDIATES_H_

#include <tiledarray.h>
#include "mpqc/chemistry/qc/lcao/factory/lcao_factory.h"

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
template <typename Array>
struct Intermediates {
  Intermediates() = default;
  ~Intermediates() = default;

  Array FIA;
  Array FIJ;
  Array FAB;

  Array Wijkl;
  Array Wijka;
  Array Wiajk;
  Array Wiabj;
  Array Wiajb;
  Array Wijab;
  Array Wabci;
  Array Waibc;
  Array Wabcd;

  Array Xab;
  Array Xia;
  Array Xij;
};

///
/// @section CCSD intermediates in solving T1 T2 amplitudes
///

/**
 * @param f_oo  <i|F|j>
 * @param g_ijab_bar   2<i j|G|a b> - <j i|G|a b>
 * @param tau   T2("a,b,i,j") + T1("a,i")*T1("b,j")
 */
template <typename Array>
Array compute_cs_ccsd_F_oo(const Array& f_oo, const Array& g_ijab_bar,
                           const Array& tau) {
  Array F_oo;
  F_oo("k,i") = f_oo("k,i") + g_ijab_bar("k,l,a,b") * tau("a,b,i,l");
  return F_oo;
};

/**
 *
 * @param f_vv   <a|F|b>
 * @param g_ijab_bar 2<i j|G|a b> - <j i|G|a b>
 * @param tau   T2("a,b,i,j") + T1("a,i")*T1("b,j")
 */
template <typename Array>
Array compute_cs_ccsd_F_vv(const Array& f_vv, const Array& g_ijab_bar,
                           const Array& tau) {
  Array F_vv;
  F_vv("a,c") = f_vv("a,c") - g_ijab_bar("k,l,c,d") * tau("a,d,k,l");
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
template <typename Array>
Array compute_cs_ccsd_F_ov(const Array& f_ov, const Array& g_ijab_bar,
                           const Array& t1) {
  Array F_ov;
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
template <typename Array>
Array compute_cs_ccsd_F_OO(const Array& F_oo, const Array& f_ov,
                           const Array& g_ijka_bar, const Array& t1) {
  Array F_OO;
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
template <typename Array>
Array compute_cs_ccsd_F_OO(const Array& f_oo, const Array& f_ov,
                           const Array& g_ijab_bar, const Array& g_ijka_bar,
                           const Array& tau, const Array& t1) {
  Array F_OO;
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
template <typename Array>
Array compute_cs_ccsd_F_VV(const Array& F_vv, const Array& f_ov,
                           const Array& g_iabc_bar, const Array& t1) {
  Array F_VV;
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
template <typename Array>
Array compute_cs_ccsd_F_VV(const Array& f_vv, const Array& f_ov,
                           const Array& g_ijab_bar, const Array& g_iabc_bar,
                           const Array& tau, const Array& t1) {
  Array F_VV;
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
template <typename Array>
Array compute_cs_ccsd_W_oooo(const Array& g_ijkl, const Array& g_ijka,
                             const Array& g_ijab, const Array& tau,
                             const Array& t1) {
  Array W_oooo;

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
template <typename Array>
Array compute_cs_ccsd_W_ooov(const Array& g_ijka, const Array& g_ijab,
                             const Array& t1) {
  Array W_ooov;

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
template <typename Array>
Array compute_cs_ccsd_W_ovov1(const Array& g_iajb, const Array& g_ijab,
                              const Array& t2) {
  Array W_ovov1;
  W_ovov1("i,a,j,b") = g_iajb("i,a,j,b") - g_ijab("i,k,c,b") * t2("c,a,j,k");
  return W_ovov1;
};

/**
 *
 * @param W_ooov ccsd intermediate W_ooov
 * @param g_iabc  <i a|G|b c>
 * @param t1  CCSD T1 amplitude T1("a,i")
 * @return second term in W_ovov intermediate
 */
template <typename Array>
Array compute_cs_ccsd_W_ovov2(const Array& W_ooov, const Array& g_iabc,
                              const Array& t1) {
  Array W_ovov2;
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
template <typename Array>
Array compute_cs_ccsd_W_ovov(const Array& W_ooov, const Array& g_iabc,
                             const Array& g_ijab, const Array& g_iajb,
                             const Array& t2, const Array& t1) {
  Array W_ovov1 = compute_cs_ccsd_W_ovov1(g_iajb, g_ijab, t2);
  Array W_ovov2 = compute_cs_ccsd_W_ovov2(W_ooov, g_iabc, t1);
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
template <typename Array>
Array compute_cs_ccsd_W_ovvo1(const Array& g_ijab_bar, const Array& g_ijab,
                              const Array& g_iajb, const Array& t2) {
  Array W_ovvo1;
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
template <typename Array>
Array compute_cs_ccsd_W_ovvo2(const Array& W_ooov, const Array& g_iabc,
                              const Array& t1) {
  Array W_ovvo2;
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
template <typename Array>
Array compute_cs_ccsd_W_ovvo(const Array& W_ooov, const Array& g_iabc,
                             const Array& g_ijab_bar, const Array& g_ijab,
                             const Array& g_iajb, const Array& t2,
                             const Array& t1) {
  Array W_ovvo1 = compute_cs_ccsd_W_ovvo1(g_ijab_bar, g_ijab, g_iajb, t2);
  Array W_ovvo2 = compute_cs_ccsd_W_ovvo2(W_ooov, g_iabc, t1);
  W_ovvo1("i,a,b,j") += W_ovvo2("i,a,b,j");
  return W_ovvo1;
};

}  // end of namespace cc
}  // end of namespace lcao
}  // end of namespace mpqc

#endif  // SRC_MPQC_CHEMISTRY_QC_LCAO_CC_CCSD_INTERMEDIATES_H_
