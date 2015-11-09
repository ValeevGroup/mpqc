//
// sr_r12intermediates_VXB_diag.h
//
// Copyright (C) 2013 Edward Valeev
//
// Author: Edward Valeev <evaleev@vt.edu>
// Maintainer: EV
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#ifndef _mpqc_src_lib_chemistry_qc_mbptr12_srr12intermediatesVXBdiag_h
#define _mpqc_src_lib_chemistry_qc_mbptr12_srr12intermediatesVXBdiag_h

#include <tiledarray.h>
#include <cmath>

namespace sc {
  inline double get_element(const TA::Array<double, 4 >& array, const std::vector<std::size_t>& ele_idx)
  {
    const std::vector<std::size_t> tile_idx = array.trange().element_to_tile(ele_idx);
    return (array.find(tile_idx).get()[ele_idx]);
  }

  inline double get_element(const TA::Array<double, 2 >& array, const std::vector<std::size_t>& ele_idx)
  {
    const std::vector<std::size_t> tile_idx = array.trange().element_to_tile(ele_idx);
    return (array.find(tile_idx).get()[ele_idx]);
  }

  template <typename T>
  typename SingleReference_R12Intermediates<T>::TArray2
  SingleReference_R12Intermediates<T>::XaiAddToXam(const TA::Array<double, 2 >& Xam,
                                                   const TA::Array<double, 2 >& Xai) {
//    MPQC_ASSERT(Xam.size() == 1);
//
//    typename TArray2::value_type tile_X =
//    static_cast<const typename TArray2::value_type&>(Xam.find(0)).clone();
//    typename TArray2::value_type tile_Xai = Xai.find(0);
//    const std::size_t zero_cols = tile_X.range().extent()[1] - tile_Xai.range().extent()[1];
//
//    std::array<std::size_t, 2> i = {{0, 0}};
//    std::size_t ix = 0;
//    for(i[0] = tile_X.range().lobound()[0];
//      i[0] < tile_X.range().upbound()[0]; ++i[0]) {
//      for(i[1] = zero_cols + tile_X.range().lobound()[1];
//        i[1] < tile_X.range().upbound()[1]; ++i[1], ++ix) {
//        tile_X[i] +=  tile_Xai[ix];
//      }
//    }
//
//    TArray2 X(Xam.get_world(), Xam.trange());
//    X.set(0, tile_X);
    TArray2 X;
    X("a,m") = Xam("a,m") + Xai("a,i") * _2("<i|I|m>");
    return X;
  }

  // compute F^A'_C' R^Pk_A'B' R^Qk_C'B' (P & Q in alpha or beta space)
  template <typename T>
  typename SingleReference_R12Intermediates<T>::TArray2
  SingleReference_R12Intermediates<T>::BPk_Qk(const char* p, const char* q,
                                              const double C_0, const double C_1) {

    const double B_C1 = 0.5 * C_0 * C_0 + 1.5 * C_1 * C_1;
    const double B_C2 = 0.5 * C_0 * C_0 - 1.5 * C_1 * C_1;

    std::string rTr_pkql = std::string("<") + p + " k|rTr|" + q + " l>";
    std::string rTr_kpql = std::string("<k ") + p + "|rTr|" + q + " l>";

    std::string r2_phJk_ql = std::string("<") + p + "_hJ(p') k|r2|" + q + " l>";
    std::string r2_pkhJ_ql = std::string("<") + p + " k_hJ(p')|r2|" + q + " l>";
    //
    std::string r2_kphJ_ql = std::string("<k ") + p + "_hJ(p')|r2|" + q + " l>";
    std::string r2_khJp_ql = std::string("<k_hJ(p') ") + p + "|r2|" + q + " l>";
    //
    std::string r2_pk_qhJl = std::string("<") + p + " k|r2|" + q + "_hJ(p') l>";
    std::string r2_pk_qlhJ = std::string("<") + p + " k|r2|" + q + " l_hJ(p')>";
    //
    std::string r2_pk_lqhJ = std::string("<") + p + " k|r2|l " + q + "_hJ(p')>";
    std::string r2_pk_lhJq = std::string("<") + p + " k|r2|l_hJ(p') " + q + ">";

    std::string pk_PQ = std::string("<") + p + " k|r|p' q'>";
    std::string kp_PQ = std::string("<k ") + p + "|r|p' q'>";
    std::string  qk_PKQ = std::string("<") + q + " k|r|p'_K(r') q'>";
    std::string  kq_PKQ = std::string("<k ") + q + "|r|p'_K(r') q'>";

    std::string pk_Pn = std::string("<") + p + " k|r|p' n>";
    std::string kp_Pn = std::string("<k ") + p + "|r|p' n>";
    std::string  qk_PFn = std::string("<") + q + " k|r|p'_F(r') n>";
    std::string  kq_PFn = std::string("<k ") + q + "|r|p'_F(r') n>";

    std::string pk_mA = std::string("<") + p + " k|r|m a'>";
    std::string kp_mA = std::string("<k ") + p + "|r|m a'>";
    std::string  qk_mFA = std::string("<") + q + " k|r|m_F(n) a'>";
    std::string  kq_mFA = std::string("<k ") + q + "|r|m_F(n) a'>";

    std::string pk_pq = std::string("<") + p + " k|r|p b>";
    std::string kp_pq = std::string("<k ") + p + "|r|p b>";
    std::string  qk_pFb = std::string("<") + q + " k|r|p_F(r) b>";
    std::string  kq_pFb = std::string("<k ") + q + "|r|p_F(r) b>";

    std::string  qk_mFPA = std::string("<") + q + " k|r|m_F(p') a'>";
    std::string  kq_mFPA = std::string("<k ") + q + "|r|m_F(p') a'>";
    //
    std::string qk_mA = std::string("<") + q + " k|r|m a'>";
    std::string kq_mA = std::string("<k ") + q + "|r|m a'>";
    std::string  pk_mFPA = std::string("<") + p + " k|r|m_F(p') a'>";
    std::string  kp_mFPA = std::string("<k ") + p + "|r|m_F(p') a'>";

    std::string pk_Ab = std::string("<") + p + " k|r|a' b>";
    std::string kp_Ab = std::string("<k ") + p + "|r|a' b>";
    std::string  qk_AFb = std::string("<") + q + " k|r|a'_F(q) b>";
    std::string  kq_AFb = std::string("<k ") + q + "|r|a'_F(q) b>";
    //
    std::string qk_Ab = std::string("<") + q + " k|r|a' b>";
    std::string kq_Ab = std::string("<k ") + q + "|r|a' b>";
    std::string  pk_AFb = std::string("<") + p + " k|r|a'_F(q) b>";
    std::string  kp_AFb = std::string("<k ") + p + "|r|a'_F(q) b>";

    std::string pq =  p + std::string(",") + q;

    TArray2 B_pq;
    B_pq(pq.c_str()) =
                 //          diag
                 ( B_C1 * _4(rTr_pkql.c_str()) + B_C2 * _4(rTr_kpql.c_str())
                 //           Q
                 + 0.5 * (
                     B_C1 * (_4(r2_phJk_ql.c_str()) + _4(r2_pkhJ_ql.c_str()))
                   + B_C2 * (_4(r2_kphJ_ql.c_str()) + _4(r2_khJp_ql.c_str()))

                   + B_C1 * (_4(r2_pk_qhJl.c_str()) + _4(r2_pk_qlhJ.c_str()))
                   + B_C2 * (_4(r2_pk_lqhJ.c_str()) + _4(r2_pk_lhJq.c_str()))
                         )
                  ) * _2("<k|I|l>")
                 //           rKr_p'q'
                 - ( B_C1 * (_4(pk_PQ.c_str()) * _4(qk_PKQ.c_str())
                           + _4(kp_PQ.c_str()) * _4(kq_PKQ.c_str()))
                   + B_C2 * (_4(kp_PQ.c_str()) * _4(qk_PKQ.c_str())
                           + _4(pk_PQ.c_str()) * _4(kq_PKQ.c_str()))
                   )
                 //           rFr_p'n
                 - ( B_C1 * (_4(pk_Pn.c_str()) * _4(qk_PFn.c_str())
                           + _4(kp_Pn.c_str()) * _4(kq_PFn.c_str()))
                   + B_C2 * (_4(kp_Pn.c_str()) * _4(qk_PFn.c_str())
                           + _4(pk_Pn.c_str()) * _4(kq_PFn.c_str()))
                    )
                 //           rFr_mA
                 + ( B_C1 * (_4(pk_mA.c_str()) * _4(qk_mFA.c_str())
                           + _4(kp_mA.c_str()) * _4(kq_mFA.c_str()))
                   + B_C2 * (_4(kp_mA.c_str()) * _4(qk_mFA.c_str())
                           + _4(pk_mA.c_str()) * _4(kq_mFA.c_str()))
                   )
                 //           rFr_pb
                 - ( B_C1 * (_4(pk_pq.c_str()) * _4(qk_pFb.c_str())
                           + _4(kp_pq.c_str()) * _4(kq_pFb.c_str()))
                   + B_C2 * (_4(kp_pq.c_str()) * _4(qk_pFb.c_str())
                           + _4(pk_pq.c_str()) * _4(kq_pFb.c_str()))
                   )
                 //
                 - ( B_C1 * (_4(pk_mA.c_str()) * _4(qk_mFPA.c_str())
                           + _4(kp_mA.c_str()) * _4(kq_mFPA.c_str()))
                   + B_C2 * (_4(kp_mA.c_str()) * _4(qk_mFPA.c_str())
                           + _4(pk_mA.c_str()) * _4(kq_mFPA.c_str()))

                   + B_C1 * (_4(qk_mA.c_str()) * _4(pk_mFPA.c_str())
                           + _4(kq_mA.c_str()) * _4(kp_mFPA.c_str()))
                   + B_C2 * (_4(kq_mA.c_str()) * _4(pk_mFPA.c_str())
                           + _4(qk_mA.c_str()) * _4(kp_mFPA.c_str()))
                    )
                   //
                 - ( B_C1 * (_4(pk_Ab.c_str()) * _4(qk_AFb.c_str())
                           + _4(kp_Ab.c_str()) * _4(kq_AFb.c_str()))
                   + B_C2 * (_4(kp_Ab.c_str()) * _4(qk_AFb.c_str())
                           + _4(pk_Ab.c_str()) * _4(kq_AFb.c_str()))

                   + B_C1 * (_4(qk_Ab.c_str()) * _4(pk_AFb.c_str())
                           + _4(kq_Ab.c_str()) * _4(kp_AFb.c_str()))
                   + B_C2 * (_4(kq_Ab.c_str()) * _4(pk_AFb.c_str())
                           + _4(qk_Ab.c_str()) * _4(kp_AFb.c_str()))
                   );
    return B_pq;
  }

  template <typename T>
  typename SingleReference_R12Intermediates<T>::TArray4
  SingleReference_R12Intermediates<T>::Bpr_qs(const char* p, const char* q)
  {

    const double C_0 = 0.5;
    const double C_1 = 0.25;
    const double B_C1 = 0.5 * C_0 * C_0 + 1.5 * C_1 * C_1;
    const double B_C2 = 0.5 * C_0 * C_0 - 1.5 * C_1 * C_1;

    std::string rTr_pkql = std::string("<") + p + " k|rTr|" + q + " l>";
    std::string rTr_kpql = std::string("<k ") + p + "|rTr|" + q + " l>";

    std::string r2_phJk_ql = std::string("<") + p + "_hJ(p') k|r2|" + q + " l>";
    std::string r2_pkhJ_ql = std::string("<") + p + " k_hJ(p')|r2|" + q + " l>";
    //
    std::string r2_kphJ_ql = std::string("<k ") + p + "_hJ(p')|r2|" + q + " l>";
    std::string r2_khJp_ql = std::string("<k_hJ(p') ") + p + "|r2|" + q + " l>";
    //
    std::string r2_pk_qhJl = std::string("<") + p + " k|r2|" + q + "_hJ(p') l>";
    std::string r2_pk_qlhJ = std::string("<") + p + " k|r2|" + q + " l_hJ(p')>";
    //
    std::string r2_pk_lqhJ = std::string("<") + p + " k|r2|l " + q + "_hJ(p')>";
    std::string r2_pk_lhJq = std::string("<") + p + " k|r2|l_hJ(p') " + q + ">";

    std::string pk_PQ = std::string("<") + p + " k|r|p' q'>";
    std::string kp_PQ = std::string("<k ") + p + "|r|p' q'>";
    std::string  ql_PKQ = std::string("<") + q + " l|r|p'_K(r') q'>";
    std::string  lq_PKQ = std::string("<l ") + q + "|r|p'_K(r') q'>";

    std::string pk_Pn = std::string("<") + p + " k|r|p' n>";
    std::string kp_Pn = std::string("<k ") + p + "|r|p' n>";
    std::string  ql_PFn = std::string("<") + q + " l|r|p'_F(r') n>";
    std::string  lq_PFn = std::string("<l ") + q + "|r|p'_F(r') n>";

    std::string pk_mA = std::string("<") + p + " k|r|m a'>";
    std::string kp_mA = std::string("<k ") + p + "|r|m a'>";
    std::string  ql_mFA = std::string("<") + q + " l|r|m_F(n) a'>";
    std::string  lq_mFA = std::string("<l ") + q + "|r|m_F(n) a'>";

    std::string pk_pq = std::string("<") + p + " k|r|p b>";
    std::string kp_pq = std::string("<k ") + p + "|r|p b>";
    std::string  ql_pFb = std::string("<") + q + " l|r|p_F(r) b>";
    std::string  lq_pFb = std::string("<l ") + q + "|r|p_F(r) b>";

    std::string  ql_mFPA = std::string("<") + q + " l|r|m_F(p') a'>";
    std::string  lq_mFPA = std::string("<l ") + q + "|r|m_F(p') a'>";
    //
    std::string ql_mA = std::string("<") + q + " l|r|m a'>";
    std::string lq_mA = std::string("<l ") + q + "|r|m a'>";
    std::string  pk_mFPA = std::string("<") + p + " k|r|m_F(p') a'>";
    std::string  kp_mFPA = std::string("<k ") + p + "|r|m_F(p') a'>";

    std::string pk_Ab = std::string("<") + p + " k|r|a' b>";
    std::string kp_Ab = std::string("<k ") + p + "|r|a' b>";
    std::string  ql_AFb = std::string("<") + q + " l|r|a'_F(q) b>";
    std::string  lq_AFb = std::string("<l ") + q + "|r|a'_F(q) b>";
    //
    std::string ql_Ab = std::string("<") + q + " l|r|a' b>";
    std::string lq_Ab = std::string("<l ") + q + "|r|a' b>";
    std::string  pk_AFb = std::string("<") + p + " k|r|a'_F(q) b>";
    std::string  kp_AFb = std::string("<k ") + p + "|r|a'_F(q) b>";

    TArray4 Bpr_qs;
    Bpr_qs("p,r,q,s") =
                       //          diag
                         B_C1 * _4(rTr_pkql.c_str()) + B_C2 * _4(rTr_kpql.c_str())
                       //           Q
                       + 0.5 * (
                           B_C1 * (_4(r2_phJk_ql.c_str()) + _4(r2_pkhJ_ql.c_str()))
                         + B_C2 * (_4(r2_kphJ_ql.c_str()) + _4(r2_khJp_ql.c_str()))

                         + B_C1 * (_4(r2_pk_qhJl.c_str()) + _4(r2_pk_qlhJ.c_str()))
                         + B_C2 * (_4(r2_pk_lqhJ.c_str()) + _4(r2_pk_lhJq.c_str()))
                               )
                       //           rKr_p'q'
                       - ( B_C1 * (_4(pk_PQ.c_str()) * _4(ql_PKQ.c_str())
                                 + _4(kp_PQ.c_str()) * _4(lq_PKQ.c_str()))
                         + B_C2 * (_4(kp_PQ.c_str()) * _4(ql_PKQ.c_str())
                                 + _4(pk_PQ.c_str()) * _4(lq_PKQ.c_str()))
                         )
                       //           rFr_p'n
                       - ( B_C1 * (_4(pk_Pn.c_str()) * _4(ql_PFn.c_str())
                                 + _4(kp_Pn.c_str()) * _4(lq_PFn.c_str()))
                         + B_C2 * (_4(kp_Pn.c_str()) * _4(ql_PFn.c_str())
                                 + _4(pk_Pn.c_str()) * _4(lq_PFn.c_str()))
                          )
                       //           rFr_mA
                       + ( B_C1 * (_4(pk_mA.c_str()) * _4(ql_mFA.c_str())
                                 + _4(kp_mA.c_str()) * _4(lq_mFA.c_str()))
                         + B_C2 * (_4(kp_mA.c_str()) * _4(ql_mFA.c_str())
                                 + _4(pk_mA.c_str()) * _4(lq_mFA.c_str()))
                         )
                       //           rFr_pb
                       - ( B_C1 * (_4(pk_pq.c_str()) * _4(ql_pFb.c_str())
                                 + _4(kp_pq.c_str()) * _4(lq_pFb.c_str()))
                         + B_C2 * (_4(kp_pq.c_str()) * _4(ql_pFb.c_str())
                                 + _4(pk_pq.c_str()) * _4(lq_pFb.c_str()))
                         )
                       //
                       - ( B_C1 * (_4(pk_mA.c_str()) * _4(ql_mFPA.c_str())
                                 + _4(kp_mA.c_str()) * _4(lq_mFPA.c_str()))
                         + B_C2 * (_4(kp_mA.c_str()) * _4(ql_mFPA.c_str())
                                 + _4(pk_mA.c_str()) * _4(lq_mFPA.c_str()))

                         + B_C1 * (_4(ql_mA.c_str()) * _4(pk_mFPA.c_str())
                                 + _4(lq_mA.c_str()) * _4(kp_mFPA.c_str()))
                         + B_C2 * (_4(lq_mA.c_str()) * _4(pk_mFPA.c_str())
                                 + _4(ql_mA.c_str()) * _4(kp_mFPA.c_str()))
                          )
                      //
                      -  ( B_C1 * (_4(pk_Ab.c_str()) * _4(ql_AFb.c_str())
                                 + _4(kp_Ab.c_str()) * _4(lq_AFb.c_str()))
                         + B_C2 * (_4(kp_Ab.c_str()) * _4(ql_AFb.c_str())
                                 + _4(pk_Ab.c_str()) * _4(lq_AFb.c_str()))

                         + B_C1 * (_4(ql_Ab.c_str()) * _4(pk_AFb.c_str())
                                 + _4(lq_Ab.c_str()) * _4(kp_AFb.c_str()))
                         + B_C2 * (_4(lq_Ab.c_str()) * _4(pk_AFb.c_str())
                                 + _4(ql_Ab.c_str()) * _4(kp_AFb.c_str()))
                         );
    return Bpr_qs;
  }

  // compute V^PQ_RS = 1/2 R^PQ_A'B' g^A'B'_RS (all indices in alpha or beta space)
  // do not use p1, q1, a1', or n1 as argument
  template <typename T>
  typename SingleReference_R12Intermediates<T>::TArray4
  SingleReference_R12Intermediates<T>::VPQ_RS(const char* p, const char* q,
                                              const char* r, const char* s)
  {

    std::string gr_pqrs = std::string("<") + p + " " + q + "|gr|" + r + " " + s + ">";

    std::string rpq_pq = std::string("<") + p + " " + q + "|r|p1 q1>";
    std::string grs_pq = std::string("<p1 q1|g|") + r + " " + s + ">";

    std::string rpq_apn = std::string("<") + p + " " + q + "|r|a1' n1>";
    std::string grs_apn = std::string("<a1' n1|g|") + r + " " + s + ">";

    std::string rpq_nap = std::string("<") + p + " " + q + "|r|n1 a1'>";
    std::string grs_nap = std::string("<n1 a1'|g|") + r + " " + s + ">";

    std::string pqrs =  p + std::string(",") + q + std::string(",")
                      + r + std::string(",") + s;

    TArray4 Vpq_rs;
    Vpq_rs(pqrs.c_str()) =  _4(gr_pqrs.c_str())
                          - _4(rpq_pq.c_str()) * _4(grs_pq.c_str())
                          - _4(rpq_apn.c_str()) * _4(grs_apn.c_str())
                          - _4(rpq_nap.c_str()) * _4(grs_nap.c_str())
                          ;
    return Vpq_rs;
  }

  // compute V^Pq_Rs = 1/2 bar{R}^Pq_A'B' bar{g}^A'B'_Rs (P, R in alpha or beta space)
  // do not use p1, q1, a1', or n1 as argument
  template <typename T>
  typename SingleReference_R12Intermediates<T>::TArray4
  SingleReference_R12Intermediates<T>::VPq_Rs(const char* p, const char* q,
                                              const char* r, const char* s,
                                              const double C_0, const double C_1)
  {
    const double V_C1 = (0.5 * C_0 + 1.5 * C_1);
    const double V_C2 = (0.5 * C_0 - 1.5 * C_1);

    std::string gr_pqrs = std::string("<") + p + " " + q + "|gr|" + r + " " + s + ">";
    std::string gr_qprs = std::string("<") + q + " " + p + "|gr|" + r + " " + s + ">";

    std::string rpq_pq = std::string("<") + p + " " + q + "|r|p1 q1>";
    std::string rqp_pq = std::string("<") + q + " " + p + "|r|p1 q1>";
    std::string grs_pq = std::string("<p1 q1|g|") + r + " " + s + ">";

    std::string rpq_apn = std::string("<") + p + " " + q + "|r|a1' n1>";
    std::string rqp_apn = std::string("<") + q + " " + p + "|r|a1' n1>";
    std::string grs_apn = std::string("<a1' n1|g|") + r + " " + s + ">";

    std::string rpq_nap = std::string("<") + p + " " + q + "|r|n1 a1'>";
    std::string rqp_nap = std::string("<") + q + " " + p + "|r|n1 a1'>";
    std::string grs_nap = std::string("<n1 a1'|g|") + r + " " + s + ">";

    std::string pqrs =  p + std::string(",") + q + std::string(",")
                      + r + std::string(",") + s;

    TArray4 Vpq_rs;
    Vpq_rs(pqrs.c_str()) =  V_C1 * _4(gr_pqrs.c_str()) + V_C2 * _4(gr_qprs.c_str())
                          - (V_C1 * _4(rpq_pq.c_str()) + V_C2 * _4(rqp_pq.c_str()))
                            * _4(grs_pq.c_str())
                          - (V_C1 * _4(rpq_apn.c_str()) + V_C2 * _4(rqp_apn.c_str()))
                            * _4(grs_apn.c_str())
                          - (V_C1 * _4(rpq_nap.c_str()) + V_C2 * _4(rqp_nap.c_str()))
                            * _4(grs_nap.c_str())
                          ;
    return Vpq_rs;
  }

  // compute V^Rk_Sk = 1/2 bar{R}^Rk_A'B' bar{g}^A'B'_Sk
  // can not use k1 l2 n1 a1' for arguments
  template <typename T>
  typename SingleReference_R12Intermediates<T>::TArray2
  SingleReference_R12Intermediates<T>::VRk_Sk(const char* r, const char* s,
                                              const double C_0, const double C_1)
  {
    const double V_C1 = (0.5 * C_0 + 1.5 * C_1);
    const double V_C2 = (0.5 * C_0 - 1.5 * C_1);

    std::string gr_rksk = std::string("<") + r + " k1|gr|" + s + " l1>";
    std::string gr_krsk = std::string("<k1 ") + r + "|gr|" + s + " l1>";

    std::string r_rkpq = std::string("<") + r + " k1|r|p1 q1>";
    std::string r_krpq = std::string("<k1 ") + r + "|r|p1 q1>";
    std::string g_skpq = std::string("<p1 q1|g|") + s + " k1>";

    std::string r_rkapn = std::string("<") + r + " k1|r|a1' n1>";
    std::string r_krapn = std::string("<k1 ") + r + "|r|a1' n1>";
    std::string g_skapn = std::string("<a1' n1|g|") + s + " k1>";

    std::string r_rknap = std::string("<") + r + " k1|r|n1 a1'>";
    std::string r_krnap = std::string("<k1 ") + r + "|r|n1 a1'>";
    std::string g_sknap = std::string("<n1 a1'|g|") + s + " k1>";

    std::string rs = r + std::string(",") + s;
    TArray2 Vrk_sk;
    Vrk_sk(rs.c_str()) = (V_C1 * _4(gr_rksk.c_str()) + V_C2 * _4(gr_krsk.c_str()))
                          * _2("<l1|I|k1")
                        - (V_C1 * _4(r_rkpq.c_str()) + V_C2 * _4(r_krpq.c_str()))
                          * _4(g_skpq.c_str())
                        - (V_C1 * _4(r_rkapn.c_str()) + V_C2 * _4(r_krapn.c_str()))
                          * _4(g_skapn.c_str())
                        - (V_C1 * _4(r_rknap.c_str()) + V_C2 * _4(r_krnap.c_str()))
                          * _4(g_sknap.c_str())
                        ;
    return Vrk_sk;
  }

  template <typename T>
  std::pair<typename SingleReference_R12Intermediates<T>::TArray2,
            typename SingleReference_R12Intermediates<T>::TArray2>
  SingleReference_R12Intermediates<T>::V_diag() {

    TArray2 V_ij_ij_cabs; V_ij_ij_cabs("i,j") = dotket(ij_xy("<i j|g|m a'>"), ij_xy("<i j|r|m a'>"));
    // don't need this in closed shell!
    //TArray2 V_ij_ij_cabs1 = dotket(ij_xy("<i j|g|a' m>"), ij_xy("<i j|r|a' m>"));

    TArray2 V_ij_ij; V_ij_ij("i,j") = take(ij_xy("<i j|gr|p q>"), ij) - dotket(ij_xy("<i j|g|p q>"), ij_xy("<i j|r|p q>"))
                                      - V_ij_ij_cabs("i,j") - V_ij_ij_cabs("j,i");

    TArray2 V_ij_ji_cabs = dotket(ij_xy("<i j|g|m a'>"), ij_xy("<i j|r|m a'>"), true);
    TArray2 V_ij_ji = take(ij_xy("<i j|gr|p q>"), ji) - dotket(ij_xy("<i j|g|p q>"), ij_xy("<i j|r|p q>"), true)
                      - V_ij_ji_cabs("i,j") - V_ij_ji_cabs("j,i");

    return std::make_pair(V_ij_ij,V_ij_ji);

  }

  template <typename T>
  std::pair<typename SingleReference_R12Intermediates<T>::TArray2,
            typename SingleReference_R12Intermediates<T>::TArray2>
  SingleReference_R12Intermediates<T>::X_diag() {

    TArray2 X_ij_ij_cabs = dotket(ij_xy("<i j|r|m a'>"), ij_xy("<i j|r|m a'>"));
    // don't need this in closed shell!
    //TArray2 X_ij_ij_cabs1 = dotket(ij_xy("<i j|r|a' m>"), ij_xy("<i j|r|a' m>"));

    TArray2 X_ij_ij = take(ij_xy("<i j|r2|p q>"), ij) - dotket(ij_xy("<i j|r|p q>"), ij_xy("<i j|r|p q>"))
                      - X_ij_ij_cabs("i,j") - X_ij_ij_cabs("j,i");

    TArray2 X_ij_ji_cabs = dotket(ij_xy("<i j|r|m a'>"), ij_xy("<i j|r|m a'>"), true);
    TArray2 X_ij_ji = take(ij_xy("<i j|r2|p q>"), ji) - dotket(ij_xy("<i j|r|p q>"), ij_xy("<i j|r|p q>"), true)
                      - X_ij_ji_cabs("i,j") - X_ij_ji_cabs("j,i");

    return std::make_pair(X_ij_ij,X_ij_ji);

  }

  template <typename T>
  std::pair<typename SingleReference_R12Intermediates<T>::TArray2,
            typename SingleReference_R12Intermediates<T>::TArray2>
  SingleReference_R12Intermediates<T>::B_diag() {

    /// this is incomplete at the moment
    TArray2 B_ij_ij = take(ij_xy("<i j|rTr|p q>"), ij) + take(ij_xy("<i_hJ(p') j|r2|p q>"), ij)
                      - dotket(ij_xy("<i j|r|p' q'>"), ij_xy("<i j|r|p' q'_K(r')>"));

    TArray2 B_ij_ji = take(ij_xy("<i j|rTr|p q>"), ji) + take(ij_xy("<i_hJ(p') j|r2|p q>"), ji);

    return std::make_pair(B_ij_ij,B_ij_ji);

  }

  namespace detail {

    /** makes a pseudo-3-index (4-index, but first index is dummy) self-energy denominator:
        (se_den)_{xyz} = 1/ ( E + <x|Y|x> - <y|X1|y> - <z|X2|z>);
      */
    template<typename T>
    struct selfenergy_denom {
        typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> EigenMatrixX;
        selfenergy_denom(const T& E,
                         const EigenMatrixX& Y_mat,
                         const EigenMatrixX& X1_mat,
                         const EigenMatrixX& X2_mat) :
            E_(E), Y_mat_(Y_mat), X1_mat_(X1_mat), X2_mat_(X2_mat) {
        }
        template<typename Index> T operator()(const Index& i) {
          return 1.0 / ( E_ + Y_mat_(i[1], i[1]) - X1_mat_(i[2], i[2]) - X2_mat_(i[3], i[3]) );
        }

      private:
        T E_;
        EigenMatrixX Y_mat_;
        EigenMatrixX X1_mat_;
        EigenMatrixX X2_mat_;
    };
  }; // namespace sc::detail

  template <typename T>
  void
  SingleReference_R12Intermediates<T>::gf2_r12(int orbital) {

    enum CalcType {
      IP, EA
    };

    MPQC_ASSERT(orbital != 0);
    const CalcType type = (orbital < 0) ? IP : EA;

    const std::string X_label = (type == IP) ? "i1" : "a1"; // active orbital label
    const Ref<OrbitalSpace>& X_space = (type == IP) ? this->r12world()->refwfn()->occ_act()
                                                    : this->r12world()->refwfn()->uocc_act();

    TArray2 fij = xy("<j|F|j>");
    TArray2 fab = xy("<a|F|b>");

    // use the diagonal element of the Fock matrix as the guess (only HOMO/LUMO supported);
    const double E = (type == IP) ? X_space->evals()(X_space->rank() + orbital) : X_space->evals()(orbital - 1);

    // compute Delta_ijaE =  1 / (- E - <a|F|a> + <i|F|i> + <j|F|j>)
    typedef detail::selfenergy_denom<double> sedenom_eval_type;
    sedenom_eval_type Delta_jab_gen(E,
                                    TA::array_to_eigen(fij),
                                    TA::array_to_eigen(fab),
                                    TA::array_to_eigen(fab) );
    sedenom_eval_type Delta_ajk_gen(E,
                                    TA::array_to_eigen(fab),
                                    TA::array_to_eigen(fij),
                                    TA::array_to_eigen(fij) );

    std::ostringstream Xjab; Xjab << "<" << X_label << " j|g|a b>";
    TArray4d g_Xjab = ijxy(Xjab.str());
    typedef TA::Array<T, 4, LazyTensor<T, 4, sedenom_eval_type > > TArray4dLazy;
    TArray4dLazy Delta_Xjab(g_Xjab.get_world(), g_Xjab.trange());

    // construct local tiles
    for(auto t = Delta_Xjab.trange().tiles().begin();
        t != Delta_Xjab.trange().tiles().end(); ++t)
      if (Delta_Xjab.is_local(*t)) {
        std::array<std::size_t, 4> index;
        std::copy(t->begin(), t->end(), index.begin());
        madness::Future < typename TArray4dLazy::value_type >
          tile((LazyTensor<T, 4, sedenom_eval_type >(&Delta_Xjab, index, &Delta_jab_gen)
              ));

        // Insert the tile into the array
        Delta_Xjab.set(*t, tile);
      }

    std::ostringstream Xajk; Xajk << "<" << X_label << " a|g|j k>";
    TArray4d g_Xajk = ijxy(Xajk.str());
    typedef TA::Array<T, 4, LazyTensor<T, 4, sedenom_eval_type > > TArray4dLazy;
    TArray4dLazy Delta_Xajk(g_Xajk.get_world(), g_Xajk.trange());

    // construct local tiles
    for(auto t = Delta_Xajk.trange().tiles().begin();
        t != Delta_Xajk.trange().tiles().end(); ++t)
      if (Delta_Xajk.is_local(*t)) {
        std::array<std::size_t, 4> index;
        std::copy(t->begin(), t->end(), index.begin());
        madness::Future < typename TArray4dLazy::value_type >
          tile((LazyTensor<T, 4, sedenom_eval_type >(&Delta_Xajk, index, &Delta_ajk_gen)
              ));

        // Insert the tile into the array
        Delta_Xajk.set(*t, tile);
      }

    TArray4 dg_Xjab; dg_Xjab("Y,j,a,b") = Delta_Xjab("Y,j,a,b") * (4 * g_Xjab("Y,j,a,b") - 2 * g_Xjab("Y,j,b,a"));
    TArray4 dg_Xajk; dg_Xajk("Y,a,j,m") = Delta_Xajk("Y,a,j,m") * (4 * g_Xajk("Y,a,j,m") - 2 * g_Xajk("Y,a,m,j"));

    // second-order self-energy, Eq. 11.16 in Dickoff-VanNeck, with spin integrated out
    // \sigma_{X}^{Y} =
//    TArray2 sigma2; sigma2("X,Y") = (1./2) * ( g_Xjab("X,j,a,b") * dg_Xjab("Y,j,a,b") +
//                                               g_Xajk("X,a,j,m") * dg_Xajk("Y,a,j,m") );
    TArray2 sigma2_hpp; sigma2_hpp("X,Y") = (1./2) * g_Xjab("X,j,a,b") * dg_Xjab("Y,j,a,b");
    TArray2 sigma2_phh; sigma2_phh("X,Y") = (1./2) * g_Xajk("X,a,j,m") * dg_Xajk("Y,a,j,m");

    ExEnv::out0() << "Fock:" << std::endl << fij << std::endl;
    ExEnv::out0() << "SE2_hpp(e=" << E << "):" << std::endl << sigma2_hpp << std::endl;
    ExEnv::out0() << "SE2_phh(e=" << E << "):" << std::endl << sigma2_phh << std::endl;
  }

  namespace detail {
    /** this functor helps to implement conjugate gradient CABS singles solver
     */
    template<typename T>
    struct _CABS_singles_h0t1 {

        typedef TiledArray::Array<T, 2> Array;

        /**
         * @param h0_AB allvirt/allvirt Fock operator
         * @param h0_ij occ/occ Fock operator
         */
        _CABS_singles_h0t1(const Array& h0_AB, const Array& h0_ij) :
            H0_AB(h0_AB), H0_IJ(h0_ij) {
        }

        const Array& H0_AB;
        const Array& H0_IJ;

        /**
         * @param[in] T1 t_i^A
         * @param[out] R1 R_i^A
         */
        void operator()(const Array& T1, Array& R1) {
          R1 = T1("i,b") * H0_AB("b,a") - H0_IJ("i,j") * T1("j,a");
        }
    };

    /** this functor helps to implement orbital response
     */
    template<typename T>
    struct _OrbResponse {

        typedef TiledArray::Array<T, 2> Array2;
        typedef TiledArray::Array<T, 4> Array4;

        /**
         * @param f_AB virt/virt Fock operator
         * @param f_ij occ/occ Fock operator
         * @param g_ij_ab <ij|ab>
         * @param g_ia_jb <ia|jb>
         */
        _OrbResponse(const Array2& f_AB, const Array2& f_ij,
                     const Array4& g_ij_ab, const Array4& g_ia_jb) :
            F_AB(f_AB), F_IJ(f_ij), G_IJ_AB(g_ij_ab), G_IA_JB(g_ia_jb) {
        }

        const Array2& F_AB;
        const Array2& F_IJ;
        const Array4& G_IJ_AB;
        const Array4& G_IA_JB;

        /**
         * @param[in] kappa kappa_i^a
         * @param[out] residual residual_i^a
         */
        void operator()(const Array2& kappa, Array2& residual) {
          residual =  kappa("i,b") * F_AB("b,a") - F_IJ("i,j") * kappa("j,a")
          + 4.0 * G_IJ_AB("i,j,a,b") *  kappa("j,b")
                   -       G_IJ_AB("i,j,b,a") *  kappa("j,b")
                   -       G_IA_JB("i,a,j,b") *  kappa("j,b");
        }
    };

    /// makes a diagonal 2-index preconditioner: pc_x^y = -1/ ( <x|O1|x> - <y|O2|y> )
    template <typename T>
    struct diag_precond2 {
        typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> EigenMatrixX;
        diag_precond2(const EigenMatrixX& O1_mat,
                      const EigenMatrixX& O2_mat) :
                          O1_mat_(O1_mat), O2_mat_(O2_mat) {
        }
        template <typename Index> T operator()(const Index& i) {
          return 1.0 / (- O1_mat_(i[0], i[0]) + O2_mat_(i[1], i[1]));
        }

      private:
        EigenMatrixX O1_mat_;
        EigenMatrixX O2_mat_;
    };

    /// makes a diagonal 4-index preconditioner: pc_xy^zw = -1/ ( <x|O1|x> + <y|O2|y> - <z|O3|z> - <w|O4|w> )
    template <typename T>
    struct diag_precond4 {
        typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> EigenMatrixX;
        diag_precond4(const EigenMatrixX& O1_mat,
                      const EigenMatrixX& O2_mat,
                      const EigenMatrixX& O3_mat,
                      const EigenMatrixX& O4_mat) :
                          O1_mat_(O1_mat), O2_mat_(O2_mat),
                          O3_mat_(O3_mat), O4_mat_(O4_mat) {
        }
        template <typename Index> T operator()(const Index& i) {
          return 1.0 / (- O1_mat_(i[0], i[0]) - O2_mat_(i[1], i[1]) + O3_mat_(i[2], i[2]) + O4_mat_(i[3], i[3]));
        }

      private:
        EigenMatrixX O1_mat_;
        EigenMatrixX O2_mat_;
        EigenMatrixX O3_mat_;
        EigenMatrixX O4_mat_;
    };

    template<typename T>
    struct Orbital_relaxation_Abjai {

        typedef TiledArray::Array<T, 2> Array2;
        typedef TiledArray::Array<T, 4> Array4;

        /**
         * @param A_bjai
         */
        Orbital_relaxation_Abjai(const Array4& a_bjai) :
          A_bjai(a_bjai) {
        }

        const Array4& A_bjai;

        /**
         * @param[in] K_bj
         * @param[out] R1 R_i^a
         */
        void operator()(const Array2& K_bj, Array2& R1) {
          R1("a,i") = K_bj("b,j") * A_bjai("b,j,a,i");
        }
    };

    // e_ij = (e_i + e_j)
    template <typename T>
    struct e_ij {
        typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> EigenMatrixX;
        e_ij(const EigenMatrixX& O1_mat, const EigenMatrixX& O2_mat):
             O1_mat_(O1_mat), O2_mat_(O2_mat) {}

        template <typename Index> T operator()(const Index& i) {
          return (O1_mat_(i[0], i[0]) + O2_mat_(i[1], i[1]));
        }

      private:
        EigenMatrixX O1_mat_;
        EigenMatrixX O2_mat_;
    };

  } // namespace sc::detail

  // Xam contribution from CABS Singles
  template <typename T>
  typename SingleReference_R12Intermediates<T>::TArray2
  SingleReference_R12Intermediates<T>::Xam_CabsSingles(const TArray2& TmA,
                                                       const TArray2& Tma) {
    // density from CABS Singles contribution
    TArray2 D_e2_mn, D_e2_AB;
    // D^m_n =  t^m_A' * t^A'_n
    D_e2_mn("m,n")= TmA("m,A'") * TmA("n,A'");
    // D^A'_B' = t^A'_m * t^m_B'
    D_e2_AB("A',B'") = TmA("m,A'") * TmA("m,B'");

    // compute intermediate to reduce memory consumption
    TArray4 g_aAmB;
    g_aAmB("a,A',m,B'") = 2.0 * _4("<a A'|g|m B'>") - _4("<a A'|g|B' m>");
    TArray4::wait_for_lazy_cleanup(g_aAmB.get_world());

    TArray4d g_aAmn = ijxy("<a A'|g|m n>");
    TArray4d g_ammn = ijxy("<a m|g|m1 n>");

    TArray2 Xam_E2;
    Xam_E2("a,m") = 2.0 * (
                          - _2("<A'|F|a>") * TmA("m,A'")
                          + Tma("n,a") * _2("<m|F|n>")
                          + (TmA("n,A'") * Tma("n,a")) * _2("<m|F|A'>")
                          //
                          - ( 2.0 * _4("<a n|g|m A'>") - _4("<a n|g|A' m>")
                            + 2.0 * g_aAmn("a,A',m,n") - g_aAmn("a,A',n,m")
                            ) * TmA("n,A'")
                          //
                          - g_aAmB("a,A',m,B'") * D_e2_AB("A',B'")
                          //
                          + (2.0 * g_ammn("a,n1,m,n2") - g_ammn("a,n1,n2,m"))
                            * D_e2_mn("n1,n2")
                          );
    return Xam_E2;
  }

  // Xam contribution from MP2 correlation
  template <typename T>
  typename SingleReference_R12Intermediates<T>::TArray2
  SingleReference_R12Intermediates<T>::Xam_mp2(const TArray4& T2_ijab,
                                               const TArray2& Dij,
                                               const TArray2& Dab) {
    TArray2 Xai_mp2, Xam_mp2;
    Xai_mp2("a,i") = 2.0 * (
                           // derivatives of 1/4 g^kl_cd T^cd_kl
                           // g^al_cd T^cd_il
                           - _4("<a l|g|c d>")
                             * (2.0 * T2_ijab("i,l,c,d") - T2_ijab("l,i,c,d"))
                           );
    Xam_mp2("a,m") = 2.0 * (
                           // derivatives of 1/4 g^kl_cd T^cd_kl
                           // g_mc^kl T^ac_kl
                            (2.0 * T2_ijab("k,l,a,c") - T2_ijab("k,l,c,a"))
                            * _4("<m c|g|k l>")

                           // derivatives of 1/2 F^c_d D^d_c - 1/2 F^l_k D^k_l
                           // g^ac_md D^d_c
                           - (2.0 * _4("<a c|g|m d>") - _4("<a c|g|d m>"))
                             * Dab("d,c")
                           // g^al_mk D^k_l
                           + (2.0 * _4("<a l|g|m k>") - _4("<a l|g|k m>"))
                             * Dij("k,l")
                           );

    // test code: print out I intermediates
//    TArray2 Iam, Iia, Ima;
//    Iam("a,m") = - (2.0 * T2_ijab("k,l,a,c") - T2_ijab("k,l,c,a"))
//                   * _4("<m c|g|k l>")
//                   ;
//    Iia("a,i") = _4("<a l|g|c d>")
//                 * (2.0 * T2_ijab("i,l,c,d") - T2_ijab("l,i,c,d"));
//    Ima("a,m") =   (2.0 * _4("<a c|g|m d>") - _4("<a c|g|d m>"))
//                   * Dab("d,c")
//                 - (2.0 * _4("<a l|g|m k>") - _4("<a l|g|k m>"))
//                   * Dij("k,l")
//                 ;
//    TArray2 Iima = XaiAddToXam(Ima, Iia);
//    std::cout << "MP2 Iam: " << std::endl << Iam << std::endl;
//    std::cout << "MP2 Ima: " << std::endl << Iima << std::endl;

    return XaiAddToXam(Xam_mp2, Xai_mp2);
  }

  // Xai contribution from MP2 F12 coupling part
  template <typename T>
  typename SingleReference_R12Intermediates<T>::TArray2
  SingleReference_R12Intermediates<T>::Xam_Cmp2f12(const double C_0, const double C_1,
                                                   const TArray4& T2_ijab, const TArray4& A_ijab,
                                                   const TArray2& Dij, const TArray2& Dab,
                                                   const TArray2& RT_apb) {
    const double R_C1 = (0.5 * C_0 + 1.5 * C_1);
    const double R_C2 = (0.5 * C_0 - 1.5 * C_1);
    const double RR_C1 = 0.5 * C_0 * C_0 + 1.5 * C_1 * C_1;
    const double RR_C2 = 0.5 * C_0 * C_0 - 1.5 * C_1 * C_1;

    TArray4 T2pA_temp;
    T2pA_temp("i,j,a,b") =  R_C1 * T2_ijab("i,j,a,b") + R_C2 * T2_ijab("j,i,a,b")
                          + RR_C1 * A_ijab("i,j,a,b") + RR_C2 * A_ijab("j,i,a,b");

    TArray2 Xai_mp2f12, Xam_mp2f12;
    // Xai
    Xai_mp2f12("a,i") =  // derivatives of 1/4 g^kl_cd \tilde{T}^cd_kl
                         // g^al_cd \tilde{T}^cd_il
                       -  _4("<a l|g|c d>")
                          * (R_C1 * A_ijab("i,l,c,d") + R_C2 * A_ijab("l,i,c,d"));
    Xai_mp2f12("a,i") =  Xai_mp2f12("a,i")
                         // F^a'_b R^al_a'c \tilde{T}^bc_il
                         // 1st part
                       - _4("<a l|r|b_F(a') c>") * T2pA_temp("i,l,b,c");
    Xai_mp2f12("a,i") =  Xai_mp2f12("a,i")
                         // 2nd part
                       - _4("<a l|r|c b_F(a')>") * T2pA_temp("i,l,c,b");

    // Xam
    Xam_mp2f12("a,m") =  // derivatives of 1/4 g^kl_cd \tilde{T}^cd_kl
                         // g_ic^kl \tilde{T}^Ac_kl
                         (R_C1 * A_ijab("k,l,a,c") + R_C2 * A_ijab("k,l,c,a"))
                         * _4("<m c|g|k l>")

                         // derivatives of 1/2 F^c_d D^d_c - 1/2 F^l_k D^k_l
                         // g^ac_md D^d_c
                       - (2.0 * _4("<a c|g|m d>") - _4("<a c|g|d m>"))
                         * Dab("d,c")
                         // g^al_mk D^k_l
                       + (2.0 * _4("<a l|g|m k>") - _4("<a l|g|k m>"))
                         * Dij("k,l");
    Xam_mp2f12("a,m") =  Xam_mp2f12("a,m")
                         // derivatives of 1/2 F^a'_b R^kl_a'c \tilde{T}^bc_kl
                         //   1/2 F^a'_m R^kl_a'b \tilde{T}^ab_kl
                         // & 1/2 F^a'_b R^kl_ma' \tilde{T}^ab_kl
                       + T2pA_temp("k,l,a,b")
                         * (_4("<k l|r|m_F(a') b>") + _4("<k l|r|m b_F(a')>"));
    TArray4 g_temp;
    {
    TArray4d gab_map = ijxy("<a b|g|m a'>");
    g_temp("a,a',m,b") =  2.0 * _4("<a a'|g|m b>") - _4("<a a'|g|b m>")
                        + 2.0 * gab_map("a,b,m,a'") - gab_map("b,a,m,a'");
    }
    Xam_mp2f12("a,m") =  Xam_mp2f12("a,m")
                         //   g^aa'_mb 1/2 R^kl_a'c \tilde{T}^bc_kl
                         // & g^ab_ma' 1/2 R^kl_a'c \tilde{T}^bc_kl
                       - g_temp("a,a',m,b") * RT_apb("a',b");

    return XaiAddToXam(Xam_mp2f12, Xai_mp2f12);
  }

  // compute F12 density from X and B terms
  template <typename T>
  void SingleReference_R12Intermediates<T>::compute_Df12_XB(const double C_0, const double C_1,
                                                            TArray2& D_f12_ij, TArray2& D_f12_ab,
                                                            TArray2& D_f12_apbp, TArray2& D_f12_apb) {

    const double RR_C1 = 0.5 * C_0 * C_0 + 1.5 * C_1 * C_1;
    const double RR_C2 = 0.5 * C_0 * C_0 - 1.5 * C_1 * C_1;

    // Dij = 1/2 R^ik_A'B' R^A'B'_kl (A': all virtual)
    {
    TArray4d r2_ijkl = ijxy("<i j|r2|k l>");
    TArray4d r_ijpq = ijxy("<i j|r|p q>");
    TArray4d r_ijapn = ijxy("<i j|r|a' n>");
    D_f12_ij("i,j") =  (RR_C1 * r2_ijkl("i,k,j,l") + RR_C2 * r2_ijkl("k,i,j,l"))
                       * _2("<k|I|l>")
                     - (RR_C1 * r_ijpq("i,k,p,q") + RR_C2 * r_ijpq("k,i,p,q"))
                       * r_ijpq("j,k,p,q")
                     - (RR_C1 * r_ijapn("i,k,a',n") + RR_C2 * r_ijapn("k,i,a',n"))
                       * r_ijapn("j,k,a',n")
                     - (RR_C1 * r_ijapn("k,i,a',n") + RR_C2 * r_ijapn("i,k,a',n"))
                       * r_ijapn("k,j,a',n");
    }
    TArray4d::wait_for_lazy_cleanup(D_f12_ij.get_world());

    // DA'B' = 1/2 R^A'C'_kl R^kl_B'C' (A': all virtual)
    TArray4d r_acpkl = ijxy("<a c'|r|k l>");
    {
    TArray4d r_apcpkl = ijxy("<a' c'|r|k l>");
    D_f12_apb("a',b") = (RR_C1 * r_apcpkl("a',c',k,l") + RR_C2 * r_apcpkl("a',c',l,k"))
                        * r_acpkl("b,c',k,l");

    D_f12_apbp("a',b'") =  (RR_C1 * r_acpkl("c,a',l,k") + RR_C2 * r_acpkl("c,a',k,l"))
                           * r_acpkl("c,b',l,k")
                         + (RR_C1 * r_apcpkl("a',c',k,l") + RR_C2 * r_apcpkl("a',c',l,k"))
                           * r_apcpkl("b',c',k,l");
//    D_f12_apbp("a',b'") =  (RR_C1 * _4("<a' p'|r|k l>") + RR_C2 * _4("<a' p'|r|l k>"))
//                           * _4("<b' p'|r|k l>")
//                         - (RR_C1 * _4("<a' m|r|k l>") + RR_C2 * _4("<a' m|r|l k>"))
//                           * _4("<b' m|r|k l>");
    }
    TArray4d::wait_for_lazy_cleanup(D_f12_ij.get_world());

    D_f12_ab("a,b") = (RR_C1 * r_acpkl("a,c',k,l") + RR_C2 * r_acpkl("a,c',l,k"))
                      * r_acpkl("b,c',k,l");
  }

  // compute X and B density contribution to Xam
  template <typename T>
  typename SingleReference_R12Intermediates<T>::TArray2
  SingleReference_R12Intermediates<T>::Xam_Df12_XB(const double C_0, const double C_1,
                                                   const TArray2& Df12_ij, const TArray2& Df12_ab,
                                                   const TArray2& Df12_apbp, const TArray2& Df12_apb) {
    TArray2 gdf12_am;
    gdf12_am("a,m") =  // X related contribution
                       (2.0 * _4("<a k|g|m l>") - _4("<a k|g|l m>"))
                       * Df12_ij("k,l");

    TArray4d g_abmc = ijxy("<a b|g|m c>");
    gdf12_am("a,m") =  gdf12_am("a,m")
                       // B related contribution
                       // 1st part
                     - (2.0 * g_abmc("a,b,m,c") - g_abmc("b,a,m,c"))
                       * Df12_ab("b,c");
    {
    TArray4 g_temp;
    g_temp("a,b',m,c'") = 2.0 * _4("<a b'|g|m c'>") - _4("<a b'|g|c' m>");
    TArray4d::wait_for_lazy_cleanup(g_temp.get_world());

    gdf12_am("a,m") =  gdf12_am("a,m")
                       // 2nd part
                     - g_temp("a,b',m,c'") * Df12_apbp("b',c'");
    }
    {
    TArray4 g_temp;
    g_temp("a,b',m,c") = 2.0 * _4("<a b'|g|m c>") - _4("<a b'|g|c m>");
    gdf12_am("a,m") =  gdf12_am("a,m")
                       // 3rd part
                     - g_temp("a,b',m,c") * Df12_apb("b',c");
    }
    TArray4d g_abmcp = ijxy("<a b|g|m c'>");
    gdf12_am("a,m") =  gdf12_am("a,m")
                       // 4th part
                     - (2.0 * g_abmcp("a,b,m,c'") - g_abmcp("b,a,m,c'"))
                       * Df12_apb("c',b");
    return gdf12_am;
  }

  // Xam contribution from F12 V part
  template <typename T>
  typename SingleReference_R12Intermediates<T>::TArray2
  SingleReference_R12Intermediates<T>::Xam_V(const double C_0, const double C_1) {

    const double R_C1 = (0.5 * C_0 + 1.5 * C_1);
    const double R_C2 = (0.5 * C_0 - 1.5 * C_1);

    TArray2 Xai_V, Xam_V;
    {
    TArray4d gr_ak_ij = ijxy("<a k|gr|i j>");
    TArray4d rpq_kl = ijxy("<p q|r|k l>");
    TArray4d rapn_kl = ijxy("<a' n|r|k l>");
    TArray4d rpq_ak = ijxy("<p q|r|a k>");

    TArray4d gij_apm = ijxy("<i j|g|a' m>");
    TArray4d rapn_ak = ijxy("<a' n|r|a k>");
    TArray4d rapn_ka = ijxy("<a' n|r|k a>");

    TArray2 Ikl = xy("<k|I|l>");

    Xai_V("a,i") =  // 1/2 R^ik_AC g^AC_ak
                  - ( (R_C1 * gr_ak_ij("a,k,i,l") + R_C2 * gr_ak_ij("a,k,l,i"))
                      * Ikl("k,l")
                    - _4("<a k|g|p q>")
                       * ( R_C1 * rpq_kl("p,q,i,k") + R_C2 * rpq_kl("p,q,k,i") )
                    - _4("<a k|g|a' n>")
                       * ( R_C1 * rapn_kl("a',n,i,k") + R_C2 * rapn_kl("a',n,k,i"))
                    - _4("<k a|g|a' n>")
                       * ( R_C1 * rapn_kl("a',n,k,i") + R_C2 * rapn_kl("a',n,i,k"))
                    )
                    // 1/2 R^ak_AC g^AC_ik
                  - ( (R_C1 * gr_ak_ij("a,l,i,k") + R_C2 * gr_ak_ij("a,l,k,i"))
                      * Ikl("k,l")
                    - _4("<i k|g|p q>")
                      * ( R_C1 * rpq_ak("p,q,a,k") + R_C2 * rpq_ak("q,p,a,k"))
                    - gij_apm("i,k,a',n")
                      * ( R_C1 * rapn_ak("a',n,a,k") + R_C2 * rapn_ka("a',n,k,a"))
                    - gij_apm("k,i,a',n")
                      * ( R_C1 * rapn_ka("a',n,k,a") + R_C2 * rapn_ak("a',n,a,k"))
                    );
    }
    TArray4d::wait_for_lazy_cleanup(Xai_V.get_world());

    TArray4d raap_kl = ijxy("<a a'|r|k l>");
    TArray4d rmap_kl = ijxy("<m a'|r|k l>");

    Xam_V("a,m") =  // 1/2 R^kl_ac' g^mc'_kl
                    ( R_C1 * raap_kl("a,a',k,l") + R_C2 * raap_kl("a,a',l,k"))
                    * _4("<k l|g|m a'>")
                    // 1/2 R^kl_mc' g^ac'_kl
                  + ( R_C1 * rmap_kl("m,a',k,l") + R_C2 * rmap_kl("m,a',l,k"))
                    * _4("<k l|g|a a'>");

    return XaiAddToXam(Xam_V, Xai_V);
  }

  // Xam contribution from F12 X part
  template <typename T>
  typename SingleReference_R12Intermediates<T>::TArray2
  SingleReference_R12Intermediates<T>::Xam_X(const double C_0, const double C_1) {

    const double RR_C1 = 0.5 * C_0 * C_0 + 1.5 * C_1 * C_1;
    const double RR_C2 = 0.5 * C_0 * C_0 - 1.5 * C_1 * C_1;

    TArray4d r2_akjl = ijxy("<a k|r2|j l>");
    TArray4d r2_kajl = ijxy("<k a|r2|j l>");
    TArray2 F_ij = xy("<i|F|j>");

    TArray4d r_pqjk = ijxy("<p q|r|j k>");
    TArray4d r_apnjk = ijxy("<a' n|r|j k>");
    TArray4d r_napjk = ijxy("<n a'|r|j k>");

    TArray2 X_ai, X_am;

    // - 1/2 (F^i_i + F^k_k) R^ak_alpha beta R_ik^alpha beta
    X_ai("a,i") = - (RR_C1 * r2_akjl("a,k,j,l") + RR_C2 * r2_kajl("k,a,j,l"))
                    * _2("<l|I|k>") * F_ij("j,i")

                  - (RR_C1 * r2_akjl("a,k,i,l") + RR_C2 * r2_kajl("k,a,i,l"))
                    * F_ij("l,k")

                  + (RR_C1 * _4("<a k|r|p q>") + RR_C2 * _4("<k a|r|p q>"))
                    * ( r_pqjk("p,q,j,k") * F_ij("j,i")
                      + r_pqjk("p,q,i,l") * F_ij("l,k"))

                  + (RR_C1 * _4("<a k|r|a' n>") + RR_C2 * _4("<k a|r|a' n>"))
                    * ( r_apnjk("a',n,j,k") * F_ij("j,i")
                      + r_apnjk("a',n,i,l") * F_ij("l,k"))

                  + (RR_C1 * _4("<a k|r|n a'>") + RR_C2 * _4("<k a|r|n a'>"))
                    * ( r_napjk("n,a',j,k") * F_ij("j,i")
                      + r_napjk("n,a',i,l") * F_ij("l,k"));

    // 1/2 (F^k_k + F^l_l) R^ab'_kl R_mb'^kl
    TArray4d r_aapkl = ijxy("<a a'|r|k l>");
    TArray4d r_klmap = ijxy("<k l|r|m a'>");
    X_am("a,m") = (RR_C1 * r_aapkl("a,b',k,l") + RR_C2 * r_aapkl("a,b',l,k"))
                 * ( r_klmap("j,l,m,b'") * F_ij("j,k")
                   + r_klmap("k,j,m,b'") * F_ij("j,l"));

    return XaiAddToXam(X_am, X_ai);
  }

  // Xam contribution from F12 B part
  template <typename T>
  typename SingleReference_R12Intermediates<T>::TArray2
  SingleReference_R12Intermediates<T>::Xam_B(const double C_0, const double C_1) {

    const double B_C1 = 0.5 * C_0 * C_0 + 1.5 * C_1 * C_1;
    const double B_C2 = 0.5 * C_0 * C_0 - 1.5 * C_1 * C_1;

    TArray2 B_ai, B_am;

    // 1st parts of B_am:
    // - R^ak_A'B' F^A'_C' R^ik_C'B' (A': all virtual orbital index)
    TArray4d r_akPQ = ijxy("<a k|r|p' q'>");
    TArray4d r_kaPQ = ijxy("<k a|r|p' q'>");
    TArray4d rik_PKQ = ijxy("<i k|r|p'_K(r') q'>");

    TArray4d r_akPn = ijxy("<a k|r|p' n>");
    TArray4d r_kaPn = ijxy("<k a|r|p' n>");
    TArray4d rik_PFn = ijxy("<i k|r|p'_F(r') n>");

    TArray4d r_akmA = ijxy("<a k|r|m a'>");
    TArray4d r_kamA = ijxy("<k a|r|m a'>");
    TArray4d rik_mFA = ijxy("<i k|r|m_F(n) a'>");

    TArray4d r_akpq = ijxy("<a k|r|p b>");
    TArray4d r_kapq = ijxy("<k a|r|p b>");
    TArray4d rik_pFb = ijxy("<i k|r|p_F(r) b>");

    TArray4d rik_nFA = ijxy("<i k|r|n_F(p') a'>");
    //
    TArray4d r_ikmA = ijxy("<i k|r|m a'>");
    TArray4d rak_nFA = ijxy("<a k|r|n_F(p') a'>");
    TArray4d rka_nFA = ijxy("<k a|r|n_F(p') a'>");
    TArray4d r_akAb = ijxy("<a k|r|a' b>");
    TArray4d r_kaAb = ijxy("<k a|r|a' b>");
    TArray4d rik_AFb = ijxy("<i k|r|a'_F(q) b>");
    //
    TArray4d r_ikAb = ijxy("<i k|r|a' b>");
    TArray4d rak_AFb = ijxy("<a k|r|a'_F(q) b>");
    TArray4d rka_AFb = ijxy("<k a|r|a'_F(q) b>");

    B_ai("a,i") =  //          diag
                 - ( B_C1 * _4("<a k|rTr|i l>") + B_C2 * _4("<k a|rTr|i l>")
                   //           Q
                   + ( B_C1 * (_4("<a_hJ(p') l|r2|i k>") + _4("<a l_hJ(p')|r2|i k>")
                             + _4("<a l|r2|i_hJ(p') k>") + _4("<a l|r2|i k_hJ(p')>"))

                     + B_C2 * (_4("<l a_hJ(p')|r2|i k>") + _4("<l_hJ(p') a|r2|i k>")
                             + _4("<l a|r2|i_hJ(p') k>") + _4("<l a|r2|i k_hJ(p')>"))
                     ) * 0.5
                   ) * _2("<k|I|l>")
                   //           rKr_p'q'
                 + ( B_C1 * (r_akPQ("a,k,p',q'") * rik_PKQ("i,k,p',q'")
                           + r_kaPQ("k,a,p',q'") * rik_PKQ("k,i,p',q'"))

                   + B_C2 * (r_kaPQ("k,a,p',q'") * rik_PKQ("i,k,p',q'")
                           + r_akPQ("a,k,p',q'") * rik_PKQ("k,i,p',q'"))
                   )
                   //           rFr_p'n
                 + ( B_C1 * (r_akPn("a,k,p',n") * rik_PFn("i,k,p',n")
                           + r_kaPn("k,a,p',n") * rik_PFn("k,i,p',n"))

                   + B_C2 * (r_kaPn("k,a,p',n") * rik_PFn("i,k,p',n")
                           + r_akPn("a,k,p',n") * rik_PFn("k,i,p',n"))
                   )
                   //           rFr_mA
                 - ( B_C1 * (r_akmA("a,k,m,a'") * rik_mFA("i,k,m,a'")
                           + r_kamA("k,a,m,a'") * rik_mFA("k,i,m,a'"))

                   + B_C2 * (r_kamA("k,a,m,a'") * rik_mFA("i,k,m,a'")
                           + r_akmA("a,k,m,a'") * rik_mFA("k,i,m,a'"))
                   )
                   //           rFr_pb
                 + ( B_C1 * (r_akpq("a,k,p,b") * rik_pFb("i,k,p,b")
                           + r_kapq("k,a,p,b") * rik_pFb("k,i,p,b"))

                   + B_C2 * (r_kapq("k,a,p,b") * rik_pFb("i,k,p,b")
                           + r_akpq("a,k,p,b") * rik_pFb("k,i,p,b"))
                   )
                   //
                 + ( B_C1 * (r_akmA("a,k,n,a'") * rik_nFA("i,k,n,a'")
                           + r_kamA("k,a,n,a'") * rik_nFA("k,i,n,a'"))

                   + B_C2 * (r_kamA("k,a,n,a'") * rik_nFA("i,k,n,a'")
                           + r_akmA("a,k,n,a'") * rik_nFA("k,i,n,a'"))

                   + B_C1 * (r_ikmA("i,k,n,a'") * rak_nFA("a,k,n,a'")
                           + r_ikmA("k,i,n,a'") * rka_nFA("k,a,n,a'"))

                   + B_C2 * (r_ikmA("k,i,n,a'") * rak_nFA("a,k,n,a'")
                           + r_ikmA("i,k,n,a'") * rka_nFA("k,a,n,a'"))
                    )
                   //
                 + ( B_C1 * (r_akAb("a,k,a',b") * rik_AFb("i,k,a',b")
                           + r_kaAb("k,a,a',b") * rik_AFb("k,i,a',b"))

                   + B_C2 * (r_kaAb("k,a,a',b") * rik_AFb("i,k,a',b")
                           + r_akAb("a,k,a',b") * rik_AFb("k,i,a',b"))

                  + B_C1 * (r_ikAb("i,k,a',b") * rak_AFb("a,k,a',b")
                          + r_ikAb("k,i,a',b") * rka_AFb("k,a,a',b"))

                  + B_C2 * (r_ikAb("k,i,a',b") * rak_AFb("a,k,a',b")
                          + r_ikAb("i,k,a',b") * rka_AFb("k,a,a',b"))
                  );

//    // test codes for computing B_ai for H2O molecule
//   const char* a = "a";
//   const char* i = "i";
//   TArray2 B_ia = Bpk_qk(i,a);
//   TArray2 B_ai = Bpk_qk(a,i);
//   TArray2 B_ai2 = B_ai("a,i") + B_ia("i,a");
//
//   TArray4 B_akil = Bpr_qs(a,i);
//   TArray2 B_ai2 = B_akil("a,k,i,l") * _2("<k|I|l");
//
//   const char* i = "i";
//   const char* j = "j";
//   TArray4 B_ijkl = Bpr_qs(i,j);
//
//   double sum_Bijij = 0;
//   std::cout << "B_ijkl" << std::endl;
//   for (std::size_t i = 0; i < 5; ++i) {
//     for (std::size_t j = 0; j < 5; ++j) {
//       std::vector<std::size_t> indices(4);
//       indices[0] = indices[2] = i;
//       indices[1] = indices[3] = j;
//       sum_Bijij += get_element(B_ijkl, indices);
//     }
//   }
//   std::cout << "Bijkl sum: " << sum_Bijij << std::endl;

    // 2nd parts of B_am:
    //    1/2 R^ab'_kl F^c'_m R^c'b'_kl
    // + (1/2 R^mb'_kl F^C'_a R^C'b'_kl + 1/2 R^mb'_kl F^c'_b' R^ab'_kl)
    TArray4d r_abpkl = ijxy("<a b'|r|k l>");
    TArray4d raFpbp_kl = ijxy("<a_F(c') b'|r|k l>");
    TArray4d raFbp_kl = ijxy("<a_F(c) b'|r|k l>");
    TArray4d rabpFp_kl = ijxy("<a b'_F(c')|r|k l>") ;

    B_am("a,m") =  (B_C1 * r_abpkl("a,b',k,l") + B_C2 * r_abpkl("a,b',l,k"))
                   * _4("<k l|r|m_F(c') b'>")
                   //
                   //+ (B_C1 * _4("<a_F(C') b'|r|k l>") + B_C2 *_4("<b' a_F(C')|r|k l>"))
                   //   * _4("<k l|r|i b'>") // do not work for CCR12
                 + (  B_C1 * raFpbp_kl("a,b',k,l") + B_C2 * raFpbp_kl("a,b',l,k")
                    + B_C1 * raFbp_kl("a,b',k,l") + B_C2 * raFbp_kl("a,b',l,k"))
                   * _4("<k l|r|m b'>")
                   //
                 + (B_C1 * rabpFp_kl("a,b',k,l") + B_C2 * rabpFp_kl("a,b',l,k"))
                   * _4("<k l|r|m b'>")
                 ;

    return XaiAddToXam(B_am, B_ai);
  }

  // frozen-core F12 Xii' contribution
  template <typename T>
  typename SingleReference_R12Intermediates<T>::TArray2
  SingleReference_R12Intermediates<T>::Xiip_VBX(const double C_0, const double C_1) {

    const double RR_C1 = 0.5 * C_0 * C_0 + 1.5 * C_1 * C_1;
    const double RR_C2 = 0.5 * C_0 * C_0 - 1.5 * C_1 * C_1;

    TArray2 Xiip_V;
    const char* ip = "i'";
    const char* i = "i";
    TArray2 V_ipi = VRk_Sk(ip,i,C_0, C_1);
    TArray2 V_iip = VRk_Sk(i,ip,C_0, C_1);
    Xiip_V("i,i'") = // 2 * (1/2 R^i'k_A'C' g^A'C'_ik + 1/2 R^ik_A'C' g^A'C'_i'k)
                     2.0 * (V_ipi("i',i") + V_iip("i,i'"));
      // test
//      TArray4d gr_ipk_ij = ijxy("<i' k|gr|i j>");
//      TArray4d rpq_kl = ijxy("<p q|r|k l>");
//      TArray4d rapn_kl = ijxy("<a' n|r|k l>");
//      TArray4d rnap_kl = ijxy("<n a'|r|k l>");
//      TArray4d rpq_ipk = ijxy("<p q|r|i' k>");
//      TArray2 Ikl = xy("<k|I|l>");
//      Xiip_V("i,i'") =  // 1/2 R^ik_AC g^AC_i'k
//                      ( (R_C1 * gr_ipk_ij("i',k,i,l") + R_C2 * gr_ipk_ij("i',k,l,i"))
//                        * Ikl("k,l")
//                      - _4("<i' k|g|p q>")
//                         * ( R_C1 * rpq_kl("p,q,i,k") + R_C2 * rpq_kl("p,q,k,i") )
//                      - _4("<i' k|g|a' n>")
//                         * ( R_C1 * rapn_kl("a',n,i,k") + R_C2 * rapn_kl("a',n,k,i"))
//                      - _4("<i' k|g|n a'>")
//                         * ( R_C1 * rnap_kl("n,a',i,k") + R_C2 * rnap_kl("n,a',k,i"))
//                      ) * 2.0
//                      // 1/2 R^i'k_AC g^AC_ik
//                    + ( (R_C1 * gr_ipk_ij("i',l,i,k") + R_C2 * gr_ipk_ij("i',l,k,i"))
//                        * Ikl("k,l")
//                      - _4("<i k|g|p q>")
//                        * ( R_C1 * rpq_ipk("p,q,i',k") + R_C2 * rpq_ipk("q,p,i',k"))
//                      - _4("<i k|g|a' n>")
//                        * ( R_C1 * _4("<a' n|r|i' k>") + R_C2 * _4("<a' n|r|k i'>"))
//                      - _4("<i k|g|n a'>")
//                        * ( R_C1 * _4("<n a'|r|i' k>") + R_C2 * _4("<n a'|r|k i'>"))
//                      ) * 2.0;

    TArray2 Xiip_X;
    TArray2 F_ij = xy("<i|F|j>");
    TArray4d r2_ipk_jl = ijxy("<i' k|r2|j l>");
    TArray4d r_ipk_pq = ijxy("<i' k|r|p q>");
    TArray4d r_pqkl = ijxy("<p q|r|k l>");
    TArray4d r_ipk_apn = ijxy("<i' k|r|a' n>");
    TArray4d r_apn_kl = ijxy("<a' n|r|k l>");
    TArray4d r_kip_apn = ijxy("<k i'|r|a' n>");
    Xiip_X("i,i'") = // 2 * (F^i_i + F^k_k) 1/2 R^i'k_A'B' R_ik^A'B'
                     2.0 * (  (RR_C1 * r2_ipk_jl("i',k,j,l") + RR_C2 * r2_ipk_jl("i',k,l,j"))
                              * _2("<l|I|k>") * F_ij("j,i")
                            + (RR_C1 * r2_ipk_jl("i',k,i,l") + RR_C2 * r2_ipk_jl("i',k,l,i"))
                              * F_ij("l,k")

                            - (RR_C1 * r_ipk_pq("i',k,p,q") + RR_C2 * r_ipk_pq("i',k,q,p"))
                              * ( r_pqkl("p,q,j,k") * F_ij("j,i")
                                + r_pqkl("p,q,i,l") * F_ij("l,k"))

                            - (RR_C1 * r_ipk_apn("i',k,a',n") + RR_C2 * r_kip_apn("k,i',a',n"))
                              * ( r_apn_kl("a',n,j,k") * F_ij("j,i")
                                + r_apn_kl("a',n,i,l") * F_ij("l,k"))

                            - (RR_C1 * r_kip_apn("k,i',a',n") + RR_C2 * r_ipk_apn("i',k,a',n"))
                              * ( r_apn_kl("a',n,k,j") * F_ij("j,i")
                                + r_apn_kl("a',n,l,i") * F_ij("l,k"))
                           );

    TArray2 Xiip_B;
    TArray2 B_ipi = BPk_Qk(ip,i,C_0,C_1);
    Xiip_B("i,i'") = // 2 R^i'k_A'B' F^A'_C' R^ik_C'B'
                     2.0 * B_ipi("i',i");

    TArray2 Xiip_f12;
    Xiip_f12("i,i'") = Xiip_V("i,i'") - Xiip_X("i,i'") + Xiip_B("i,i'");
    return Xiip_f12;
  }

  // frozen-core Xii' contribution from CCSD F12 coupling
  template <typename T>
  typename SingleReference_R12Intermediates<T>::TArray2
  SingleReference_R12Intermediates<T>::Xiip_CVT(const double C_0, const double C_1,
                                                   const TArray2& T1, const TArray4& T2){

    const double R_C1 = (0.5 * C_0 + 1.5 * C_1);
    const double R_C2 = (0.5 * C_0 - 1.5 * C_1);

    TArray2 Xiip_CT2;
    Xiip_CT2("i,i'") = // 2 F^a'_b R^i'l_a'c T2^bc_il
                       2.0 * ( _4("<i' l|r|b_F(a') c>")
                               * (R_C1 * T2("b,c,i,l") + R_C2 * T2("b,c,l,i"))
                              +  _4("<i' l|r|c b_F(a')>")
                               * (R_C1 * T2("c,b,i,l") + R_C2 * T2("c,b,l,i"))
                             );

    TArray2 Xiip_VT2;
    const char* ip = "i'";
    const char* i = "i";
    const char* l = "l";
    const char* c = "c";
    const char* d = "d";
    TArray4 V_ipl_cd = VPq_Rs(ip,l,c,d,C_0, C_1);
    Xiip_VT2("i,i'") = // 2 * 1/4 R^i'l_A'B' g^A'B'_cd T^cd_il
                      2.0 * V_ipl_cd("i',l,c,d") * T2("c,d,i,l");

    TArray2 Xiip_VT1;
    TArray4 V_ipl_ic = VPq_Rs(ip,l,i,c, C_0, C_1);
    TArray4 V_il_ipc = VPq_Rs(i,l,ip,c,C_0, C_1);
    TArray2 V_ipc = VRk_Sk(ip,c,C_0, C_1);
    Xiip_VT1("i,i'") =         // 2 *
                       2.0 * ( // 1/2 R^i'l_A'B' g^A'B'_ic t1^c_l
                               // 1/2 R^il_A'B' g^A'B'_i'c t1^c_l
                               (V_ipl_ic("i',l,i,c") + V_il_ipc("i,l,i',c")) * T1("c,l")

                                // 1/2 R^i'k_A'B' g^A'B'_ck t1^c_i
                              + V_ipc("i',c") * T1("c,i")
                             );

    TArray2 Xiip_C;
    Xiip_C("i,i'") = Xiip_CT2("i,i'") + Xiip_VT2("i,i'") + Xiip_VT1("i,i'");
    return Xiip_C;
  }

  // Xam contribution of CT2 part resulted from CCSD F12 coupling
  template <typename T>
  typename SingleReference_R12Intermediates<T>::TArray2
  SingleReference_R12Intermediates<T>::Xam_CT2_ccsd(const double C_0, const double C_1,
                                                    const TArray4& T2, const TArray2& RT2_aPb) {
    const double R_C1 = (0.5 * C_0 + 1.5 * C_1);
    const double R_C2 = (0.5 * C_0 - 1.5 * C_1);

    TArray4d g_abmaP = ijxy("<a b|g|m a'>");
    TArray2 Xam_CT2, Xai_CT2;
    Xam_CT2("a,m") = 2.0 * ( // + 1/2 F^a'_m R^kl_a'b T2^ab_kl
                             // + 1/2 F^a'_b R^kl_ma' T2^ab_kl
                             (R_C1 * T2("a,b,k,l") + R_C2 * T2("b,a,k,l"))
                             * (_4("<k l|r|m_F(a') b>") + _4("<k l|r|m b_F(a')>"))
                             //RT2_aPb("a',a") * _2("<a'|F|m>")

                              // - g^aa'_mb 1/2 R^kl_a'c T2^bc_kl
                              // - g^ab_ma' 1/2 R^kl_a'c T2^bc_kl
                            - (  2.0 * _4("<a a'|g|m b>") - _4("<a a'|g|b m>")
                               + 2.0 * g_abmaP("a,b,m,a'") - g_abmaP("b,a,m,a'")
                              ) * RT2_aPb("a',b")
                           );
    Xai_CT2("a,i") = 2.0 * (  // - F^a'_b R^al_a'c T2^bc_il
                            - _4("<a l|r|b_F(a') c>")
                              * (R_C1 * T2("b,c,i,l") + R_C2 * T2("b,c,l,i"))
                            - _4("<a l|r|c b_F(a')>")
                              * (R_C1 * T2("c,b,i,l") + R_C2 * T2("c,b,l,i"))
                           );

    return XaiAddToXam(Xam_CT2, Xai_CT2);
  }

  // Xam contribution of VT(T1&T2) part resulted from CCSD F12 coupling
  template <typename T>
  typename SingleReference_R12Intermediates<T>::TArray2
  SingleReference_R12Intermediates<T>::Xam_VT_ccsd(const double C_0, const double C_1,
                                                   const TArray2& T1, const TArray4& T2) {

    const double R_C1 = (0.5 * C_0 + 1.5 * C_1);
    const double R_C2 = (0.5 * C_0 - 1.5 * C_1);

    const char* i = "i";
    const char* j = "j";
    const char* k = "k";
    const char* l = "l";
    const char* m = "m";

    const char* a = "a";
    const char* c = "c";
    const char* d = "d";

    TArray4 V_alic = VPq_Rs(a,l,i,c, C_0, C_1);
    TArray4 V_ilac = VPq_Rs(i,l,a,c,C_0, C_1);
    TArray2 V_ac = VRk_Sk(a,c,C_0, C_1);
    TArray2 V_jm = VRk_Sk(j,m,C_0, C_1);
    TArray4d r_klmaP = ijxy("<k l|r|m a'>");
    TArray4d r_abPkl = ijxy("<a b'|r|k l>");
    TArray2 Xam_VT1, Xai_VT1;
    Xam_VT1("a,m") = 2.0 * (  // + 1/2 R^jk_A'B' g^A'B'_mk t1^a_j
                              T1("a,j") * V_jm("j,m")
//                             // test
//                             ( (R_C1 * _4("<j k|gr|m l>") + R_C2 * _4("<k j|gr|m l>"))
//                               * _2("<k|I|l>")
//                             - (R_C1 * _4("<j k|r|p q>") + R_C2 * _4("<k j|r|p q>"))
//                               * _4("<p q|g|m k>")
//                             - (R_C1 * _4("<j k|r|a' n>") + R_C2 * _4("<k j|r|a' n>"))
//                               * _4("<a' n|g|m k>")
//                             - (R_C1 * _4("<j k|r|n a'>") + R_C2 * _4("<k j|r|n a'>"))
//                               * _4("<n a'|g|m k>")
//                             )

                            + (  // + R^kl_mb' g^ab'_kc t1^c_l
                                 _4("<a b'|g|k c>")
                                 * (R_C1 * r_klmaP("k,l,m,b'") + R_C2 * r_klmaP("l,k,m,b'"))
                               + _4("<a b'|g|c k>")
                                 * (R_C1 * r_klmaP("l,k,m,b'") + R_C2 * r_klmaP("k,l,m,b'"))
                               // + R^kl_ab' g^mb'_kc t1^c_l
                               + (R_C1 * r_abPkl("a,b',k,l") + R_C2 * r_abPkl("a,b',l,k"))
                                 * _4("<k c|g|m b'>")
                               + (R_C1 * r_abPkl("a,b',l,k") + R_C2 * r_abPkl("a,b',k,l"))
                                 * _4("<c k|g|m b'>")
                              ) * T1("c,l")
                           );
    Xai_VT1("a,i") = 2.0 * (  // - 1/2 R^al_A'B' g^A'B'_ic t1^c_l
                              // - 1/2 R^il_A'B' g^A'B'_ac t1^c_l
                            - (V_alic("a,l,i,c") + V_ilac("i,l,a,c")) * T1("c,l")

                              // - 1/2 R^ak_A'B' g^A'B'_ck t1^c_i
                            - V_ac("a,c") * T1("c,i")
//                             // test
//                            ( (R_C1 * _4("<a k|gr|c l>") + R_C2 * _4("<k a|gr|c l>"))
//                              * _2("<k|I|l>")
//                            - (R_C1 * _4("<a k|r|p q>") + R_C2 * _4("<k a|r|p q>"))
//                              * _4("<p q|g|c k>")
//                            - (R_C1 * _4("<a k|r|a' n>") + R_C2 * _4("<k a|r|a' n>"))
//                              * _4("<a' n|g|c k>")
//                            - (R_C1 * _4("<a k|r|n a'>") + R_C2 * _4("<k a|r|n a'>"))
//                              * _4("<n a'|g|c k>")
//                            )
                           );
    TArray2 Xam_VT1_tot = XaiAddToXam(Xam_VT1, Xai_VT1);

    TArray4 V_alcd = VPq_Rs(a,l,c,d,C_0, C_1);
    TArray4 V_klmd = VPq_Rs(k,l,m,d,C_0, C_1);
    TArray2 Xam_VT2, Xai_VT2;
    Xam_VT2("a,m") = 2.0 * (  //   1/4 R^kl_A'B' g^A'B'_md T^ad_kl
                              T2("a,d,k,l") * V_klmd("k,l,m,d")

                            + (  // + 1/4 R^kl_ab' g^mb'_cd T^cd_kl
                                 (R_C1 * r_abPkl("a,b',k,l") + R_C2 * r_abPkl("a,b',l,k"))
                                 * _4("<c d|g|m b'>")
                                 // + 1/4 R^kl_mb' g^ab'_cd T^cd_kl
                               + _4("<a b'|g|c d>")
                                 * (R_C1 * r_klmaP("k,l,m,b'") + R_C2 * r_klmaP("l,k,m,b'"))
                              ) * T2("c,d,k,l")
                           );
    Xai_VT2("a,i") = 2.0 * (  // - 1/4 R^al_A'B' g^A'B'_cd T^cd_il
                            - V_alcd("a,l,c,d") * T2("c,d,i,l")
                           );
    TArray2 Xam_VT2_tot = XaiAddToXam(Xam_VT2, Xai_VT2);

    TArray2 Xam_VT_tot;
    Xam_VT_tot("a,m") = Xam_VT1_tot("a,m") + Xam_VT2_tot("a,m");
    return Xam_VT_tot;
  }

  // compute T1 & T2 amplitudes of CC2
  template <typename T>
  void SingleReference_R12Intermediates<T>::compute_T_cc2(TArray2& t1, TArray4& t2) {

    TArray2 fij = xy("<j|F|j>");
    TArray2 fab = xy("<a|F|b>");
    TArray2 fai = xy("<a|F|i>");

    // compute Delta_ai = 1 / (- <a|F|a> + <i|F|i>)
    typedef detail::diag_precond2<double> pceval_type;
    pceval_type Delta_ai_gen(TA::array_to_eigen(fab), TA::array_to_eigen(fij));

    typedef TA::Array<T, 2, LazyTensor<T, 2, pceval_type > > TArray2dLazy;
    TArray2dLazy Delta_ai(fai.get_world(), fai.trange());

    typedef TA::Array<T, 2, LazyTensor<T, 2, pceval_type > > TArray2d;
    // construct local tiles
    for(auto t = Delta_ai.trange().tiles().begin();
        t != Delta_ai.trange().tiles().end(); ++t)
      if (Delta_ai.is_local(*t)) {
        std::array<std::size_t, 2> index;
        std::copy(t->begin(), t->end(), index.begin());
        madness::Future < typename TArray2dLazy::value_type >
          tile((LazyTensor<T, 2, pceval_type >(&Delta_ai, index, &Delta_ai_gen)
              ));

        // Insert the tile into the array
        Delta_ai.set(*t, tile);
      }

    // compute Delta_ijab =  1 / (- <a|F|a> - <b|F|b> + <i|F|i> + <j|F|j>)
    typedef detail::diag_precond4<double> pc4eval_type;
    pc4eval_type Delta_abij_gen(TA::array_to_eigen(fab), TA::array_to_eigen(fab),
                                TA::array_to_eigen(fij),TA::array_to_eigen(fij));

    TArray4d g_abij = ijxy("<a b|g|i j>");
    typedef TA::Array<T, 4, LazyTensor<T, 4, pc4eval_type > > TArray4dLazy;
    TArray4dLazy Delta_abij(g_abij.get_world(), g_abij.trange());

    // construct local tiles
    for(auto t = Delta_abij.trange().tiles().begin();
        t != Delta_abij.trange().tiles().end(); ++t)
      if (Delta_abij.is_local(*t)) {
        std::array<std::size_t, 4> index;
        std::copy(t->begin(), t->end(), index.begin());
        madness::Future < typename TArray4dLazy::value_type >
          tile((LazyTensor<T, 4, pc4eval_type >(&Delta_abij, index, &Delta_abij_gen)
              ));

        // Insert the tile into the array
        Delta_abij.set(*t, tile);
      }

    // compute initial T1 & T2 amplitudes and energy
    t1("a,i") = fai("a,i") * Delta_ai("a,i");
    t2("a,b,i,j") = g_abij("a,b,i,j") * Delta_abij("a,b,i,j");
    TArray4 tau;
    tau("a,b,i,j") = t2("a,b,i,j") + t1("a,i") * t1("b,j");

    double E_0 = 0.0;
    double E_1 = 2.0 * dot(fai("a,i"), t1("a,i"))
                + dot((2.0 * g_abij("a,b,i,j") - g_abij("b,a,i,j")), tau("a,b,i,j") );
    double Delta_E = std::abs(E_0 - E_1);

    TArray2 Iij = xy("<i|I|j>");
    TArray2 Iab = xy("<a|I|b>");

    TArray4d g_iabj = ijxy("<i a|g|b j>");
    TArray4d g_iajb = ijxy("<i a|g|j b>");
    TArray4d g_iacd = ijxy("<i a|g|c d>");
    TArray4d g_aicd = ijxy("<a i|g|c d>");
    TArray4d g_klai = ijxy("<k l|g|a i>");
    TArray4d g_klia = ijxy("<k l|g|i a>");
    TArray4d g_ijkl = ijxy("<i j|g|k l>");
    TArray4d g_abcd = ijxy("<a b|g|c d>");

    TArray2 hac, hki, hck;
    TArray4 a_klij, b_abcd;

    double iter = 0;
    std::cout << indent << scprintf("%-5s", "Iter") << scprintf("%-20s", "Delta_E")
              << scprintf("%-20s", "E(CC2)") << std::endl;
    std::cout << indent << scprintf("%-5.0f", iter)   << scprintf("%-20.10f", Delta_E)
              << scprintf("%-15.10f", E_1) << std::endl;

    while (Delta_E >= 1.0e-12) {

      // recompute T1 and T2 amplitudes
      hac("a,c") = //- fab("a,c") * (Iab("a,c") - 1.0)
                  - (2.0 * g_abij("c,d,k,l") - g_abij("c,d,l,k")) * tau("a,d,k,l");
      hki("k,i") = //fij("k,i") * (1.0 - Iab("k,i"))
                  //+
                   (2.0 * g_abij("c,d,k,l") - g_abij("d,c,k,l")) * tau("c,d,i,l");
      hck("c,k") = fai("c,k")
                  + (2.0 * g_abij("c,d,k,l") - g_abij("d,c,k,l")) * t1("d,l");

      t1("a,i") = Delta_ai("a,i") * (
                   //
                   fai("a,i") - 2.0 * fai("c,k") * t1("a,k") * t1("c,i")
                   //
                 + t1("c,i") * hac("a,c") - t1("a,k") * hki("k,i")
                   //
                 + hck("c,k")
                   * (2.0 * t2("c,a,k,i") - t2("c,a,i,k") + t1("c,i") * t1("a,k"))
                   //
                 + (2.0 * g_iabj("k,a,c,i") - g_iajb("k,a,i,c")) * t1("c,k")
                   //
                 + (2.0 * g_iacd("k,a,c,d") - g_iacd("k,a,d,c")) * tau("c,d,k,i")
                   //
                 - (2.0 * g_klai("k,l,c,i") - g_klai("l,k,c,i")) * tau("c,a,k,l")
                  );

      a_klij("k,l,i,j") =  g_ijkl("k,l,i,j")
                         + g_klia("k,l,i,c") * t1("c,j") + g_klai("k,l,c,j") * t1("c,i")
                         + g_abij("c,d,k,l") * t1("c,i") * t1("d,j");
      b_abcd("a,b,c,d") =  g_abcd("a,b,c,d")
                         - g_aicd("a,k,c,d") * t1("b,k") - g_iacd("k,b,c,d") * t1("a,k");

      t2("a,b,i,j") = Delta_abij("a,b,i,j") * (
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
                     + (g_iacd("i,c,a,b") - g_iajb("k,b,i,c") * t1("a,k")) * t1("c,j")
                     + (g_iacd("j,c,b,a") - g_iajb("k,a,j,c") * t1("b,k")) * t1("c,i")
                       //
                     - (g_klai("i,j,a,k") + g_iabj("i,c,a,k") * t1("c,j")) * t1("b,k")
                     - (g_klai("j,i,b,k") + g_iabj("j,c,b,k") * t1("c,i")) * t1("a,k")
                      );
      tau("a,b,i,j") = t2("a,b,i,j") + t1("a,i") * t1("b,j");

      // recompute energy
      E_0 = E_1;
      E_1 = 2.0 * dot(fai("a,i"), t1("a,i"))
           + dot((2.0 * g_abij("a,b,i,j") - g_abij("b,a,i,j")), tau("a,b,i,j") );
      Delta_E = std::abs(E_0 - E_1);
      iter += 1;
      std::cout << indent << scprintf("%-5.0f", iter) << scprintf("%-20.10f", Delta_E)
                << scprintf("%-15.10f", E_1) << std::endl;
    }
    //std::cout << "CC2 T1 amplitudes: " << std::endl << t1 << std::endl;

  }

  // compute lambda amplitudes of CC2
  // using formula from JCP 87, 5361 (1987)
  template <typename T>
  void SingleReference_R12Intermediates<T>::compute_lambda_cc2(const TArray2& t1, const TArray4& t2,
                                                               TArray2& L1, TArray4& L2) {

    TArray2 fij = xy("<j|F|j>");
    TArray2 fab = xy("<a|F|b>");
    TArray2 fai = xy("<a|F|i>");

    // compute Delta_ai = 1 / (- <a|F|a> + <i|F|i>)
    typedef detail::diag_precond2<double> pceval_type;
    pceval_type Delta_ai_gen(TA::array_to_eigen(fab), TA::array_to_eigen(fij));

    typedef TA::Array<T, 2, LazyTensor<T, 2, pceval_type > > TArray2dLazy;
    TArray2dLazy Delta_ai(fai.get_world(), fai.trange());

    typedef TA::Array<T, 2, LazyTensor<T, 2, pceval_type > > TArray2d;
    // construct local tiles
    for(auto t = Delta_ai.trange().tiles().begin();
        t != Delta_ai.trange().tiles().end(); ++t)
      if (Delta_ai.is_local(*t)) {
        std::array<std::size_t, 2> index;
        std::copy(t->begin(), t->end(), index.begin());
        madness::Future < typename TArray2dLazy::value_type >
          tile((LazyTensor<T, 2, pceval_type >(&Delta_ai, index, &Delta_ai_gen)
              ));

        // Insert the tile into the array
        Delta_ai.set(*t, tile);
      }

    // compute Delta_ijab =  1 / (- <a|F|a> - <b|F|b> + <i|F|i> + <j|F|j>)
    typedef detail::diag_precond4<double> pc4eval_type;
    pc4eval_type Delta_abij_gen(TA::array_to_eigen(fab), TA::array_to_eigen(fab),
                                TA::array_to_eigen(fij),TA::array_to_eigen(fij));

    TArray4d g_abij = ijxy("<a b|g|i j>");
    typedef TA::Array<T, 4, LazyTensor<T, 4, pc4eval_type > > TArray4dLazy;
    TArray4dLazy Delta_abij(g_abij.get_world(), g_abij.trange());

    // construct local tiles
    for(auto t = Delta_abij.trange().tiles().begin();
        t != Delta_abij.trange().tiles().end(); ++t)
      if (Delta_abij.is_local(*t)) {
        std::array<std::size_t, 4> index;
        std::copy(t->begin(), t->end(), index.begin());
        madness::Future < typename TArray4dLazy::value_type >
          tile((LazyTensor<T, 4, pc4eval_type >(&Delta_abij, index, &Delta_abij_gen)
              ));

        // Insert the tile into the array
        Delta_abij.set(*t, tile);
      }

    L1("a,i") = t1("a,i");
    L2("a,b,i,j") = t2("a,b,i,j");

    TArray4 tau;
    tau("a,b,i,j") = t2("a,b,i,j") + t1("a,i") * t1("b,j");

    // pseudo energy
    double E_0_L = 0.0;
    double E_1_L = dot((2.0 * g_abij("a,b,i,j") - g_abij("b,a,i,j")), L2("a,b,i,j") );
    double Delta_E_L = std::abs(E_0_L - E_1_L);


    TArray4d g_abci = ijxy("<a b|g|c i>"); // g_aicd = ijxy("<a i|g|c d>");
    TArray4d g_aikl = ijxy("<a i|g|k l>"); // g_klia = ijxy("<k l|g|i a>");
    TArray4d g_aibj = ijxy("<a i|g|b j>"); // g_iajb = ijxy("<i a|g|j b>");
    TArray4d g_aijb = ijxy("<a i|g|j b>"); // g_iabj = ijxy("<i a|g|b j>");
    TArray4d g_ijkl = ijxy("<i j|g|k l>");
    TArray4d g_abcd = ijxy("<a b|g|c d>");

    double iter = 0;
    std::cout << indent << scprintf("%-5s", "Iter") << scprintf("%-20s", "Delta_E_L")
              << scprintf("%-20s", "E_L(CC2)") << std::endl;
    std::cout << indent << scprintf("%-5.0f", iter) << scprintf("%-20.10f", Delta_E_L)
              << scprintf("%-15.10f", E_1_L) << std::endl;

    // intermediates for L1 & L2
    TArray2 Iai, Gac, Kki, Iai_L2;
    TArray4 R_akic, Q_akic, A_klim, E_klic, B_aecd, F_akcd;
    TArray4 H_abcj, J_kbij;

    Iai("a,i") =  fai("a,i")
                + (2.0 * g_abij("a,c,i,k") - g_abij("a,c,k,i")) * t1("c,k");
    //TArray4 I_abij; I_abij("a,b,i,j") = g_abij("a,b,i,j");

    Gac("a,c") = //- fab("a,c") * (1.0 - Iab("a,c")) +
                 // fai("a,k") * t1("c,k")
                - (2.0 * g_abci("a,d,c,k") - g_abci("d,a,c,k")) * t1("d,k")
                + (2.0 * g_abij("a,d,k,l") - g_abij("a,d,l,k")) * tau("c,d,k,l");

    Kki("k,i") =  //fij("k,i") * (1.0 - Iij("k,i")) +
                  fai("c,i") * t1("c,k")
                + (2.0 * g_aikl("c,k,l,i") - g_aikl("c,k,i,l")) * t1("c,l")
                + (2.0 * g_abij("c,d,i,l") - g_abij("c,d,l,i")) * tau("c,d,k,l");

    R_akic("a,k,i,c") =  g_aibj("a,k,c,i") + g_abci("a,d,c,i") * t1("d,k")
                       - g_aikl("a,k,l,i") * t1("c,l")
                       - g_abij("a,d,l,i") * tau("c,d,l,k")
                       ;

    Q_akic("a,k,i,c") =  g_aijb("a,k,i,c") + g_abci("d,a,c,i") * t1("d,k")
                       - g_aikl("a,k,i,l") * t1("c,l")
                       - g_abij("a,d,l,i") * t2("c,d,k,l")
                       + g_abij("a,d,i,l") * (2.0 * t2("c,d,k,l") - tau("c,d,l,k"))
                       ;

    A_klim("k,l,i,m") = - (2.0 * g_ijkl("k,l,i,m") - g_ijkl("k,l,m,i"))
                        - (2.0 * g_abij("c,d,i,m") - g_abij("c,d,m,i")) * t1("c,k") * t1("d,l")
                        - (2.0 * g_aikl("c,l,i,m") - g_aikl("c,l,m,i")) * t1("c,k")
                        - (2.0 * g_aikl("c,k,m,i") - g_aikl("c,k,i,m")) * t1("c,l");

    E_klic("k,l,i,c") =  2.0 * g_aikl("c,i,l,k") - g_aikl("c,i,k,l")
                         // a^kl_im tcm
                       + A_klim("k,l,i,m")  * t1("c,m")
                        //
                       + (2.0 * g_abci("e,d,c,i") - g_abci("d,e,c,i")) * t1("d,k") * t1("e,l")
                        //
                       + (2.0 * g_aijb("d,l,i,c") - g_aibj("d,l,c,i")) * t1("d,k")
                       + (2.0 * g_aibj("d,k,c,i") - g_aijb("d,k,i,c")) * t1("d,l")
                       ;

    B_aecd("a,e,c,d") = - (2.0 * g_abcd("a,e,c,d") - g_abcd("e,a,c,d"))
                        - (2.0 * g_abij("a,e,k,l") - g_abij("e,a,k,l")) * t1("c,k") * t1("d,l")
                        + (2.0 * g_abci("e,a,d,k") - g_abci("a,e,d,k")) * t1("c,k")
                        + (2.0 * g_abci("a,e,c,k") - g_abci("e,a,c,k")) * t1("d,k");

    F_akcd("a,k,c,d") = - (2.0 * g_abci("c,d,a,k") - g_abci("d,c,a,k"))
                          // b^ae_cd tek
                        + B_aecd("a,e,c,d") * t1("e,k")
                          //
                        - (2.0 * g_aikl("a,k,m,l") - g_aikl("a,k,l,m")) * t1("d,l") * t1("c,m")
                          //
                        + (2.0 * g_aibj("a,k,c,l") - g_aijb("a,k,l,c")) * t1("d,l")
                        + (2.0 * g_aijb("a,k,l,d") - g_aibj("a,k,d,l")) * t1("c,l")
                        ;

    H_abcj("a,b,c,j") = - g_abci("a,b,c,j")
                        + g_abij("a,b,k,j") * t1("c,k");

    J_kbij("k,b,i,j") =  g_aikl("b,k,j,i")
                       + g_abij("b,c,j,i")  * t1("c,k");

    while (Delta_E_L >= 1.0e-12) {

      // close-shell L^a_i & alpha-beta case L^ab_ij
      L1("a,i") = Delta_ai("a,i") * (
                  // Iai
                  Iai("a,i")
                  // Gac Lci
                - Gac("a,c") * L1("c,i")
                  // Kki Lak
                - Kki("k,i") * L1("a,k")
                  // (R - 2Q)ak_ic Lck
                - (R_akic("a,k,i,c") - 2.0 * Q_akic("a,k,i,c")) * L1("c,k")
                  // E^kl_ic L2^ac_kl
                - E_klic("k,l,i,c") * L2("a,c,k,l")
                  // F^ak_cd L2^cd_ik
                - F_akcd("a,k,c,d") * L2("c,d,i,k")
                  );

      L2("a,b,i,j") = Delta_abij("a,b,i,j") * (
                      // I^ab_ij
                      g_abij("a,b,i,j")

                      // Ibi Laj + Iaj Lbi
                    //- L1("a,j") * Iai("b,i")
                    //- L1("b,i") * Iai("a,j")

                      // Ibj Lai + Iai Lbj
                    + L1("a,i") * Iai("b,j")
                    + L1("b,j") * Iai("a,i")

                      // Gbc L^ac_ij + Gac L^bc_ji
                      // CC2: Gbc = - fbc (1 - delta_bc) = 0
                    //+ fab("b,c") * (1.0 - Iab("b,c")) * L2("a,c,i,j")
                    //+ fab("a,c") * (1.0 - Iab("a,c")) * L2("b,c,j,i")

                      // Kkj L^ab_ik + Kki L^ba_jk
                      // CC2: Kkj = fkj (1 - delta_kj) = 0
                    //- fij("k,j") * (1.0 - Iij("k,j")) * L2("a,b,i,k")
                    //- fij("k,i") * (1.0 - Iij("k,i")) * L2("b,a,j,k")

                      // H^ab_cj Lci + H^ba_ci Lcj
                    - H_abcj("a,b,c,j") * L1("c,i")
                    - H_abcj("b,a,c,i") * L1("c,j")

                      // J^kb_ij Lak - j^ka_ji Lbk
                    - J_kbij("k,b,i,j") * L1("a,k")
                    - J_kbij("k,a,j,i") * L1("b,k")
                    );

#if 0  // close-shell L^a_i & alpha-beta case L^ab_ij in explicit form

      L1("a,i") = Delta_ai("a,i") * (
                  // Iai
                  fai("a,i")
                + (2.0 * g_abij("a,b,i,j") - g_abij("a,b,j,i")) * t1("b,j")

                  // Gac Lci
                - (  fai("a,j") * t1("c,j")
                   - (2.0 * g_abci("a,b,c,j") - g_abci("b,a,c,j")) * t1("b,j")
                   + (2.0 * g_abij("a,d,k,l") - g_abij("a,d,l,k")) * tau("c,d,k,l"))
                  * L1("c,i")

                  // Kki Lak
                - (  fai("c,i") * t1("c,k")
                   + (2.0 * g_aikl("b,k,j,i") - g_aikl("b,k,i,j")) * t1("b,j")
                   + (2.0 * g_abij("c,d,i,l") - g_abij("c,d,l,i")) * tau("c,d,k,l"))
                  * L1("a,k")

                  // (R - 2Q)ak_ic Lck
                + (  2.0 * g_aijb("a,k,i,c") - g_aibj("a,k,c,i")
                  + (2.0 * g_abci("b,a,c,i") - g_abci("a,b,c,i")) * t1("b,k")
                  - (2.0 * g_aikl("a,k,i,l") - g_aikl("a,k,l,i")) * t1("c,l")
                  - (2.0 * g_abij("a,b,i,j") - g_abij("a,b,j,i")) * t1("c,j") * t1("b,k")
                  + (2.0 * g_abij("a,b,i,j") - g_abij("a,b,j,i"))
                    * (2.0 *  t2("c,b,k,j") - t2("c,b,j,k"))
                  ) * L1("c,k")

                  // E^kl_ic L2^ac_kl
                - ( 2.0 * g_aikl("b,i,l,k") - g_aikl("b,i,k,l")
                    // a^kl_im tcm
                  + ( - (2.0 * g_ijkl("k,l,i,j") - g_ijkl("k,l,j,i"))
                      - (2.0 * g_abij("c,d,i,j") - g_abij("c,d,j,i")) * t1("c,k") * t1("d,l")

                      - (2.0 * g_aikl("c,l,i,j") - g_aikl("c,l,j,i")) * t1("c,k")
                      - (2.0 * g_aikl("c,k,j,i") - g_aikl("c,k,i,j")) * t1("c,l")
                    ) * t1("b,j")
                   //
                  + (2.0 * g_abci("d,c,b,i") - g_abci("c,d,b,i")) * t1("d,l") * t1("c,k")
                   //
                  + (2.0 * g_aijb("c,l,i,b") - g_aibj("c,l,b,i")) * t1("c,k")
                  + (2.0 * g_aibj("c,k,b,i") - g_aijb("c,k,i,b")) * t1("c,l")
                  ) * L2("a,b,k,l")

                  // F^ak_cd L2^cd_ik
                - ( - 2.0 * g_abci("c,d,a,j") - g_abci("d,c,a,j")
                    // b^ae_cd tek
                  + ( - (2.0 * g_abcd("a,b,c,d") - g_abcd("b,a,c,d"))
                      - (2.0 * g_abij("a,b,k,l") - g_abij("b,a,k,l")) * t1("c,k") * t1("d,l")

                      + (2.0 * g_abci("b,a,d,k") - g_abci("a,b,d,k")) * t1("c,k")
                      + (2.0 * g_abci("a,b,c,k") - g_abci("b,a,c,k")) * t1("d,k")
                    ) * t1("b,j")
                    //
                  - (2.0 * g_aikl("a,j,k,l") - g_aikl("a,j,l,k")) * t1("d,l") * t1("c,k")
                    //
                  + (2.0 * g_aibj("a,j,c,k") - g_aijb("a,j,k,c")) * t1("d,k")
                  + (2.0 * g_aijb("a,j,k,d") - g_aibj("a,j,d,k")) * t1("c,k")
                  ) * L2("c,d,i,j")
                  );

      L2("a,b,i,j") = Delta_abij("a,b,i,j") * (
                      // I^ab_ij
                      g_abij("a,b,i,j")

                      // Ibi Laj + Iaj Lbi
                      // Ibj Lai + Iai Lbj

                    + (  fai("a,i")
                       + (2.0 * g_abij("a,f,i,n") - g_abij("a,f,n,i")) * t1("f,n"))
                      * L1("b,j")

                    + (  fai("b,j")
                       + (2.0 * g_abij("b,f,j,n") - g_abij("b,f,n,j")) * t1("f,n"))
                      * L1("a,i")

                      // H^ab_cj Lci + H^ba_ci Lcj
                    - (- g_abci("a,b,e,j")
                       + g_abij("a,b,n,j") * t1("e,n")) * L1("e,i")
                    - (- g_abci("b,a,e,i")
                       + g_abij("b,a,n,i") * t1("e,n")) * L1("e,j")

                      // J^kb_ij Lak - j^ka_ji Lbk
                    - ( g_aikl("b,m,j,i")
                      + g_abij("b,f,j,i")  * t1("f,m")) * L1("a,m")
                    - ( g_aikl("a,m,i,j")
                      + g_abij("a,f,i,j")  * t1("f,m")) * L1("b,m")
                    );
#endif

      E_0_L = E_1_L;
      E_1_L = dot((2.0 * g_abij("a,b,i,j") - g_abij("b,a,i,j")), L2("a,b,i,j") );
      Delta_E_L = std::abs(E_0_L - E_1_L);

      iter += 1;
      std::cout << indent << scprintf("%-5.0f", iter) << scprintf("%-20.10f", Delta_E_L)
                << scprintf("%-15.10f", E_1_L) << std::endl;
    }
    //std::cout << "CC2 L1 amplitudes: " << std::endl << L1 << std::endl;
  }

  // compute lambda amplitudes of CC2
  // use formula from Gauss and Stanton, JCP, 103 (1995)
  template <typename T>
  void SingleReference_R12Intermediates<T>::compute_lambda_cc2_2(const TArray2& t1, const TArray4& t2,
                                                               TArray2& L1, TArray4& L2) {

    TArray2 fij = xy("<j|F|j>");
    TArray2 fab = xy("<a|F|b>");
    TArray2 fai = xy("<a|F|i>");

    // compute Delta_ai = 1 / (- <a|F|a> + <i|F|i>)
    typedef detail::diag_precond2<double> pceval_type;
    pceval_type Delta_ai_gen(TA::array_to_eigen(fab), TA::array_to_eigen(fij));

    typedef TA::Array<T, 2, LazyTensor<T, 2, pceval_type > > TArray2dLazy;
    TArray2dLazy Delta_ai(fai.get_world(), fai.trange());

    typedef TA::Array<T, 2, LazyTensor<T, 2, pceval_type > > TArray2d;
    // construct local tiles
    for(auto t = Delta_ai.trange().tiles().begin();
        t != Delta_ai.trange().tiles().end(); ++t)
      if (Delta_ai.is_local(*t)) {
        std::array<std::size_t, 2> index;
        std::copy(t->begin(), t->end(), index.begin());
        madness::Future < typename TArray2dLazy::value_type >
          tile((LazyTensor<T, 2, pceval_type >(&Delta_ai, index, &Delta_ai_gen)
              ));

        // Insert the tile into the array
        Delta_ai.set(*t, tile);
      }

    // compute Delta_ijab =  1 / (- <a|F|a> - <b|F|b> + <i|F|i> + <j|F|j>)
    typedef detail::diag_precond4<double> pc4eval_type;
    pc4eval_type Delta_abij_gen(TA::array_to_eigen(fab), TA::array_to_eigen(fab),
                                TA::array_to_eigen(fij),TA::array_to_eigen(fij));

    TArray4d g_abij = ijxy("<a b|g|i j>");
    typedef TA::Array<T, 4, LazyTensor<T, 4, pc4eval_type > > TArray4dLazy;
    TArray4dLazy Delta_abij(g_abij.get_world(), g_abij.trange());

    // construct local tiles
    for(auto t = Delta_abij.trange().tiles().begin();
        t != Delta_abij.trange().tiles().end(); ++t)
      if (Delta_abij.is_local(*t)) {
        std::array<std::size_t, 4> index;
        std::copy(t->begin(), t->end(), index.begin());
        madness::Future < typename TArray4dLazy::value_type >
          tile((LazyTensor<T, 4, pc4eval_type >(&Delta_abij, index, &Delta_abij_gen)
              ));

        // Insert the tile into the array
        Delta_abij.set(*t, tile);
      }

    // pseudo energy
    L1("a,i") = t1("a,i");
    L2("a,b,i,j") = t2("a,b,i,j");

    double E_0_L = 0.0;
    double E_1_L = dot((2.0 * g_abij("a,b,i,j") - g_abij("b,a,i,j")), L2("a,b,i,j") );
    double Delta_E_L = std::abs(E_0_L - E_1_L);

    double iter = 0;
    std::cout << indent << scprintf("%-5s", "Iter") << scprintf("%-20s", "Delta_E_L")
              << scprintf("%-20s", "E_L(CC2)") << std::endl;
    std::cout << indent << scprintf("%-5.0f", iter) << scprintf("%-20.10f", Delta_E_L)
              << scprintf("%-15.10f", E_1_L) << std::endl;

    // intermediates for computing L1 & L2
    TArray4d g_abci = ijxy("<a b|g|c i>");
    TArray4d g_aikl = ijxy("<a i|g|k l>");
    TArray4d g_aibj = ijxy("<a i|g|b j>");
    TArray4d g_aijb = ijxy("<a i|g|j b>");
    TArray4d g_ijkl = ijxy("<i j|g|k l>");
    TArray4d g_abcd = ijxy("<a b|g|c d>");

    TArray4 tau2_aa, tau2_ab, Ttau_aa, Ttau_ab;
    tau2_aa("a,b,i,j") =  t2("a,b,i,j") - t2("b,a,i,j")
                        + (t1("a,i") * t1("b,j") - t1("b,i") * t1("a,j"));
    tau2_ab("a,b,i,j") =  t2("a,b,i,j")
                        + (t1("a,i") * t1("b,j"));

    Ttau_aa("a,b,i,j") =  t2("a,b,i,j") - t2("b,a,i,j")
                        + 0.5 * (t1("a,i") * t1("b,j")  - t1("b,i") * t1("a,j"));
    Ttau_ab("a,b,i,j") =  t2("a,b,i,j")
                        + 0.5 * t1("a,i") * t1("b,j");

    TArray4 tau2t1_aa, tau2t1_ab, Ttaut1_aa, Ttaut1_ab;
    tau2t1_aa("a,b,i,j") =  t1("a,i") * t1("b,j") - t1("b,i") * t1("a,j");
    tau2t1_ab("a,b,i,j") =  t1("a,i") * t1("b,j");
    Ttaut1_aa("a,b,i,j") =  0.5 * t1("a,i") * t1("b,j")  - t1("b,i") * t1("a,j");
    Ttaut1_ab("a,b,i,j") =  0.5 * t1("a,i") * t1("b,j");

    TArray2 TFme, TFae, TFmi;
    TFme("m,e") =  //   fme
                   // + t^f_n g^mn_ef
                  (2.0 * g_abij("e,f,m,n") - g_abij("e,f,n,m")) * t1("f,n");
    TFae("a,e") =  //   fae - 1/2 fme t^a_m
                   // + t^f_m g^am_ef
                   t1("f,m") * (2.0 * g_abci("e,f,a,m") - g_abci("f,e,a,m"))
                   // - 1/2 Ttau^af_mn g^mn_ef
                 - 0.5 * Ttau_aa("a,f,m,n") * (g_abij("e,f,m,n") - g_abij("f,e,m,n"))
                 - Ttau_ab("a,f,m,n") * g_abij("e,f,m,n")
                 ;
    TFmi("m,i") =  //   fmi + 1/2 fme t^e_i
                   // + t^e_n g^mn_ie
                   t1("e,n") * (2.0 * g_aikl("e,i,n,m") - g_aikl("e,i,m,n"))
                   // + 1/2 Ttau^ef_in g^mn_ef
                 + 0.5 * Ttau_aa("e,f,i,n") * (g_abij("e,f,m,n") - g_abij("e,f,n,m"))
                 + Ttau_ab("e,f,i,n") * g_abij("e,f,m,n")
                 ;

    TArray2 CFme, CFae, CFmi;
    CFme("m,e") = TFme("m,e");
    CFae("a,e") = TFae("a,e") - 0.5 * t1("a,m") * TFme("m,e");
    CFmi("m,i") = TFmi("m,i") + 0.5 * t1("e,i") * TFme("m,e");

    // compute \cal{W}mbej
    TArray4 TW_MBEJ_aa, TW_MbEj_ab,CW_MBEJ_aa, CW_MbEj_ab;

    // \tilde{W}mbej:
    // alpha-alpha case
    TW_MBEJ_aa("m,b,e,j") =  //   g^mb_ej
                             g_aijb("e,j,m,b") - g_aibj("e,j,b,m")

                             // + t^f_j g^mb_ef
                           + t1("f,j") * (g_abci("f,e,b,m") - g_abci("e,f,b,m"))

                             // - t^b_n g^mn_ej
                           - t1("b,n") * (g_aikl("e,j,m,n") - g_aikl("e,j,n,m"))

                             // - 1/2 (t^fb_jn + 2 t^f_j t^b_n) g^mn_ef:
                             // + 1/2 t^bf_jn g^mn_ef
                           + 0.5 * (  (2.0 * t2("b,f,j,n") - t2("f,b,j,n"))
                                      * g_abij("e,f,m,n")
                                    - (t2("b,f,j,n") - t2("f,b,j,n"))
                                      * g_abij("f,e,m,n")
                                   )
                             // - t^f_j t^b_n g^mn_ef
                           - t1("f,j") * t1("b,n")
                             * (g_abij("e,f,m,n") - g_abij("f,e,m,n"))
                           ;
    // alpha-beta case
    TW_MbEj_ab("m,b,e,j") = // g^mb_ej
                            g_aijb("e,j,m,b")
                            // + t^f_j g^mb_ef
                          + t1("f,j") * g_abci("f,e,b,m")
                            // - t^b_n g^mn_ej
                          - t1("b,n") * g_aikl("e,j,m,n")

                            // - 1/2 (t^fb_jn + 2 t^f_j t^b_n) g^mn_ef:
                            // + 1/2 t^bf_jn g^mn_ef
                          + 0.5 * (  (2.0 * t2("b,f,j,n") - t2("f,b,j,n"))
                                     * g_abij("e,f,m,n")
                                   - t2("b,f,j,n") * g_abij("f,e,m,n")
                                   )
                            // - t^f_j t^b_n g^mn_ef
                          -  t1("f,j") * t1("b,n") * g_abij("e,f,m,n")
                          ;

    // \cal{W}mbej
    CW_MBEJ_aa("m,b,e,j") =  //   \tilde{W}mbej
                             TW_MBEJ_aa("m,b,e,j")

                             // - 1/2 t^fb_jn g^mn_ef = + 1/2 t^fb_nj g^mn_ef
                           + 0.5 * (t2("f,b,n,j") - t2("b,f,n,j"))
                                 * (g_abij("e,f,m,n") - g_abij("f,e,m,n"))
                           + 0.5 * t2("f,b,n,j") * g_abij("e,f,m,n")
                           ;
    CW_MbEj_ab("m,b,e,j") =  //   \tilde{W}mbej
                             TW_MbEj_ab("m,b,e,j")

                             // - 1/2 t^fb_jn g^mn_ef = + 1/2 t^fb_nj g^mn_ef
                           + 0.5 * (t2("f,b,n,j") - t2("b,f,n,j"))
                                 * g_abij("e,f,m,n")
                           + 0.5 * t2("f,b,n,j")
                                 * (g_abij("e,f,m,n") - g_abij("f,e,m,n"))
                           ;


    // compute \cal{W}abei
    TArray4 TW_ABEF_aa, TW_AbEf_ab, CW_ABEF_aa, CW_AbEf_ab, CW_ABEI_aa, CW_AbEi_ab;

    // \tilde{W}abef
    TW_ABEF_aa("a,b,e,f") =  //  g^ab_ef
                             g_abcd("a,b,e,f") - g_abcd("b,a,e,f")

                             // - P(ab) t^b_m g^am_ef:
                             // - t^b_m g^am_ef
                           - t1("b,m") * (g_abci("e,f,a,m") - g_abci("f,e,a,m"))
                             // + t^a_m g^bm_ef
                           + t1("a,m") * (g_abci("e,f,b,m") - g_abci("f,e,b,m"))

                             // + 1/4 tau^ab_mn g^mn_ef
                           + 0.25 * tau2t1_aa("a,b,m,n") * (g_abij("e,f,m,n") - g_abij("f,e,m,n"))
                           ;

    TW_AbEf_ab("a,b,e,f") =  //  g^ab_ef
                             g_abcd("a,b,e,f")

                             // - P(ab) t^b_m g^am_ef:
                             // - t^b_m g^am_ef
                           - t1("b,m") * g_abci("e,f,a,m")
                             // - t^a_m g^mb_ef
                           - t1("a,m") * g_abci("f,e,b,m")

                             // + 1/4 tau^ab_mn g^mn_ef
                           + 0.5 * tau2t1_ab("a,b,m,n") * g_abij("e,f,m,n")
                           ;

    // \cal{W}abef
    CW_ABEF_aa("a,b,e,f") =  // \tilde{W}abef
                             TW_ABEF_aa("a,b,e,f")
                             // + 1/4 tau^ab_mn g^mn_ef
                           + 0.25 * tau2t1_aa("a,b,m,n") * (g_abij("e,f,m,n") - g_abij("f,e,m,n"));

    CW_AbEf_ab("a,b,e,f") =  // \tilde{W}abef
                             TW_AbEf_ab("a,b,e,f")
                             // + 1/4 tau^ab_mn g^mn_ef
                           + 0.5 * tau2t1_ab("a,b,m,n") * g_abij("e,f,m,n");

    // \cal{W}abei
    CW_ABEI_aa("a,b,e,i") =  //   g^ab_ei
                             g_abci("a,b,e,i") - g_abci("b,a,e,i")

                             // - \cal{F}me t^ab_mi
                           //- CFme("m,e") * (t2("a,b,m,i") - t2("b,a,m,i"))

                             // + t^f_i \cal{W}abef
                           + t1("f,i") * CW_ABEF_aa("a,b,e,f")

                             // + 1/2 g^mn_ei tau^ab_mn
                           + 0.5 * (g_aikl("e,i,m,n") - g_aikl("e,i,n,m")) * tau2t1_aa("a,b,m,n")

                             // - P(ab) g^mb_ef t^af_mi (=> 0 for cc2):
                             // - g^mb_ef t^af_mi = - g^bm_ef t^fa_mi
                           //- (g_abci("f,e,b,m") - g_abci("e,f,b,m")) * (t2("a,f,m,i") - t2("f,a,m,i")) // m,f in alpha space
                           //- g_abci("e,f,b,m") * t2("f,a,m,i") // m,f in beta space
                             // + g^ma_ef t^bf_mi = + g^am_ef t^fb_mi
                           //+ g_abci("f,e,a,m") * t2("b,f,m,i") - g_abci("e,f,a,m") * t2("b,f,m,i")
                           //- g_abci("f,e,a,m") * t2("f,b,m,i") + 2.0 * g_abci("e,f,a,m") * t2("f,b,m,i")

                             // - P(ab) t^a_m (g^mb_ei - t^bf_ni g^mn_ef):
                             // - t^a_m (g^mb_ei - t^bf_ni g^mn_ef) = - t^a_m (g^mb_ei + t^fb_ni g^mn_ef)
                           - t1("a,m")
                             * (  g_aijb("b,m,i,e") - g_aibj("b,m,e,i")
                                //- (  t2("b,f,n,i") * g_abij("e,f,m,n") - 2.0 * t2("f,b,n,i") * g_abij("e,f,m,n")
                                //   - t2("b,f,n,i") * g_abij("f,e,m,n") + t2("f,b,n,i") * g_abij("f,e,m,n"))
                               )
                             // + t^b_m (g^ma_ei - t^af_ni g^mn_ef) = + t^b_m (- g^am_ei + t^af_in g^mn_ef)
                           + t1("b,m")
                             * (- g_aibj("a,m,e,i") + g_aijb("a,m,i,e")
                                //+ (  2.0 * t2("f,a,n,i") * g_abij("e,f,m,n") - t2("a,f,n,i") * g_abij("e,f,m,n")
                                //   - t2("f,a,n,i") * g_abij("f,e,m,n") + t2("a,f,n,i") * g_abij("f,e,m,n")
                                //   )
                                )
                           ;

    CW_AbEi_ab("a,b,e,i") =  //   g^ab_ei
                             g_abci("a,b,e,i")

                             // - \cal{F}me t^ab_mi
                           //- CFme("m,e") * t2("a,b,m,i")

                             // + t^f_i \cal{W}abef
                           + t1("f,i") * CW_AbEf_ab("a,b,e,f")

                             // + 1/2 g^mn_ei tau^ab_mn
                           + g_aikl("e,i,m,n") * tau2t1_ab("a,b,m,n")

                             // - P(ab) g^mb_ef t^af_mi (=> 0 for cc2):
                             // - g^mb_ef t^af_mi
                           //- g_abci("f,e,b,m") * t2("a,f,m,i")
                             // + g^ma_ef t^bf_mi = + g^am_ef t^fb_mi
                           //+ (g_abci("e,f,a,m") - g_abci("f,e,a,m")) * t2("f,b,m,i")
                           //+ g_abci("e,f,a,m") * (t2("f,b,m,i") - t2("b,f,m,i"))

                             // - P(ab) t^a_m (g^mb_ei - t^bf_ni g^mn_ef):
                             // - t^a_m (g^mb_ei - t^bf_ni g^mn_ef) = - t^a_m (g^mb_ei + t^fb_ni g^mn_ef)
                           - t1("a,m")
                             * (  g_aijb("b,m,i,e")
                                //+ t2("f,b,n,i") * (g_abij("e,f,m,n") - g_abij("f,e,m,n"))
                                //+ (t2("f,b,n,i") - t2("b,f,n,i")) * g_abij("e,f,m,n")
                                )
                              // + t^b_m (g^ma_ei - t^af_ni g^mn_ef) = + t^b_m (- g^am_ei + t^af_ni g^nm_ef)
                           + t1("b,m")
                             * (- g_aibj("a,m,e,i")
                                //+ t2("a,f,n,i") * g_abij("e,f,n,m")
                               )
                           ;

    // compute \cal{W}mbij
    TArray4 TW_MNIJ_aa, TW_MnIj_ab, CW_MNIJ_aa, CW_MnIj_ab, CW_MBIJ_aa, CW_MbIj_ab;

    // \tilde{W}mnij
    TW_MNIJ_aa("m,n,i,j") =  // g^mn_ij
                             g_ijkl("m,n,i,j") - g_ijkl("n,m,i,j")

                             // + P(ij) t^e_j g^mn_ie:
                             // + t^e_j g^mn_ie
                           + t1("e,j") * (g_aikl("e,i,n,m") - g_aikl("e,i,m,n"))
                             // + t^e_i g^mn_ej
                           + t1("e,i") * (g_aikl("e,j,m,n") - g_aikl("e,j,n,m"))

                             // + 0.25 * tau^ef_ij * g^mn_ef
                           + 0.25 * tau2t1_aa("e,f,i,j") * (g_abij("e,f,m,n") - g_abij("e,f,n,m"))
                            ;
    TW_MnIj_ab("m,n,i,j") =  // g^mn_ij
                             g_ijkl("m,n,i,j")

                             // + P(ij) t^e_j g^mn_ie:
                             // + t^e_j g^mn_ie
                           + t1("e,j") * g_aikl("e,i,n,m")
                             // + t^e_i g^mn_ej
                           + t1("e,i") * g_aikl("e,j,m,n")

                             // + 0.25 * tau^ef_ij * g^mn_ef
                           + 0.5 * tau2t1_ab("e,f,i,j") * g_abij("e,f,m,n")
                           ;

    // \cal{W}mnij
    CW_MNIJ_aa("m,n,i,j") =  //   \tilde{W}mnij
                             TW_MNIJ_aa("m,n,i,j")
                             // + 1/4 tau^ef_ij g^mn_ef
                           + 0.25 * tau2t1_aa("e,f,i,j") * (g_abij("e,f,m,n") - g_abij("e,f,n,m"))
                           ;
    CW_MnIj_ab("m,n,i,j") =  // \tilde{W}mnij
                             TW_MnIj_ab("m,n,i,j")
                             // + 1/4 tau^ef_ij g^mn_ef
                           + 0.5 * tau2t1_ab("e,f,i,j") * g_abij("e,f,m,n")
                           ;

    // \cal{W}mbij
    CW_MBIJ_aa("m,b,i,j") =  //   g^mb_ij
                             g_aikl("b,m,j,i") - g_aikl("b,m,i,j")

                             // - \cal{F}me t^be_ij = + CFme t^eb_ij
                           //+ CFme("m,e") * (t2("e,b,i,j") - t2("b,e,i,j"))

                             // - t^b_n \cal{W}mnij
                           - t1("b,n") * CW_MNIJ_aa("m,n,i,j")

                             // + 0.5 g^mb_ef tau^ef_ij
                           + 0.5 * (g_abci("f,e,b,m") - g_abci("e,f,b,m")) * tau2t1_aa("e,f,i,j")

                             // + P(ij) g^mn_ie t^be_jn (=> 0 for cc2)：
                             // + g^mn_ie t^be_jn
                           //+ 2.0 * g_aikl("e,i,n,m") * t2("b,e,j,n") - g_aikl("e,i,m,n") * t2("b,e,j,n")
                           //- g_aikl("e,i,n,m") * t2("b,e,n,j") + g_aikl("e,i,m,n") * t2("b,e,n,j")
                             // - g^mn_je t^be_in
                           //+ 2.0 * g_aikl("e,j,n,m") * t2("b,e,i,n") - g_aikl("e,j,m,n") * t2("b,e,i,n")
                           //- g_aikl("e,j,n,m") * t2("e,b,i,n") + g_aikl("e,j,m,n") * t2("e,b,i,n")

                             // + P(IJ) t^e_i (g^mb_ej - t^bf_nj g^mn_ef)
                             // + t^e_i (g^mb_ej - t^bf_nj g^mn_ef) = + t^e_i (g^mb_ej + t^bf_jn g^mn_ef)
                           + t1("e,i")
                             * (  g_aijb("e,j,m,b") - g_aibj("e,j,b,m")

                                //+ 2.0 * t2("b,f,j,n") * g_abij("e,f,m,n") - t2("f,b,j,n") * g_abij("e,f,m,n")
                                //- t2("b,f,j,n") * g_abij("f,e,m,n") + t2("f,b,j,n") * g_abij("f,e,m,n")
                               )
                             // - t^e_j (g^mb_ei - t^bf_ni g^mn_ef) = - t^e_j (g^mb_ei + t^bf_in g^mn_ef)
                           - t1("e,j")
                             * (  g_aijb("e,i,m,b") - g_aibj("e,i,b,m")

                                //+ 2.0 * t2("b,f,i,n") * g_abij("e,f,m,n") - t2("f,b,i,n") * g_abij("e,f,m,n")
                                //- t2("b,f,i,n") * g_abij("f,e,m,n") + t2("f,b,i,n") * g_abij("f,e,m,n")
                                )
                           ;

    CW_MbIj_ab("m,b,i,j") =  //   g^mb_ij
                             g_aikl("b,m,j,i")

                             // - \cal{F}me t^be_ij = + CFme t^eb_ij
                           //+ CFme("m,e") * t2("e,b,i,j")

                             // - t^b_n \cal{W}mnij
                           - t1("b,n") * CW_MnIj_ab("m,n,i,j")

                             // + 0.5 g^mb_ef tau^ef_ij
                           + g_abci("f,e,b,m") * tau2t1_ab("e,f,i,j")

                             // + P(ij) g^mn_ie t^be_jn (=> 0 for cc2):
                             // + g^mn_ie t^be_jn
                           //+ 2.0 * g_aikl("e,i,n,m") * t2("b,e,j,n") - g_aikl("e,i,n,m") * t2("b,e,n,j")
                           //- g_aikl("e,i,m,n") * t2("b,e,j,n")
                             // - g^mn_je t^be_in = - g^nm_je t^eb_in
                           //+ g_aikl("e,j,m,n") * t2("e,b,i,n")

                             // + P(IJ) t^e_i (g^mb_ej - t^bf_nj g^mn_ef)
                             // + t^e_i (g^mb_ej - t^bf_nj g^mn_ef) = + t^e_i (g^mb_ej + t^bf_jn g^mn_ef)
                           + t1("e,i")
                             * (  g_aijb("e,j,m,b")

                                //+ 2.0 * t2("b,f,j,n") * g_abij("e,f,m,n") - t2("f,b,j,n") * g_abij("e,f,m,n")
                                //- t2("b,f,j,n") * g_abij("f,e,m,n")
                               )
                              // - t^e_j (g^mb_ei - t^bf_ni g^mn_ef) = - t^e_j (- g^mb_ie + t^bf_ni g^mn_fe)
                           - t1("e,j")
                             * (- g_aibj("e,i,b,m")
                                //+ t2("b,f,n,i") * g_abij("f,e,m,n")
                               )
                           ;

    // intermediates for L2
    TArray4 CW_MnIe_ab, CW_MniE_ab, CW_AmEf_ab, CW_aMEf_ab;

    CW_MnIe_ab("m,n,i,e") =  //   g^mn_ie
                             g_aikl("e,i,n,m")
                             // + t^f_i g^mn_fe
                           + t1("f,i") * g_abij("f,e,m,n")
                           ;
    CW_MniE_ab("m,n,i,e") =  //   g^mn_ie
                           - g_aikl("e,i,m,n")
                             // + t^f_i g^mn_fe
                           - t1("f,i") * g_abij("f,e,n,m")
                           ;

    CW_AmEf_ab("a,m,e,f") =  //   g^am_ef
                             g_abci("e,f,a,m")
                             // - t^a_n g^nm_ef
                           - t1("a,n") * g_abij("e,f,n,m")
                           ;
    CW_aMEf_ab("a,m,e,f") =  //   g^am_ef
                            - g_abci("f,e,a,m")
                             // - t^a_n g^nm_ef
                           + t1("a,n") * g_abij("f,e,n,m")
                           ;

    while (Delta_E_L >= 1.0e-9) {

      L1("a,i") =  Delta_ai("a,i") * (
                   //   \cal{F}ia
                   CFme("i,a")
                   // + \lambda^i_e \cal{F}ea
                 + L1("e,i") * CFae("e,a")
                   // - \lambda^m_a \cal{F}im
                 - L1("a,m") * CFmi("i,m")
                   // + \lambda^m_e \cal{W}ieam
                 + L1("e,m") * (CW_MBEJ_aa("i,e,a,m") + CW_MbEj_ab("i,e,a,m"))
                   // + 1/2 \lambda^im_ef \cal{W}efam
                 + 0.5 * (L2("e,f,i,m") - L2("e,f,m,i") ) * CW_ABEI_aa("e,f,a,m")
                 + L2("e,f,i,m") * CW_AbEi_ab("e,f,a,m")
                   // - 1/2 \lambda^mn_ae \cal{W}iemn
                 - 0.5 * (L2("a,e,m,n") - L2("e,a,m,n") ) * CW_MBIJ_aa("i,e,m,n")
                 - L2("a,e,m,n") * CW_MbIj_ab("i,e,m,n")
                  );

      L2("a,b,i,j") =  Delta_abij("a,b,i,j") * (
                       //   g^ij_ab
                       g_abij("a,b,i,j")

                       // + P(ab) \lambda^ij_ae \cal{F}eb
                     //+ L2("a,e,i,j") * CFae("e,b")
                     //+ L2("e,b,i,j") * CFae("e,a")

                       // - P(ij) \lambda^im_ab \cal{F}jm
                     //- L2("a,b,i,m") * CFmi("j,m")
                     //- L2("a,b,m,j") * CFmi("i,m")

                     // + P(ij) \lambda^i_e \cal{W}ejab
                     + L1("e,i") * CW_AmEf_ab("e,j,a,b")
                     - L1("e,j") * CW_aMEf_ab("e,i,a,b")

                     // - P(ab) \lambda^m_a \cal{W}ijmb
                     - L1("a,m") * CW_MnIe_ab("i,j,m,b")
                     + L1("b,m") * CW_MniE_ab("i,j,m,a")

                     // + P(ij) P(ab) \lambda^i_a \cal{F}jb
                     + L1("a,i") * CFme("j,b")
                     //- L1("b,i") * CFme("j,a")
                     //- L1("a,j") * CFme("i,b")
                     + L1("b,j") * CFme("i,a")
                     );

      E_0_L = E_1_L;
      E_1_L = dot((2.0 * g_abij("a,b,i,j") - g_abij("b,a,i,j")), L2("a,b,i,j") );
      Delta_E_L = std::abs(E_0_L - E_1_L);

      iter += 1;
      std::cout << indent << scprintf("%-5.0f", iter) << scprintf("%-20.10f", Delta_E_L)
                << scprintf("%-15.10f", E_1_L) << std::endl;
    }
    //std::cout << "CC2 L1 amplitudes: " << std::endl << L1 << std::endl;
  }

  // compute CC2 one-electron density from amplitudes
  template <typename T>
  void SingleReference_R12Intermediates<T>::compute_cc2_1rdm_amp(
                           const TArray2& T1_cc2, const TArray4& T2_cc2,
                           const TArray2& L1_cc2, const TArray4& L2_cc2,
                           TArray2& Dij_cc2, TArray2& Dab_cc2,
                           TArray2& Dia_cc2, TArray2& Dai_cc2) {

    Dij_cc2("i,j") =  // - 1/4 P+(ij) t^ef_im lambda^jm_ef
                    - 0.5 * (  (2.0 *  T2_cc2("c,d,i,k") - T2_cc2("c,d,k,i"))
                               * L2_cc2("c,d,j,k")
                             + (2.0 *  T2_cc2("c,d,j,k") - T2_cc2("c,d,k,j"))
                               * L2_cc2("c,d,i,k")
                            )
                      // - 1/2 P+(ij) t^e_i lambda^j_e
                    - 0.5 * (  T1_cc2("c,i") * L1_cc2("c,j")
                             + T1_cc2("c,j") * L1_cc2("c,i")
                            )
                    ;
    //std::cout << "Dij_cc2: " << std::endl << Dij_cc2 << std::endl;

    Dab_cc2("a,b") =  //   1/4 P+(ab) t^mn_ae lambda^be_mn
                      0.5 * (  (2.0 *  T2_cc2("a,c,k,l") - T2_cc2("c,a,k,l"))
                               * L2_cc2("b,c,k,l")
                             + (2.0 *  T2_cc2("b,c,k,l") - T2_cc2("c,b,k,l"))
                               * L2_cc2("a,c,k,l")
                            )
                      // + 1/2 P+(ab) lambda^m_a t^b_m
                    + 0.5 * (  L1_cc2("a,k") * T1_cc2("b,k")
                             + L1_cc2("b,k") * T1_cc2("a,k")
                            )
                    ;

    Dia_cc2("i,a") = //   t^a_i
                     T1_cc2("a,i")
                     // + (t^ae_im - t^e_i t^a_m) \lambda^m_e
                   + (2.0 * T2_cc2("a,c,i,k") - T2_cc2("a,c,k,i")
                      - T1_cc2("c,i") * T1_cc2("a,k"))
                     * L1_cc2("c,k")
                     // CCSD term
                     // - 1/2 lambda^mn_ef (t^ef_in t^a_m + t^e_i t^af_mn)
//                   - (2.0 * L2_cc2("c,d,k,l") - L2_cc2("d,c,k,l"))
//                     * (  T2_cc2("c,d,i,l") * T1_cc2("a,k")
//                        + T2_cc2("a,d,k,l") * T1_cc2("c,i")
//                       )
                   ;

    Dai_cc2("a,i") = // lambda^i_a
                     L1_cc2("a,i");
  }

  // compute lambda amplitudes of CCSD
  // use formula from Gauss and Stanton, JCP, 103 (1995)
  // compute Delta_ai = 1 / (- <a|F|a> + <i|F|i>)
  //       & Delta_ijab =  1 / (- <a|F|a> - <b|F|b> + <i|F|i> + <j|F|j>)
  template <typename T>
    void SingleReference_R12Intermediates<T>::compute_Delta_cc(const TA::TiledRange& TR_ai,
                                                               const TA::TiledRange& TR_abij,
                                                               TArray2& D_ai, TArray4& D_abij) {
    TArray2 fij = xy("<j|F|j>");
    TArray2 fab = xy("<a|F|b>");

    // compute Delta_ai = 1 / (- <a|F|a> + <i|F|i>)
    typedef detail::diag_precond2<double> pceval_type;
    pceval_type Delta_ai_gen(TA::array_to_eigen(fab), TA::array_to_eigen(fij));

    typedef TA::Array<T, 2, LazyTensor<T, 2, pceval_type > > TArray2dLazy;
    TArray2dLazy Delta_ai(fij.get_world(), TR_ai);

    typedef TA::Array<T, 2, LazyTensor<T, 2, pceval_type > > TArray2d;
    // construct local tiles
    for(auto t = Delta_ai.trange().tiles().begin();
        t != Delta_ai.trange().tiles().end(); ++t)
      if (Delta_ai.is_local(*t)) {
        std::array<std::size_t, 2> index;
        std::copy(t->begin(), t->end(), index.begin());
        madness::Future < typename TArray2dLazy::value_type >
          tile((LazyTensor<T, 2, pceval_type >(&Delta_ai, index, &Delta_ai_gen)
              ));

        // Insert the tile into the array
        Delta_ai.set(*t, tile);
      }
    D_ai("a,i") = Delta_ai("a,i");

    // compute Delta_ijab =  1 / (- <a|F|a> - <b|F|b> + <i|F|i> + <j|F|j>)
    typedef detail::diag_precond4<double> pc4eval_type;
    pc4eval_type Delta_abij_gen(TA::array_to_eigen(fab), TA::array_to_eigen(fab),
                                TA::array_to_eigen(fij),TA::array_to_eigen(fij));

    typedef TA::Array<T, 4, LazyTensor<T, 4, pc4eval_type > > TArray4dLazy;
    TArray4dLazy Delta_abij(fij.get_world(), TR_abij);

    // construct local tiles
    for(auto t = Delta_abij.trange().tiles().begin();
        t != Delta_abij.trange().tiles().end(); ++t)
      if (Delta_abij.is_local(*t)) {
        std::array<std::size_t, 4> index;
        std::copy(t->begin(), t->end(), index.begin());
        madness::Future < typename TArray4dLazy::value_type >
          tile((LazyTensor<T, 4, pc4eval_type >(&Delta_abij, index, &Delta_abij_gen)
              ));

        // Insert the tile into the array
        Delta_abij.set(*t, tile);
      }
    D_abij("a,b,i,j") = Delta_abij("a,b,i,j");
    fij.get_world().gop.fence();
  }

  // compute intermediates needed for T amplitudes
  template <typename T>
  void SingleReference_R12Intermediates<T>::compute_intermediates_TFW_ccsd(
                                              const TArray2& t1, const TArray4& t2,
                                              const TArray4d& g_abij, const TArray4d& g_aikl,
                                              const TArray4d& g_aibj, const TArray4d& g_aijb,
                                              const TArray4d& g_abci,
                                              const TArray4d& g_ijkl, const TArray4d& g_abcd,
                                              TArray2& TFkc, TArray2& TFac, TArray2& TFki,
                                              TArray4& TW_KbCj_ab,TArray4& TW_KbcJ_ba,
                                              TArray4& TW_AbCd_ab, TArray4& TW_KlIj_ab) {
    TArray4 tau_ab, Ttau_ab;
    tau_ab("a,b,i,j") = t2("a,b,i,j") + (t1("a,i") * t1("b,j"));
    Ttau_ab("a,b,i,j") = t2("a,b,i,j") + 0.5 * t1("a,i") * t1("b,j");

    TArray4 gabij_temp;
    gabij_temp("a,b,i,j") = 2.0 * g_abij("a,b,i,j") - g_abij("a,b,j,i");

    TFkc("k,c") =  //   fkc
                   // + t^d_l g^kl_cd
                   gabij_temp("c,d,k,l") * t1("d,l");
    TFac("a,c") =  //   fac (1 - delta_ac) - 1/2 fkc t^a_k
                   // + t^d_k g^ak_cd
                   t1("d,k") * (2.0 * g_abci("c,d,a,k") - g_abci("d,c,a,k"))
                   // - 1/2 Ttau^ad_kl g^kl_cd
                 - Ttau_ab("a,d,k,l") * gabij_temp("c,d,k,l")
                 ;
    TFki("k,i") =  //   fki (1 - delta_ki) + 1/2 fkc t^c_i
                   // + t^c_l g^kl_ic
                   t1("c,l") * (2.0 * g_aikl("c,i,l,k") - g_aikl("c,i,k,l"))
                   // + 1/2 Ttau^cd_il g^kl_cd
                 + Ttau_ab("c,d,i,l") * gabij_temp("c,d,k,l")
                 ;

    // \tilde{W}mbej:
    // alpha-beta case
    TW_KbCj_ab("k,b,c,j") = // g^kb_cj
                            g_aijb("c,j,k,b")
                            // + t^d_j g^kb_cd
                          + t1("d,j") * g_abci("d,c,b,k")
                            // - t^b_l g^kl_cj
                          - t1("b,l") * g_aikl("c,j,k,l")

                            // - 1/2 (t^db_jl + 2 t^d_j t^b_l) g^kl_cd:
                            // + 1/2 t^bd_jl g^kl_cd
                          + 0.5 * (  t2("b,d,j,l") * gabij_temp("c,d,k,l")
                                   - t2("d,b,j,l") * g_abij("c,d,k,l")
                                   )
                            // - t^d_j t^b_l g^kl_cd
                          -  t1("d,j") * (t1("b,l") * g_abij("c,d,k,l"))
                          ;

    ////// CCSD only term
    TW_KbcJ_ba("k,b,c,j") = // g^kb_cj
                          - g_aibj("c,j,b,k")
                            // + t^d_j g^kb_cd
                          - t1("d,j") * g_abci("c,d,b,k")
                            // - t^b_l g^kl_cj
                          + t1("b,l") * g_aikl("c,j,l,k")

                            // - 1/2 (t^db_jl + 2 t^d_j t^b_l) g^kl_cd:
                            // + 1/2 t^bd_jl g^kl_cd
                          + 0.5 * t2("d,b,j,l") * g_abij("d,c,k,l")
                            // - t^d_j t^b_l g^kl_cd
                          + t1("d,j") * (t1("b,l") * g_abij("d,c,k,l"))
                          ;

    // \tilde{W}abef
    TW_AbCd_ab("a,b,c,d") =  //  g^ab_cd
                             g_abcd("a,b,c,d")

                             // - P(ab) t^b_k g^ak_cd:
                             // - t^b_k g^ak_cd
                           - t1("b,k") * g_abci("c,d,a,k")
                             // - t^a_k g^kb_cd
                           - t1("a,k") * g_abci("d,c,b,k")

                             ////// CCSD modified terms
                             // + 1/4 tau^ab_kl g^kl_cd
                           + 0.5 * tau_ab("a,b,k,l") * g_abij("c,d,k,l")
                           ;

    // \tilde{W}mnij
    TW_KlIj_ab("k,l,i,j") =  // g^kl_ij
                             g_ijkl("k,l,i,j")

                             // + P(ij) t^c_j g^kl_ic:
                             // + t^c_j g^kl_ic
                           + t1("c,j") * g_aikl("c,i,l,k")
                             // + t^e_i g^mn_ej
                           + t1("c,i") * g_aikl("c,j,k,l")

                             ////// CCSD modified terms
                             // + 0.25 * tau^cd_ij * g^kl_cd
                           + 0.5 * tau_ab("c,d,i,j") * g_abij("c,d,k,l")
                           ;
  }

  // compute intermediates for lambda amplitudes
  template <typename T>
  void SingleReference_R12Intermediates<T>::compute_intermediates_CW_ccsd(
                                              const TArray2& t1, const TArray4& t2,
                                              TArray2& CFkc, TArray2& CFac, TArray2& CFki,
                                              TArray4& CW_KbCj_ab, TArray4& CW_KbcJ_ba,
                                              TArray4& CW_AbCd_ab, TArray4& CW_AbCi_ab,
                                              TArray4& CW_KlIj_ab, TArray4& CW_KbIj_ab,
                                              TArray4& CW_KlIc_ab, TArray4& CW_KliC_ab,
                                              TArray4& CW_AkCd_ab) {

    TArray4 tau_ab, Ttau_ab;
    tau_ab("a,b,i,j") = t2("a,b,i,j") + (t1("a,i") * t1("b,j"));

    Ttau_ab("a,b,i,j") = t2("a,b,i,j") + 0.5 * t1("a,i") * t1("b,j");

    TArray4d g_abij = ijxy("<a b|g|i j>");
    TArray4d g_abci = ijxy("<a b|g|c i>");
    TArray4d g_aikl = ijxy("<a i|g|k l>");
    TArray4d g_aibj = ijxy("<a i|g|b j>");
    TArray4d g_aijb = ijxy("<a i|g|j b>");
    TArray4d g_ijkl = ijxy("<i j|g|k l>");
    TArray4d g_abcd = ijxy("<a b|g|c d>");

    // compute intermediates
    TArray2 TFkc, TFac, TFki;
    TArray4 TW_KbCj_ab,
            ////// CCSD only term
            TW_KbcJ_ba,
            TW_AbCd_ab, TW_KlIj_ab;

    compute_intermediates_TFW_ccsd(t1, t2,
                                   g_abij, g_aikl,
                                   g_aibj, g_aijb, g_abci,
                                   g_ijkl, g_abcd,
                                   TFkc, TFac, TFki,
                                   TW_KbCj_ab, TW_KbcJ_ba,
                                   TW_AbCd_ab, TW_KlIj_ab);

     // \cal{F}
     CFkc("k,c") = TFkc("k,c");
     CFac("a,c") = TFac("a,c") - 0.5 * t1("a,k") * TFkc("k,c");
     CFki("k,i") = TFki("k,i") + 0.5 * t1("c,i") * TFkc("k,c");

    // \cal{W}mbej
    CW_KbCj_ab("k,b,c,j") =  //   \tilde{W}kbcj
                             TW_KbCj_ab("k,b,c,j")

                             // - 1/2 t^db_jl g^kn_cd = + 1/2 t^db_lj g^kl_cd
                           + 0.5 * (t2("d,b,l,j") - t2("b,d,l,j"))
                                 * g_abij("c,d,k,l")
                           + 0.5 * t2("d,b,l,j")
                                 * (g_abij("c,d,k,l") - g_abij("d,c,k,l"))
                           ;

    CW_KbcJ_ba("k,b,c,j") =  //   \tilde{W}kbcj
                             TW_KbcJ_ba("k,b,c,j")

                             // - 1/2 t^db_jl g^kn_cd = + 1/2 t^db_lj g^kl_cd
                           + 0.5 * t2("b,d,l,j") * g_abij("d,c,k,l")
                           ;

    // \cal{W}abef
    CW_AbCd_ab("a,b,c,d") =  // \tilde{W}abcd
                             TW_AbCd_ab("a,b,c,d")

                             ////// CCSD modified terms
                             // + 1/4 tau^ab_kl g^kl_cd
                           + 0.5 * tau_ab("a,b,k,l") * g_abij("c,d,k,l");

    TArray4 t2_temp;
    t2_temp("a,b,i,j") = 2.0 * t2("a,b,i,j") - t2("a,b,j,i");
    // \cal{W}abei
    CW_AbCi_ab("a,b,c,i") =  //   g^ab_ci
                             g_abci("a,b,c,i")

                             // - \cal{F}kc t^ab_ki
                           - CFkc("k,c") * t2("a,b,k,i")

                             // + t^d_i \cal{W}abcd
                           + t1("d,i") * CW_AbCd_ab("a,b,c,d")

                             ////// CCSD modified terms
                             // + 1/2 g^kl_ci tau^ab_kl
                           + g_aikl("c,i,k,l") * tau_ab("a,b,k,l")

                             ////// CCSD only terms
                             // - P(ab) g^kb_cd t^ad_ki:
                             // - g^kb_cd t^ad_ki
                           - g_abci("d,c,b,k") * t2("a,d,k,i")
                             // + g^ka_cd t^bd_ki = + g^ak_cd t^db_ki
                           + (  g_abci("c,d,a,k") * t2_temp("b,d,i,k")
                              - g_abci("d,c,a,k") * t2("b,d,i,k"))

                             // - P(ab) t^a_k (g^kb_ci - t^bd_li g^kl_cd):
                             // - t^a_k (g^kb_ci - t^bd_li g^kl_cd) = - t^a_k (g^kb_ci + t^db_li g^kl_cd)
                           - t1("a,k")
                             * (  g_aijb("b,k,i,c")
                                  ////// CCSD only terms
                                 + (  t2_temp("d,b,l,i") * g_abij("c,d,k,l")
                                    - t2("d,b,l,i") * g_abij("d,c,k,l"))
                                )
                              // + t^b_k (g^ka_ci - t^ad_li g^kl_cd) = + t^b_k (- g^ak_ci + t^ad_li g^lk_cd)
                           + t1("b,k")
                             * (- g_aibj("a,k,c,i")
                                  ////// CCSD only terms
                                + t2("a,d,l,i") * g_abij("c,d,l,k")
                               )
                           ;

    // \cal{W}mnij
    CW_KlIj_ab("k,l,i,j") =  // \tilde{W}klij
                             TW_KlIj_ab("k,l,i,j")
                             ////// CCSD modified terms
                             // + 1/4 tau^cd_ij g^kl_cd
                           + 0.5 * tau_ab("c,d,i,j") * g_abij("c,d,k,l")
                           ;

    // \cal{W}mbij
    CW_KbIj_ab("k,b,i,j") =  //   g^kb_ij
                             g_aikl("b,k,j,i")

                             // - \cal{F}kc t^bc_ij = + CFkc t^cb_ij
                           + CFkc("k,c") * t2("c,b,i,j")

                             // - t^b_l \cal{W}klij
                           - t1("b,l") * CW_KlIj_ab("k,l,i,j")

                             ////// CCSD modified terms
                             // + 0.5 g^kb_cd tau^cd_ij
                           + g_abci("d,c,b,k") * tau_ab("c,d,i,j")

                             ////// CCSD only terms
                             // + P(ij) g^kl_ic t^bc_jl
                             // + g^kl_ic t^bc_jl
                           + g_aikl("c,i,l,k") * t2_temp("b,c,j,l")
                           - g_aikl("c,i,k,l") * t2("b,c,j,l")
                             // - g^kl_jc t^bc_il = - g^lk_jc t^cb_il
                           - g_aikl("c,j,k,l") * t2("c,b,i,l")

                             // + P(IJ) t^c_i (g^kb_cj - t^bd_lj g^kl_cd)
                             // + t^c_i (g^kb_cj - t^bd_lj g^kl_cd) = + t^c_i (g^kb_cj + t^bd_jl g^kl_cd)
                           + t1("c,i")
                             * (  g_aijb("c,j,k,b")
                                  ////// CCSD only terms
                                + t2_temp("b,d,j,l") * g_abij("c,d,k,l")
                                - t2("b,d,j,l") * g_abij("d,c,k,l")
                               )
                              // - t^c_j (g^kb_ci - t^bd_li g^kl_cd) = - t^c_j (- g^kb_ic + t^bd_li g^kl_dc)
                           - t1("c,j")
                             * (- g_aibj("c,i,b,k")
                                  ////// CCSD only terms
                                + t2("b,d,l,i") * g_abij("d,c,k,l")
                               )
                           ;

    // CCSD only term: CW_MNIE_ab
    CW_KlIc_ab("k,l,i,c") =  //   g^kl_ic
                             g_aikl("c,i,l,k")
                             // + t^d_i g^kl_dc
                           + t1("d,i") * g_abij("d,c,k,l")
                           ;
    CW_KliC_ab("k,l,i,c") =  //   g^kl_ic
                           - g_aikl("c,i,k,l")
                             // + t^d_i g^kl_dc
                           - t1("d,i") * g_abij("d,c,l,k")
                           ;

    // CCSD only term: CW_AMEF_aa
    CW_AkCd_ab("a,k,c,d") =  //   g^ak_cd
                             g_abci("c,d,a,k")
                             // - t^a_l g^lk_cd
                           - t1("a,l") * g_abij("c,d,l,k")
                           ;

  }

  // compute T1 & T2 amplitudes of CCSD
  template <typename T>
  void SingleReference_R12Intermediates<T>::compute_T_ccsd(TArray2& t1, TArray4& t2,
                                                          const std::string method) {

    // compute Delta_ai = 1 / (- <a|F|a> + <i|F|i>)
    //       & Delta_ijab =  1 / (- <a|F|a> - <b|F|b> + <i|F|i> + <j|F|j>)
    TArray2 fai = xy("<a|F|i>");
    TArray4d g_abij = ijxy("<a b|g|i j>");

    TArray2 Delta_ai;
    TArray4 Delta_abij;
    compute_Delta_cc(fai.trange(), g_abij.trange(), Delta_ai, Delta_abij);

    // compute initial T1 & T2 amplitudes and energy
    t1("a,i") = fai("a,i") * Delta_ai("a,i");
    t2("a,b,i,j") = g_abij("a,b,i,j") * Delta_abij("a,b,i,j");

    TArray2 T1_f12contri;
    TArray4 T2_f12contri,
            T2_f12contri_T1, V_KLAJ_temp, V_IJKL_temp; // F12 contributions involing T1
    if (method == "f12b" || method == "F12b") {
      const double C_0 = 1.0 / 2.0;
      const double C_1 = 1.0 / 4.0;
      const double RC1 = 0.5 * (C_0 + C_1);
      const double RC2 = 0.5 * (C_0 - C_1);

      const char* i = "i";
      const char* j = "j";
      const char* k = "k";
      const char* l = "l";
      const char* a = "a";
      const char* b = "b";

      // F12 contributions to amplitude equations
      //
      // contributions to T1 equation
      TArray2 Via = VRk_Sk(i, a, C_0, C_1);
      T1_f12contri("a,i") = Delta_ai("a,i") * Via("i,a"); // 1/2 R^ik_AB g^AB_ak
      t1("a,i") = t1("a,i") + T1_f12contri("a,i");

      // contributions to T2 equation
      TArray4d r_ijaFb = ijxy("<i j|r|a_F(c') b>");
      TArray4 V_IJAB = VPQ_RS(i, j, a, b);
      T2_f12contri("a,b,i,j") = Delta_abij("a,b,i,j") * (
                                  // F^a'_a R^ij_a'b + F^a'_b R^ij_aa'
                                  RC1 * r_ijaFb("i,j,a,b") + RC2 * r_ijaFb("j,i,a,b")
                                + RC1 * r_ijaFb("j,i,b,a") + RC2 * r_ijaFb("i,j,b,a")
                                  // + 1/2 R^ij_AB g^AB_ab
                                + RC1 * V_IJAB("i,j,a,b") + RC2 * V_IJAB("j,i,a,b")
                                );
      // compute terms needed for constructing contributions that depends on T amplitudes
      // - 1/2 R^ij_AB g^AB_ak t^k_b - 1/2 R^ij_AB g^AB_kb t^k_a
      // + 1/2 R^ij_AB g^AB_kl t^k_a t^l_b
      TArray4 V_KLAJ = VPQ_RS(k, l, a, j);
      TArray4 V_IJKL = VPQ_RS(i, j, k, l);
      V_KLAJ_temp("i,j,a,k") = RC1 * V_KLAJ("i,j,a,k") + RC2 * V_KLAJ("j,i,a,k");
      V_IJKL_temp("i,j,k,l") = RC1 * V_IJKL("i,j,k,l") + RC2 * V_IJKL("j,i,k,l");

      t2("a,b,i,j") =  t2("a,b,i,j")
                     + T2_f12contri("a,b,i,j")
                     ;
    }

    TArray4 tau_ab;
    tau_ab("a,b,i,j") = t2("a,b,i,j") + t1("a,i") * t1("b,j");

    TArray4 gabij_temp;
    gabij_temp("a,b,i,j") = 2.0 * g_abij("a,b,i,j") - g_abij("b,a,i,j");

    double E_0 = 0.0;
    double E_1 =  //2.0 * dot(fai("a,i"), t1("a,i")) +
                 dot(gabij_temp("a,b,i,j"), tau_ab("a,b,i,j") );
    double Delta_E = std::abs(E_0 - E_1);

    double iter = 0;
    std::cout << indent << scprintf("%-5s", "Iter") << scprintf("%-20s", "Delta_E")
              << scprintf("%-20s", "E(CCSD)") << std::endl;
    std::cout << indent << scprintf("%-5.0f", iter)   << scprintf("%-15.12f", Delta_E)
              << scprintf("%-15.12f", E_1) << std::endl;

    TArray2 TFkc, TFac, TFki;
    TArray4 TW_KbCj_ab,
            ////// CCSD only term
            TW_KbcJ_ba,
            TW_AbCd_ab, TW_KlIj_ab;

    TArray4d g_aikl = ijxy("<a i|g|k l>");
    TArray4d g_aibj = ijxy("<a i|g|b j>");
    TArray4d g_aijb = ijxy("<a i|g|j b>");
    TArray4d g_abci = ijxy("<a b|g|c i>");
    TArray4d g_ijkl = ijxy("<i j|g|k l>");
    TArray4d g_abcd = ijxy("<a b|g|c d>");

    TArray4 t2_unsymm, t2_interim;

    while (Delta_E >= 1.0e-12) {
      // compute intermediates
      compute_intermediates_TFW_ccsd(t1, t2,
                                     g_abij, g_aikl,
                                     g_aibj, g_aijb, g_abci,
                                     g_ijkl, g_abcd,
                                     TFkc, TFac, TFki,
                                     TW_KbCj_ab, TW_KbcJ_ba,
                                     TW_AbCd_ab, TW_KlIj_ab);

      t2_interim("a,b,i,j") = 2.0 * t2("a,b,i,j") - t2("b,a,i,j");

      auto t1_temp = Delta_ai("a,i") * (
                   //   fai
                   //fai("a,i")
                   // + TFac t^c_i
                   TFac("a,c") * t1("c,i")
                   // - TFki t^a_k
                 - TFki("k,i") * t1("a,k")
                   // + TFkc t^ac_ik
                 + TFkc("k,c") * t2_interim("a,c,i,k")
                   // + t^c_k g^ak_ic
                 + t1("c,k") * (2.0 * g_aijb("a,k,i,c") - g_aibj("a,k,c,i"))
                   // - 1/2 t^ac_kl g^kl_ic
                 - t2_interim("c,a,k,l") * g_aikl("c,i,k,l")
                   // + 1/2 t^cd_ik g^ak_cd
                 + t2_interim("c,d,i,k") * g_abci("c,d,a,k")
                  );

      auto t2_temp = Delta_abij("a,b,i,j") * (
                              //   g^ab_ij
                              g_abij("a,b,i,j")

                              // + P_(ab) t^ac_ij (TFbc - 1/2 t^b_k TFkc)
                            + t2("a,c,i,j") * (TFac("b,c") - 0.5 * t1("b,k") * TFkc("k,c"))
                            + t2("c,b,i,j") * (TFac("a,c") - 0.5 * t1("a,k") * TFkc("k,c"))

                              // - P_(ij)  T^ab_ik (TFkj + 1/2 t^c_j TFkc)
                            - t2("a,b,i,k") * (TFki("k,j") + 0.5 * t1("c,j") * TFkc("k,c"))
                            - t2("a,b,k,j") * (TFki("k,i") + 0.5 * t1("c,i") * TFkc("k,c"))

                              // + 1/2 tau^ab_kl TW_klij
                            + tau_ab("a,b,k,l") * TW_KlIj_ab("k,l,i,j")

                              // + 1/2 tau^cd_ij TW_abcd
                            + tau_ab("c,d,i,j") * TW_AbCd_ab("a,b,c,d")

                              // + P_(ij) P_(ab) (t^ac_ik TW_kbcj - t^c_i t^a_k g^kb_cj)
                            + (  t2("a,c,i,k") * TW_KbcJ_ba("k,b,c,j")
                               + t2_interim("a,c,i,k") * TW_KbCj_ab("k,b,c,j")
                               //
                               - t1("c,i") * (t1("a,k") * g_aijb("b,k,j,c"))
                              )
                              // - P(ab)
                            - (- t2("c,b,i,k") * TW_KbcJ_ba("k,a,c,j") // TW_kBCj_ba => TW_KbcJ_ba
                               + t1("c,i") * (t1("b,k") * g_aibj("a,k,c,j"))
                              )
                              // - P(ij)
                            - (- t2("a,c,k,j") * TW_KbcJ_ba("k,b,c,i")
                               + t1("c,j") * (t1("a,k") * g_aibj("b,k,c,i"))
                              )
                              // + P_(ij) P_(ab)
                            + (  t2("b,c,j,k") * TW_KbcJ_ba("k,a,c,i")
                               + t2_interim("b,c,j,k") * TW_KbCj_ab("k,a,c,i")
                                 //
                               - t1("c,j") * (t1("b,k") * g_aijb("a,k,i,c"))
                              )

                              // + P(ij) t^c_i g^ab_cj
                            + t1("c,i") * g_abci("a,b,c,j")
                            + t1("c,j") * g_abci("b,a,c,i")

                              // - P(ab) t^a_k g^kb_ij
                            - t1("a,k") * g_aikl("b,k,j,i")
                            - t1("b,k") * g_aikl("a,k,i,j")
                            );

     if (method == "f12b" || method == "F12b") {
        t1("a,i") = t1_temp
                   + T1_f12contri("a,i")
                   ;

        // recompute F12 contributions to T2 which depend on T amplitudes
        T2_f12contri_T1("a,b,i,j") = Delta_abij("a,b,i,j") * (
                                      // - 1/2 R^ij_AB g^AB_ak t^k_b - 1/2 R^ij_AB g^AB_kb t^k_a
                                    - V_KLAJ_temp("i,j,a,k") * t1("b,k")
                                    - V_KLAJ_temp("j,i,b,k") * t1("a,k")
                                      // + 1/2 R^ij_AB g^AB_kl t^k_a t^l_b
                                    + V_IJKL_temp("i,j,k,l") * t1("a,k") * t1("b,l")
                                    );

        t2_unsymm("a,b,i,j") =  t2_temp
                              + T2_f12contri("a,b,i,j")
                              + T2_f12contri_T1("a,b,i,j")
                              ;
      } else {
        t1("a,i") = t1_temp;
        t2_unsymm("a,b,i,j") = t2_temp;
      }

      t2("a,b,i,j") = 0.5 * (t2_unsymm("a,b,i,j") + t2_unsymm("b,a,j,i"));
      tau_ab("a,b,i,j") = t2("a,b,i,j") + t1("a,i") * t1("b,j");

      // recompute energy
      E_0 = E_1;
      E_1 = //2.0 * dot(fai("a,i"), t1("a,i")) +
            dot(gabij_temp("a,b,i,j"), tau_ab("a,b,i,j") );
      Delta_E = std::abs(E_0 - E_1);
      iter += 1;
      std::cout << indent << scprintf("%-5.0f", iter) << scprintf("%-15.12f", Delta_E)
                << scprintf("%-15.12f", E_1) << std::endl;
    }
    //std::cout << "CCSD T1 amplitudes: " << std::endl << t1 << std::endl;

  }

  template <typename T>
  void SingleReference_R12Intermediates<T>::compute_lambda_ccsd(const TArray2& t1, const TArray4& t2,
                                                                TArray2& L1, TArray4& L2,
                                                                const std::string method) {

    TArray4d g_abij = ijxy("<a b|g|i j>");
    TArray4d g_abci = ijxy("<a b|g|c i>");
    TArray2 fai = xy("<a|F|i>");

    // compute Delta_ai = 1 / (- <a|F|a> + <i|F|i>)
    //       & Delta_ijab =  1 / (- <a|F|a> - <b|F|b> + <i|F|i> + <j|F|j>)
    TArray2 Delta_ai;
    TArray4 Delta_abij;
    compute_Delta_cc(fai.trange(), g_abij.trange(), Delta_ai, Delta_abij);

    // initial guess for L1 and L2
    TArray2 L1_f12contri, L1_f12b, L1_f12b_L2;
    TArray4 L2_f12contri, Vijab, Vklib, Vklij;
    if (method == "ccsd" || method == "CCSD") {
      L1("a,i") = t1("a,i");
      L2("a,b,i,j") = t2("a,b,i,j");

    } else {
      const double C_0 = 1.0 / 2.0;
      const double C_1 = 1.0 / 4.0;
      const double RC1 = 0.5 * (C_0 + C_1);
      const double RC2 = 0.5 * (C_0 - C_1);

      TArray2 V_ia;
      const char* i = "i";
      const char* a = "a";
      V_ia = VRk_Sk(i,a, C_0, C_1);

      L1_f12contri("a,i") = Delta_ai("a,i") * V_ia("i,a");

      TArray4 rg_ijab;
      TArray4d r_ijaPn = ijxy("<i j|r|a' n>");
      TArray4d g_aPnab = ijxy("<a' n|g|a b>");
      rg_ijab("i,j,a,b") =  _4("<i j|gr|a b>")
                          - _4("<i j|r|p q>") * _4("<p q|g|a b>")
                          - r_ijaPn("i,j,a',n") * g_aPnab("a',n,a,b")
                          - r_ijaPn("j,i,a',n") * g_aPnab("a',n,b,a");

      TArray4d r_ijaFb = ijxy("<i j|r|a_F(c') b>");

      L2_f12contri("a,b,i,j") = Delta_abij("a,b,i,j")
                                * (  RC1 * r_ijaFb("i,j,a,b") + RC2 * r_ijaFb("j,i,a,b")
                                   + RC1 * r_ijaFb("j,i,b,a") + RC2 * r_ijaFb("i,j,b,a")
                                   + RC1 * rg_ijab("i,j,a,b") + RC2 * rg_ijab("j,i,a,b")
                                  );

      if (method == "f12" || method == "F12") {
        L1("a,i") = t1("a,i") + L1_f12contri("a,i") * 2.0;
        L2("a,b,i,j") = t2("a,b,i,j") + L2_f12contri("a,b,i,j") * 2.0;

      } else if (method == "f12b" || method == "F12b") {

        // contributions to L1 equation
        const char* j = "j";
        const char* k = "k";
        const char* l = "l";
        const char* b = "b";

        // 1st contribution
        Vijab = VPq_Rs(i, j, a, b, C_0, C_1);
        L1_f12b("a,i") = Delta_ai("a,i")
                         * (Vijab("i,j,a,b") * t1("b,j")); // 1/2 R_ij^AB g_AB^ab t^j_b

        // 2nd contribution that depends on L2 amplitudes
        Vklib = VPq_Rs(k, l, i, b, C_0, C_1);
        Vklij = VPq_Rs(k, l, i, j, C_0, C_1);
        //TArray4 V_KLAJ = VPQ_RS(k, l, a, j);
        //TArray4 V_IJKL = VPQ_RS(i, j, k, l);
//        L1_f12b_L2("a,i") = Delta_ai("a,i") * (
//                             // - 1/4 l^ab_kl R^kl_AB g^AB_ib
//                           - Vklib("k,l,i,b") * L2("a,b,k,l")
//
//                             // + 1/4 l^ab_kl R^kl_AB g^AB_ij t^j_b
//                           + Vklij("k,l,i,j") * L2("a,b,k,l") * t1("b,j")
//                           );

        L1("a,i") = t1("a,i")
                   + L1_f12contri("a,i")
                   + L1_f12b("a,i")// + L1_f12b_L2("a,i")
                   ;
        L2("a,b,i,j") =  t2("a,b,i,j")
                       + L2_f12contri("a,b,i,j")
                       ;
      }
    }

    TArray4 gabij_temp;
    gabij_temp("a,b,i,j") = 2.0 * g_abij("a,b,i,j") - g_abij("b,a,i,j");

    // pseudo energy
    double E_0_L = 0.0;
    double E_1_L = dot(gabij_temp("a,b,i,j"), L2("a,b,i,j") );
    double Delta_E_L = std::abs(E_0_L - E_1_L);

    double iter = 0;
    std::cout << indent << scprintf("%-5s", "Iter") << scprintf("%-20s", "Delta_E_L")
              << scprintf("%-20s", "E_L(CCSD)") << std::endl;
    std::cout << indent << scprintf("%-5.0f", iter) << scprintf("%-15.12f", Delta_E_L)
              << scprintf("%-15.12f", E_1_L) << std::endl;

    // compute intermediates
    TArray2 CFkc, CFac, CFki;
    // \cal{W}mbej
    TArray4 CW_KbCj_ab,
            ////// CCSD only term
            CW_KbcJ_ba,
    // \cal{W}abei
            CW_AbCd_ab, CW_AbCi_ab,
    // \cal{W}mbij
            CW_KlIj_ab, CW_KbIj_ab,
    // intermediates for L2
            CW_KlIc_ab, CW_KliC_ab, CW_AkCd_ab;

    compute_intermediates_CW_ccsd(t1, t2, CFkc, CFac, CFki,
                                  CW_KbCj_ab, CW_KbcJ_ba,
                                  CW_AbCd_ab, CW_AbCi_ab,
                                  CW_KlIj_ab, CW_KbIj_ab,
                                  CW_KlIc_ab, CW_KliC_ab, CW_AkCd_ab);

    TArray2 CGac, CGki;

    TArray4 L2_unsymm, L2abij_temp;

    while (Delta_E_L >= 1.0e-12) {

      L2abij_temp("a,b,i,j") = 2.0 * L2("a,b,i,j") - L2("b,a,i,j");

      CGac("a,c") =  // - 1/2 t^cd_kl lambda^kl_ad
                   - t2("c,d,k,l") * L2abij_temp("a,d,k,l");
      CGki("k,i") =  // 1/2 t^cd_kl lambda^il_cd
                   t2("c,d,k,l") * L2abij_temp("c,d,i,l");

      auto L1_temp =  Delta_ai("a,i") * (
                      //   \cal{F}ia
                      CFkc("i,a")
                      // + \lambda^i_c \cal{F}ca
                    + L1("c,i") * CFac("c,a")
                      // - \lambda^k_a \cal{F}ik
                    - L1("a,k") * CFki("i,k")
                      // + \lambda^k_c \cal{W}icak
                    + L1("c,k") * (2.0 * CW_KbCj_ab("i,c,a,k") + CW_KbcJ_ba("i,c,a,k"))
                      // + 1/2 \lambda^ik_cd \cal{W}cdak
                    + L2abij_temp("c,d,i,k") * CW_AbCi_ab("c,d,a,k")
                      // - 1/2 \lambda^kl_ac \cal{W}ickl
                    - L2abij_temp("a,c,k,l") * CW_KbIj_ab("i,c,k,l")

                      // CCSD only terms:
                      // - cal{G}cd cal{W}cida
                    - CGac("c,d") * (2.0 * CW_AkCd_ab("c,i,d,a") - CW_AkCd_ab("c,i,a,d")) // CW_aKcD => CW_AkCd_ab
                      // - cal{G}kl cal{W}kila
                    - CGki("k,l") * (2.0 * CW_KlIc_ab("k,i,l,a") + CW_KliC_ab("k,i,l,a"))
                     );

      auto L2_temp =  Delta_abij("a,b,i,j") * (
                              //   g^ij_ab
                              g_abij("a,b,i,j")

                              // + P(ij) \lambda^i_c \cal{W}cjab
                            + L1("c,i") * CW_AkCd_ab("c,j,a,b")
                            + L1("c,j") * CW_AkCd_ab("c,i,b,a")

                              // - P(ab) \lambda^k_a \cal{W}ijkb
                            - L1("a,k") * CW_KlIc_ab("i,j,k,b")
                            + L1("b,k") * CW_KliC_ab("i,j,k,a")

                              // + P(ij) P(ab) \lambda^i_a \cal{F}jb
                            + L1("a,i") * CFkc("j,b")
                              //- L1("b,i") * CFkc("j,a")
                              //- L1("a,j") * CFkc("i,b")
                            + L1("b,j") * CFkc("i,a")

                              // CCSD term
                              // + P(ab) \lambda^ij_ac \cal{F}cb
                            + L2("a,c,i,j") * CFac("c,b")
                            + L2("c,b,i,j") * CFac("c,a")

                              // CCSD term
                              // - P(ij) \lambda^ik_ab \cal{F}jk
                            - L2("a,b,i,k") * CFki("j,k")
                            - L2("a,b,k,j") * CFki("i,k")

                              ////// CCSD only terms
                              // + 1/2 \lambda^kl_ab \cal{W}ijkl
                            + L2("a,b,k,l") * CW_KlIj_ab("i,j,k,l")

                              // + 1/2 \lambda^ij_cd \cal{W}cdab
                            + L2("c,d,i,j") * CW_AbCd_ab("c,d,a,b")

                              // + P(ij) P(ab) \lambda^ik_ac \cal{W}jcbk
                            + (  L2abij_temp("a,c,i,k") * CW_KbCj_ab("j,c,b,k")
                               + L2("a,c,i,k") * CW_KbcJ_ba("j,c,b,k")) // CW_kbcj => CW_KBCJ_aa
                              // P(ab)
                            + L2("c,b,i,k") * CW_KbcJ_ba("j,c,a,k") // CW_kBCj => CW_KbcJ_ab
                              // P(ij)
                            + L2("a,c,k,j") * CW_KbcJ_ba("i,c,b,k")
                              // P(ij) P(ab)
                            + (  L2abij_temp("b,c,j,k") * CW_KbCj_ab("i,c,a,k")
                               + L2("b,c,j,k") * CW_KbcJ_ba("i,c,a,k"))

                            // + P(ab) g^ij_ac Gbc
                            + g_abij("a,c,i,j") * CGac("b,c")
                            + g_abij("c,b,i,j") * CGac("a,c")

                            // - P(ij) g^ik_ab Gkj
                            - g_abij("a,b,i,k") * CGki("k,j")
                            - g_abij("a,b,k,j") * CGki("k,i")
                            );

      if (method == "ccsd" || method == "CCSD") {
        L1("a,i") = L1_temp;
        L2_unsymm("a,b,i,j") = L2_temp;

      } else if (method == "f12" || method == "F12") {
        L1("a,i") = L1_temp + L1_f12contri("a,i") * 2.0;
        L2_unsymm("a,b,i,j") = L2_temp + L2_f12contri("a,b,i,j") * 2.0;

      } else if (method == "f12b" || method == "F12b") {

        L1_f12b_L2("a,i") = Delta_ai("a,i") * (
                             // - 1/4 l^ab_kl R^kl_AB g^AB_ib
                           - Vklib("k,l,i,b") * L2("a,b,k,l")

                             // + 1/4 l^ab_kl R^kl_AB g^AB_ij t^j_b
                           + Vklij("k,l,i,j") * L2("a,b,k,l") * t1("b,j")
                           );

        L1("a,i") =  L1_temp
                   + L1_f12contri("a,i")
                   + L1_f12b("a,i") + L1_f12b_L2("a,i")
                   ;
        L2_unsymm("a,b,i,j") =  L2_temp
                              + L2_f12contri("a,b,i,j")
                              ;

      }

      L2("a,b,i,j") = 0.5 * (L2_unsymm("a,b,i,j") + L2_unsymm("b,a,j,i"));

      E_0_L = E_1_L;
      E_1_L = dot(gabij_temp("a,b,i,j"), L2("a,b,i,j") );
      Delta_E_L = std::abs(E_0_L - E_1_L);

      iter += 1;
      std::cout << indent << scprintf("%-5.0f", iter) << scprintf("%-15.12f", Delta_E_L)
                << scprintf("%-15.12f", E_1_L) << std::endl;
    }
    //std::cout << "CC2 L1 amplitudes: " << std::endl << L1 << std::endl;
  }

  // compute intermediates for CCSD Xam (Xai)
  template <typename T>
  void SingleReference_R12Intermediates<T>::compute_Gamma_ijab_ccsd(const TArray2& T1, const TArray4& T2,
                                                                    const TArray4& tau_ab,
                                                                    const TArray2& L1, const TArray4& L2,
                                                                    TArray4& Gamma_IjAb_ab) {
//    // i,j,a,b in alpha or beta space
//     Gamma_IJAB_aa("i,j,a,b") =  //   1/4 tau^ab_ij
//                                 0.25 * tau_aa("a,b,i,j")
//
//                                 // CCSD modified term
//                                 // + 1/16 tau^cd_ij lambda^kl_cd tau^ab_kl
//                               + 0.0625 * tau_aa("c,d,i,j")
//                                        * (L2("c,d,k,l") - L2("c,d,l,k"))
//                                        * tau_aa("a,b,k,l")
//
//                                 // CCSD terms only
//                                 // - 1/8 P-(ij) t^cd_il lambda^kl_cd tau^ab_kj
//                               - 0.125 * (T2("c,d,i,l") - T2("c,d,l,i"))
//                                       * (L2("c,d,k,l") - L2("c,d,l,k")) * tau_aa("a,b,k,j")
//                               - 0.25 * T2("c,d,i,l") * L2("c,d,k,l") * tau_aa("a,b,k,j")
//                                 // - P(ij)
//                               + 0.125 * (T2("c,d,j,l") - T2("c,d,l,j"))
//                                       * (L2("c,d,k,l") - L2("c,d,l,k")) * tau_aa("a,b,k,i")
//                               + 0.25 * T2("c,d,j,l") * L2("c,d,k,l") * tau_aa("a,b,k,i")
//
//
//                                 // - 1/4 P-(ij) t^c_i lambda^k_c tau^ab_kj
//                               - 0.25  * (  T1("c,i") * tau_aa("a,b,k,j")
//                                          - T1("c,j") * tau_aa("a,b,k,i")
//                                          ) * L1("c,k")
//
//                                 // CCSD only term
//                                 // - 1/8 P-(ab) t^ad_kl lambda^kl_cd tau^cb_ij
//                               - 0.125 * (T2("a,d,k,l") - T2("d,a,k,l"))
//                                       * (L2("c,d,k,l") - L2("c,d,l,k")) * tau_aa("c,b,i,j")
//                               - 0.25 * T2("a,d,k,l") * L2("c,d,k,l") * tau_aa("c,b,i,j")
//                                 // - P(ab)
//                               + 0.125 * (T2("b,d,k,l") - T2("d,b,k,l"))
//                                       * (L2("c,d,k,l") - L2("c,d,l,k")) * tau_aa("c,a,i,j")
//                               + 0.25 * T2("b,d,k,l") * L2("c,d,k,l") * tau_aa("c,a,i,j")
//
//                                 // - 1/4 P-(ab) t^a_k lambda^k_c tau^cb_ij
//                               - 0.25 * (  T1("a,k") * tau_aa("c,b,i,j")
//                                         - T1("b,k") * tau_aa("c,a,i,j")
//                                        ) * L1("c,k")
//
//                                 // CCSD only term
//                                 // - 1/8 P-(ij) P-(ab) (t^ac_ki + 2 t^c_i t^a_k) lambda^kl_cd t^bd_jl
//                                 // + 1/8 P-(ij) P-(ab) t^ac_ik lambda^kl_cd t^bd_jl
//                               + 0.125 * (  (T2("a,c,i,k") - T2("c,a,i,k"))
//                                            * (L2("c,d,k,l") - L2("d,c,k,l"))
//                                            * (T2("b,d,j,l") - T2("d,b,j,l"))
//                                          + (T2("a,c,i,k") - T2("c,a,i,k")) * L2("c,d,k,l") * T2("b,d,j,l")
//                                          + T2("a,c,i,k") * (L2("c,d,k,l") - L2("d,c,k,l")) * T2("b,d,j,l")
//                                          + T2("a,c,i,k") * L2("c,d,k,l") * (T2("b,d,j,l") - T2("d,b,j,l"))
//
//                                            // P(ab)
//                                          - (T2("b,c,i,k") - T2("c,b,i,k"))
//                                            * (L2("c,d,k,l") - L2("d,c,k,l"))
//                                            * (T2("a,d,j,l") - T2("d,a,j,l"))
//                                          - (T2("b,c,i,k") - T2("c,b,i,k")) * L2("c,d,k,l") * T2("a,d,j,l")
//                                          - T2("b,c,i,k") * (L2("c,d,k,l") - L2("d,c,k,l")) * T2("a,d,j,l")
//                                          - T2("b,c,i,k") * L2("c,d,k,l") * (T2("a,d,j,l") - T2("d,a,j,l"))
//
//                                            // P(ij)
//                                          - (T2("a,c,j,k") - T2("c,a,j,k"))
//                                            * (L2("c,d,k,l") - L2("d,c,k,l"))
//                                            * (T2("b,d,i,l") - T2("d,b,i,l"))
//                                          - (T2("a,c,j,k") - T2("c,a,j,k")) * L2("c,d,k,l") * T2("b,d,i,l")
//                                          - T2("a,c,j,k") * (L2("c,d,k,l") - L2("d,c,k,l")) * T2("b,d,i,l")
//                                          - T2("a,c,j,k") * L2("c,d,k,l") * (T2("b,d,i,l") - T2("d,b,i,l"))
//
//                                            // P(ij)P(ab)
//                                          + (T2("b,c,j,k") - T2("c,b,j,k"))
//                                            * (L2("c,d,k,l") - L2("d,c,k,l"))
//                                            * (T2("a,d,i,l") - T2("d,a,i,l"))
//                                          + (T2("b,c,j,k") - T2("c,b,j,k")) * L2("c,d,k,l") * T2("a,d,i,l")
//                                          + T2("b,c,j,k") * (L2("c,d,k,l") - L2("d,c,k,l")) * T2("a,d,i,l")
//                                          + T2("b,c,j,k") * L2("c,d,k,l") * (T2("a,d,i,l") - T2("d,a,i,l"))
//                                         )
//                                 // - 1/4 P-(ij) P-(ab) t^c_i t^a_k lambda^kl_cd t^bd_jl
//                               - 0.25 * ( (  T1("c,i") * T1("a,k") * (T2("b,d,j,l") - T2("d,b,j,l"))
//                                           - T1("c,i") * T1("b,k") * (T2("a,d,j,l") - T2("d,a,j,l"))
//                                           - T1("c,j") * T1("a,k") * (T2("b,d,i,l") - T2("d,b,i,l"))
//                                           + T1("c,j") * T1("b,k") * (T2("a,d,i,l") - T2("d,a,i,l"))
//                                           ) * (L2("c,d,k,l") - L2("d,c,k,l"))
//
//                                         + (  T1("c,i") * T1("a,k") * T2("b,d,j,l")
//                                            - T1("c,i") * T1("b,k") * T2("a,d,j,l")
//                                            - T1("c,j") * T1("a,k") * T2("b,d,i,l")
//                                            + T1("c,j") * T1("b,k") * T2("a,d,i,l")
//                                           ) * L2("c,d,k,l")
//                                        )
//
//                                 // - 1/4 P-(ij) P-(ab) (t^ac_ki + 2 t^c_i t^a_k) lambda^k_c t^b_j
//                                 // => - 1/4 P-(ij) P-(ab) (- t^ac_ik + 2 t^c_i t^a_k) t^b_j
//                               - 0.25 * (  (- 2.0 * T2("a,c,i,k") + T2("c,a,i,k")
//                                            + 2.0 * T1("c,i") * T1("a,k")
//                                           ) * T1("b,j")
//                                           // P(ab)
//                                         - (- 2.0 * T2("b,c,i,k") + T2("c,b,i,k")
//                                            + 2.0 * T1("c,i") * T1("b,k")
//                                           ) * T1("a,j")
//                                           // P(ij)
//                                         - (- 2.0 * T2("a,c,j,k") + T2("c,a,j,k")
//                                            + 2.0 * T1("c,j") * T1("a,k")
//                                           ) * T1("b,i")
//                                           // P(ij) P(ab)
//                                         + (- 2.0 * T2("b,c,j,k") + T2("c,b,j,k")
//                                            + 2.0 * T1("c,j") * T1("b,k")
//                                           ) * T1("a,i")
//                                        ) * L1("c,k")
//
//                                 // + 3/4 P-(ij) P-(ab) t^a_i t^c_j lambda^k_c t^b_k
//                               + 0.75 * (  T1("a,i") * T1("c,j") * T1("b,k")
//                                         - T1("b,i") * T1("c,j") * T1("a,k")
//                                         - T1("a,j") * T1("c,i") * T1("b,k")
//                                         + T1("b,j") * T1("c,i") * T1("a,k")
//                                        ) * L1("c,k")
//
//                               ;

     // a,i in alpha space, b,j in beta space
     Gamma_IjAb_ab("i,j,a,b") =  //   1/4 tau^ab_ij
                                 0.25 * tau_ab("a,b,i,j")

                                 // CCSD modified term
                                 // + 1/16 tau^cd_ij lambda^kl_cd tau^ab_kl
                               + 0.25 * tau_ab("c,d,i,j") * L2("c,d,k,l")
                                      * tau_ab("a,b,k,l")

                                 // CCSD terms only
                                 // - 1/8 P-(ij) t^cd_il lambda^kl_cd tau^ab_kj
                               - 0.125 * (T2("c,d,i,l") - T2("c,d,l,i"))
                                       * (L2("c,d,k,l") - L2("c,d,l,k")) * tau_ab("a,b,k,j")
                               - 0.25 * T2("c,d,i,l") * L2("c,d,k,l") * tau_ab("a,b,k,j")
                                 // - P(ij)
                               - 0.125 * (T2("c,d,j,l") - T2("c,d,l,j"))
                                       * (L2("c,d,k,l") - L2("c,d,l,k")) * tau_ab("a,b,i,k")
                               - 0.25 * T2("c,d,j,l") * L2("c,d,k,l") * tau_ab("a,b,i,k")

                                 // - 1/4 P-(ij) t^c_i lambda^k_c tau^ab_kj
                               - 0.25 * (  T1("c,i") * tau_ab("a,b,k,j")
                                         + T1("c,j") * tau_ab("a,b,i,k")
                                        ) * L1("c,k")

                                  // CCSD only term
                                  // - 1/8 P-(ab) t^ad_kl lambda^kl_cd tau^cb_ij
                                - 0.125 * (T2("a,d,k,l") - T2("d,a,k,l"))
                                        * (L2("c,d,k,l") - L2("c,d,l,k")) * tau_ab("c,b,i,j")
                                - 0.25 * T2("a,d,k,l") * L2("c,d,k,l") * tau_ab("c,b,i,j")
                                  // - P(ab)
                                - 0.125 * (T2("b,d,k,l") - T2("d,b,k,l"))
                                        * (L2("c,d,k,l") - L2("c,d,l,k")) * tau_ab("a,c,i,j")
                                - 0.25 * T2("b,d,k,l") * L2("c,d,k,l") * tau_ab("a,c,i,j")

                                  // - 1/4 P-(ab) t^a_k lambda^k_c tau^cb_ij
                                - 0.25 *  (  T1("a,k") * tau_ab("c,b,i,j")
                                           + T1("b,k") * tau_ab("a,c,i,j")
                                          ) * L1("c,k")

                                  // - 1/8 P-(ij) P-(ab) (t^ac_ki + 2 t^c_i t^a_k) lambda^kl_cd t^bd_jl
                                  // + 1/8 P-(ij) P-(ab) t^ac_ik lambda^kl_cd t^bd_jl
                                + 0.125 * (  (T2("a,c,i,k") - T2("c,a,i,k"))
                                             * (L2("c,d,k,l") - L2("d,c,k,l")) * T2("b,d,j,l")
                                           + (T2("a,c,i,k") - T2("c,a,i,k"))
                                             * L2("c,d,k,l") * (T2("b,d,j,l") - T2("d,b,j,l"))
                                           + T2("a,c,i,k") * (L2("c,d,k,l") - L2("d,c,k,l"))
                                             * (T2("b,d,j,l") - T2("d,b,j,l"))
                                           + T2("a,c,i,k") * L2("c,d,k,l") * T2("b,d,j,l")

                                             // P(ab)
                                           + T2("c,b,i,k") * L2("c,d,l,k") * T2("a,d,l,j")

                                             // P(ij)
                                           + T2("a,c,k,j") * L2("d,c,k,l") * T2("d,b,i,l")

                                             //P(ij)P(ab)
                                           + T2("b,c,j,k") * (L2("c,d,k,l") - L2("d,c,k,l"))
                                             * (T2("a,d,i,l") - T2("d,a,i,l"))
                                           + T2("b,c,j,k") * L2("c,d,k,l") * T2("a,d,i,l")
                                           + (T2("b,c,j,k") - T2("c,b,j,k"))
                                             * L2("c,d,k,l") * (T2("a,d,i,l") - T2("d,a,i,l"))
                                           + (T2("b,c,j,k") - T2("c,b,j,k"))
                                             * (L2("c,d,k,l") - L2("d,c,k,l"))
                                             * T2("a,d,i,l")
                                          )
                                  // - 1/4 P-(ij) P-(ab) t^c_i t^a_k lambda^kl_cd t^bd_jl
                                - 0.25 * (  T1("c,i") * T1("a,k")
                                            * ( (L2("c,d,k,l") - L2("d,c,k,l")) * T2("b,d,j,l")
                                               + L2("c,d,k,l") * (T2("b,d,j,l") - T2("d,b,j,l"))
                                              )
                                            // P(ab)
                                          - T1("c,i") * T1("b,k") * L2("c,d,l,k") * T2("a,d,l,j")
                                            // P(ij)
                                          - T1("c,j") * T1("a,k") * L2("d,c,k,l") * T2("d,b,i,l")
                                            // P(ij) P(ab)
                                          + T1("c,j") * T1("b,k")
                                            * (  (L2("c,d,k,l") - L2("d,c,k,l")) * T2("a,d,i,l")
                                               + L2("c,d,k,l") * (T2("a,d,i,l") - T2("d,a,i,l"))
                                              )
                                         )

                                  // - 1/4 P-(ij) P-(ab) (t^ac_ki + 2 t^c_i t^a_k) lambda^k_c t^b_j
                                - 0.25 * ( (- 2.0 * T2("a,c,i,k") + T2("c,a,i,k")
                                            + 2.0 * T1("c,i") * T1("a,k")
                                           )  * T1("b,j")
                                            //
                                          + (- 2.0 * T2("b,c,j,k") + T2("c,b,j,k")
                                             + 2.0 * T1("c,j") * T1("b,k")
                                            )  * T1("a,i")
                                         ) * L1("c,k")

                                  // + 3/4 P-(ij) P-(ab) t^a_i t^c_j lambda^k_c t^b_k
                                + 0.75 * (  T1("a,i") * T1("c,j") * T1("b,k")
                                          + T1("b,j") * T1("c,i") * T1("a,k")
                                         ) * L1("c,k")
                                ;

  }

  // compute Gamma intermediates following JCP, 95, 2623(1991)
  template <typename T>
  void SingleReference_R12Intermediates<T>::compute_Gamma2_ccsd(const TArray2& T1, const TArray4& T2,
                                                               const TArray2& L1, const TArray4& L2,
                                                               const TArray4& tau_aa, const TArray4& tau_ab,
                                                               TArray4& Gamma_IjKa_ab,
                                                               TArray4& Gamma_AbCi_ab,
                                                               TArray4& Gamma_iBjA_ab, TArray4& Gamma_iBJa_ba,
                                                               TArray4& Gamma_AbCd_ab, TArray4& Gamma_IjKl_ab,
                                                               TArray4& Gamma_IjAb_ab) {

    TArray2 CGac, CGki;
    CGac("a,c") = - 0.5 * (T2("c,d,k,l") - T2("d,c,k,l"))
                        * (L2("a,d,k,l") - L2("d,a,k,l"))
                  - T2("c,d,k,l") * L2("a,d,k,l");

    CGki("k,i") =  0.5 * (T2("c,d,k,l") - T2("c,d,l,k"))
                       * (L2("c,d,i,l") - L2("c,d,l,i"))
                 + T2("c,d,k,l") * L2("c,d,i,l");

    TArray4 VKlIj_ab, VIBJA_aa, VIbJa_ab, ViBJa_ba, VAbCd_ab;

    VKlIj_ab("k,l,i,j") = tau_ab("c,d,k,l") * L2("c,d,i,j");

    VIBJA_aa("i,b,j,a") =  0.5 * (T2("a,c,i,k") - T2("a,c,k,i"))
                               * (L2("b,c,j,k") - L2("b,c,k,j"))
                         + 0.5 * T2("a,c,i,k") * L2("b,c,j,k");
    VIbJa_ab("i,b,j,a") = 0.5 * T2("c,a,i,k") * L2("c,b,j,k");
    ViBJa_ba("i,b,j,a") =  0.5 * (T2("a,c,i,k") - T2("a,c,k,i")) * L2("b,c,j,k")
                         + 0.5 * T2("a,c,i,k") * (L2("b,c,j,k") - L2("b,c,k,j"));
    // VIbjA_ba("i,b,j,a") = ViBJa_ba("i,b,j,a")
    // ViBjA_ab("i,b,j,a") = VIbJa_ab("i,b,j,a")
    // ViBjA_ab("i,b,j,a") = 0.5 * T2("a,c,k,i") * L2("b,c,k,j");
    // VIBja_ba("i,b,j,a") = VIbJa_ab("i,b,j,a")
    //VIBja_ba("i,b,j,a") = 0.5 * T2("c,a,i,k") * L2("b,c,k,j");

    VAbCd_ab("a,b,c,d") = tau_ab("c,d,k,l") * L2("a,b,k,l");
    // V_aBCd = - tau_ab("d,c,k,l") * L2("a,b,k,l")

    Gamma_IjKa_ab("i,j,k,a") =  // - 1/8 tau^Ca_Ij lambda^K_C
                              - 0.125 * tau_ab("c,a,i,j") * L1("c,k")

                               // - 1/8 lambda^Ca_Ij t^K_C
                             - 0.125 * L2("c,a,i,j") * T1("c,k")

                                // + 1/8 V_IjKl t^a_l
                              + 0.125 * VKlIj_ab("i,j,k,l") * T1("a,l")

                                //*** + 1/4 P-(ij) V_kcia (V_icka ???) t^c_j
                                // + 1/4 V_IcKa t^c_j
                              + 0.25 * VIbJa_ab("i,c,k,a") * T1("c,j")
                                // - 1/4 V_jCKa t^C_I
                              - 0.25 * ViBJa_ba("j,c,k,a") * T1("c,i")

                                // - 1/8 P-(ij) CG_IK t^a_j
                              - 0.125 * CGki("i,k") * T1("a,j")
                              ;

//    Gamma_iJkA_ab("i,j,k,a") =  // - 1/8 tau^cA_iJ lambda^k_c
//                              - 0.125 * tau_ab("c,a,i,j") * L1("c,k")
//
//                               // - 1/8 lambda^cA_iJ t^k_c
//                             - 0.125 * L2("c,a,i,j") * T1("c,k")
//
//                                // + 1/8 V_iJkL t^A_L
//                              + 0.125 * VKlIj_ab("i,j,k,l") * T1("a,l")
//
////                                // + 1/4 P-(ij) V_kCiA t^C_J
////                                // - 1/4 V_kcJA t^c_i
//
//                                // - 1/8 P-(ij) CG_ik t^A_J
//                              - 0.125 * CGki("i,k") * T1("a,j")
//                              ;

    Gamma_AbCi_ab("a,b,c,i") =  //   1/8 tau^ab_ki lambda^k_c
                                0.125 * tau_ab("a,b,k,i") * L1("c,k")

                                //   1/8 lambda^Ki_Ab t^C_K
                              + 0.125 * L2("a,b,k,i") * T1("c,k")

                                // - 1/8 t^d_i V_CdAb
                              - 0.125 * T1("d,i") * VAbCd_ab("c,d,a,b")

                                //*** + 1/4 P-(ab) V_kcib (V_ickb ???) t^a_k
                                // + 1/4 V_iCKb t^A_K
                              + 0.25 * ViBJa_ba("i,c,k,b") * T1("a,k")
                                // - 1/4 V_iCkA t^b_k
                              - 0.25 * VIbJa_ab("i,c,k,a") * T1("b,k")

                                // - 1/8 P-(ab) t^b_i CG_ca
                                // - 1/8 t^b_i CG_CA
                              - 0.125 * CGac("c,a") * T1("b,i")
                              ;

    // Gamma_AbcI_ab("a,b,c,i") = - Gamma_AbcI_ab("b,a,c,i")
//    Gamma_AbcI_ab("a,b,c,i") =  //   1/8 tau^Ab_kI lambda^k_c
//                              - 0.125 * tau_ab("a,b,i,k") * L1("c,k") // tau_ab("b,a,k,i")
//
//                                //   1/8 lambda^kI_Ab t^c_k
//                              - 0.125 * L2("b,a,k,i") * T1("c,k")
//
//                                // - 1/8 t^D_I V_cDAb (- V_DcAb or - V_cDbA)
//                              + 0.125 * T1("d,i") * VAbCd_ab("c,d,b,a")
//
//                                // + 1/4 P-(ab) V_kcib t^a_k
//                                // + 1/4 V_KCib t^A_K
//                                // - 1/4 V_kCiA t^b_k
//
//                                // - 1/8 P-(ab) t^b_i CG_ca
//                                // + 1/8 t^A_I CG_cb
//                              + 0.125 * CGac("c,b") * T1("a,i")
//                              ;

    // Gamma_IBJA_aa = Gamma_IbJa_ab + Gamma_iBJa_ab

    Gamma_iBJa_ba("i,b,j,a") =   // -?? 1/4 P+(ia,jb) Vibja
                                 // - 1/4 ViBJa
                               - 0.25 * ViBJa_ba("i,b,j,a")
                                 //  1/4 VJaiB (ViBJa)
                               - 0.25 * ViBJa_ba("j,a,i,b")

                                 // + 1/8 P+(ia,jb) t^c_i t^a_k lambda^jk_bc
                                 // + 1/8 t^c_i t^a_k lambda^Jk_Bc
                               + 0.125 * T1("c,i") * T1("a,k") * L2("b,c,j,k")
                                 // + 1/8 t^C_J t^B_K lambda^iK_aC
                               + 0.125 * T1("c,j") * T1("b,k") * L2("a,c,i,k")

                                 // - 1/8 P+(ia,jb) t^a_i lambda^j_b
                                 // - 1/8 t^a_i lambda^J_B
                               - 0.125 * L1("b,j") * T1("a,i")
                                 // - 1/8 t^B_J lambda^i_a
                               - 0.125 * L1("a,i") * T1("b,j")
                               ;

    // Gamma_iBjA_ab = Gamma_IbJa_ab
    Gamma_iBjA_ab("i,b,j,a") =  //  - ?? 1/4 P+(ia,jb) Vibja
                                // - 1/4 ViBjA (VIbJa)
                              - 0.25 * VIbJa_ab("i,b,j,a")
                                // - 1/4 VjAiB
                              - 0.25 * VIbJa_ab("j,a,i,b")

                                // + 1/8 P+(ia,jb) t^c_i t^a_k lambda^jk_bc
                                // + 1/8 t^c_i t^A_K lambda^jK_Bc
                              -  0.125 * T1("c,i") * T1("a,k") * L2("c,b,j,k")
                                // + 1/8 t^c_j t^B_K lambda^iK_Ac
                              - 0.125 * T1("c,j") * T1("b,k") * L2("c,a,i,k")

                                // - 1/8 P+(ia,jb) t^a_i lambda^j_b
                              ;

    // Gamma_IJKL_aa("i,j,k,l") = Gamma_IjKl_ab("i,j,k,l") - Gamma_IjKl_ab("j,i,k,l")
    Gamma_IjKl_ab("i,j,k,l") =  // 1/8 P+(ij,kl) Vijkl
                                0.125 * VKlIj_ab("k,l,i,j")
                              + 0.125 * VKlIj_ab("i,j,k,l");

    // Gamma_ABCD_ab("a,b,c,d") = Gamma_AbCd_ab("a,b,c,d") - Gamma_AbCd_ab("b,a,c,d")
    Gamma_AbCd_ab("a,b,c,d") =  // 1/8 P+(ij,kl) Vabcd
                                0.125 * VAbCd_ab("a,b,c,d")
                              + 0.125 * VAbCd_ab("c,d,a,b");

    Gamma_IjAb_ab("i,j,a,b") =  //   1/8 tau^Ab_Ij
                                0.125 * tau_ab("a,b,i,j")

                                //   1/8 lambda^Ab_Ij
                              + 0.125 * L2("a,b,i,j")

                                // + 1/16 VIjkl tau^Ab_kl
                              + 0.125 * VKlIj_ab("i,j,k,l") * tau_ab("a,b,k,l")

                                // - 1/8 P-(ij) (Gik + t^c_i lambda^k_c) tau^ab_kj
                                // - 1/8 (GIk + t^c_i lambda^k_c) tau^Ab_kj
                              - 0.125 * (CGki("i,k") + T1("c,i") * L1("c,k")) * tau_ab("a,b,k,j")
                                // + 1/8 (Gjk + t^c_j lambda^k_c) tau^Ab_kI
                              - 0.125 * (CGki("j,k") + T1("c,j") * L1("c,k")) * tau_ab("a,b,i,k")

                                // + 1/8 P-(ab) (Gca - t^a_k lambda^k_c) tau^cb_ij
                                // + 1/8 (GCA - t^A_K lambda^K_C) tau^Cb_Ij
                              + 0.125 * (CGac("c,a") - T1("a,k") * L1("c,k")) * tau_ab("c,b,i,j")
                                // - 1/8 (Gcb - t^b_k lambda^k_c) tau^cA_Ij
                              + 0.125 * (CGac("c,b") - T1("b,k") * L1("c,k")) * tau_ab("a,c,i,j")

                                //*** - 1/8 P-(ij) P-(ab) (t^ac_ki + 2 t^c_i t^a_k) Vjckb
                                // - 1/8 t^Ac_kI Vjkcb = + 1/8 t^AC_IK VjCKb + 1/8 t^Ac_Ik Vjckb
                              + 0.125 * (T2("a,c,i,k") - T2("c,a,i,k")) * ViBJa_ba("j,c,k,b")
                              + 0.125 * T2("a,c,i,k") * VIBJA_aa("j,c,k,b")
                                // + 1/8 t^bC_kI VjCkA
                              + 0.125 * T2("b,c,k,i") * VIbJa_ab("j,c,k,a")
                                // + 1/8 t^Ac_Kj VIcKb
                              + 0.125 * T2("a,c,k,j") * VIbJa_ab("i,c,k,b")
                                // - 1/8 t^bc_kj VIkcA = + 1/8 t^bC_jK VICKA + 1/8 t^bc_jk VIckA
                              + 0.125 * T2("b,c,j,k") * VIBJA_aa("i,c,k,a")
                              + 0.125 * (T2("b,c,j,k") - T2("c,b,j,k")) * ViBJa_ba("i,c,k,a")

                              //*** - 1/4 t^C_I t^A_K VjCKb
                              - 0.25 * T1("c,i") * T1("a,k") * ViBJa_ba("j,c,k,b")
                              // + 1/4 t^C_I t^b_k VjCkA
                              + 0.25 * T1("c,i") * T1("b,k") * VIbJa_ab("j,c,k,a")
                              // + 1/4 t^c_j t^A_K VIcKb
                              + 0.25 * T1("c,j") * T1("a,k") * VIbJa_ab("i,c,k,b")
                              // - 1/4 t^c_j t^b_k VIckA
                              - 0.25 * T1("c,j") * T1("b,k") * ViBJa_ba("i,c,k,a")


                                  // - 1/8 P-(ij) P-(ab) (t^ac_ki + 2 t^c_i t^a_k) lambda^k_c t^b_j
                              - 0.125 * (
                                           // (t^Ac_kI + 2 t^c_I t^A_k) lambda^k_c t^b_j
                                           (- 2.0 * T2("a,c,i,k") + T2("c,a,i,k")
                                            + 2.0 * T1("c,i") * T1("a,k")
                                           )  * T1("b,j")
                                           // (t^bc_kj + 2 t^c_j t^b_k) lambda^k_c t^A_I
                                         + (- 2.0 * T2("b,c,j,k") + T2("c,b,j,k")
                                            + 2.0 * T1("c,j") * T1("b,k")
                                           )  * T1("a,i")
                                        ) * L1("c,k")

                                // *** + 3/2 (3/8 ???) P-(ij) P-(ab) t^a_i t^c_j lambda^k_c t^b_k
                              + 0.375 * (  // t^A_I t^c_j lambda^k_c t^b_k
                                         T1("a,i") * T1("c,j") * T1("b,k")
                                         // t^b_j t^C_I lambda^K_C t^A_K
                                       + T1("b,j") * T1("c,i") * T1("a,k")
                                       ) * L1("c,k")
                              ;


  }

  // compute CCSD one-electron density from amplitudes
  template <typename T>
  void SingleReference_R12Intermediates<T>::compute_ccsd_1rdm_amp(const TArray2& T1, const TArray4& T2,
                                                                  const TArray2& L1, const TArray4& L2,
                                                                  TArray2& Dij, TArray2& Dab,
                                                                  TArray2& Dia, TArray2& Dai) {

    TArray4 T2_temp;
    T2_temp("a,b,i,j") = 2.0 *  T2("a,b,i,j") - T2("a,b,j,i");

    Dij("i,j") =  // - 1/4 P+(ij) t^cd_ik lambda^jk_cd
                - 0.5 * (  T2_temp("c,d,i,k") * L2("c,d,j,k")
                         + T2_temp("c,d,j,k") * L2("c,d,i,k")
                        )
                  // - 1/2 P+(ij) t^c_i lambda^j_c
                - 0.5 * (  T1("c,i") * L1("c,j")
                         + T1("c,j") * L1("c,i")
                        )
                ;

    Dab("a,b") =  //   1/4 P+(ab) t^kl_ac lambda^bc_kl
                  0.5 * (  T2_temp("a,c,k,l") * L2("b,c,k,l")
                         + T2_temp("b,c,k,l") * L2("a,c,k,l")
                        )
                  // + 1/2 P+(ab) lambda^k_a t^b_k
                + 0.5 * (  L1("a,k") * T1("b,k")
                         + L1("b,k") * T1("a,k")
                        )
                ;

    Dia("i,a") = //   t^a_i
                 T1("a,i")
                 // + (t^ac_ik - t^c_i t^a_k) \lambda^k_c
               + (T2_temp("a,c,i,k") - T1("c,i") * T1("a,k")) * L1("c,k")
                 // CCSD term
                 // - 1/2 lambda^kl_cd (t^cd_il t^a_k + t^c_i t^ad_kl)
               - (2.0 * L2("c,d,k,l") - L2("d,c,k,l"))
                 * (  T2("c,d,i,l") * T1("a,k")
                    + T2("a,d,k,l") * T1("c,i")
                   )
               ;

    Dai("a,i") = // lambda^i_a
                 L1("a,i");
  }

  // compute CCSD Xam (the right-hand side of Z-vector equations)
  // following formula from JCP, 103, 3561 (1995)
  template <typename T>
  void SingleReference_R12Intermediates<T>::compute_Xam_ccsd(const TArray2& T1, const TArray4& T2,
                                                             const TArray2& L1, const TArray4& L2,
                                                             TArray2& Xam_tot, TArray2& Xiip) {

    TArray4 tau_ab;
    tau_ab("a,b,i,j") = T2("a,b,i,j") + (T1("a,i") * T1("b,j"));

    // compute \Gamma(pq,rs) intermediates

    TArray4 Gamma_IjAb_ab;
    compute_Gamma_ijab_ccsd(T1, T2, tau_ab, L1, L2,
                            Gamma_IjAb_ab);

    TArray4 Gamma_AbIj_ab;
    Gamma_AbIj_ab("a,b,i,j") = // 1/4 lambda^ij_ab
                               0.25 * L2("a,b,i,j");
//    TArray4 Gamma_ABIJ_aa;
//    Gamma_ABIJ_aa("a,b,i,j") = // 1/4 lambda^ij_ab
//                               0.25 * (L2("a,b,i,j") - L2("b,a,i,j"));
//    TArray4 Gamma_aa_1, diff_1;
//    Gamma_aa_1("a,b,i,j") = Gamma_AbIj_ab("a,b,i,j") - Gamma_AbIj_ab("a,b,j,i");
//    diff_1("a,b,i,j") = Gamma_ABIJ_aa("a,b,i,j")  - Gamma_aa_1("a,b,i,j");
//    std::cout << std::endl << "Diff 1: " << std::endl
//                           << "norm = "  << static_cast<double>(diff_1("a,b,i,j").norm())
//                           << "  max = " << static_cast<double>(diff_1("a,b,i,j").max())
//                           << "  min = " << static_cast<double>(diff_1("a,b,i,j").min()) << std::endl;

    TArray4 Gamma_AbCd_ab;
    Gamma_AbCd_ab("a,b,c,d") = // 1/8 lambda^kl_ab tau^cd_kl (CCSD modified term)
                               0.25 * L2("a,b,k,l") * tau_ab("c,d,k,l");

    TArray4 Gamma_IjKl_ab;
    Gamma_IjKl_ab("i,j,k,l") = // 1/8 tau^cd_ij lambda^kl_cd
                               0.25 * tau_ab("c,d,i,j") * L2("c,d,k,l");
//    TArray4 Gamma_IJKL_aa;
//    Gamma_IJKL_aa("i,j,k,l") = // 1/8 tau^cd_ij lambda^kl_cd
//                               0.125 * tau_aa("c,d,i,j") * (L2("c,d,k,l") - L2("c,d,l,k"));
//    TArray4 Gamma_aa_1, diff_1, Gamma_aa_2, diff_2;
//    Gamma_aa_1("i,j,k,l") = Gamma_IjKl_ab("i,j,k,l") - Gamma_IjKl_ab("j,i,k,l");
//    diff_1("i,j,k,l") = Gamma_IJKL_aa("i,j,k,l")  - Gamma_aa_1("i,j,k,l");
//    std::cout << std::endl << "Diff 1: " << std::endl
//                           << "norm = "  << static_cast<double>(diff_1("i,j,k,l").norm())
//                           << "  max = " << static_cast<double>(diff_1("i,j,k,l").max())
//                           << "  min = " << static_cast<double>(diff_1("i,j,k,l").min()) << std::endl;
//
//    Gamma_aa_2("i,j,k,l") = Gamma_IjKl_ab("i,j,k,l") - Gamma_IjKl_ab("i,j,l,k");
//    diff_2("i,j,k,l") = Gamma_IJKL_aa("i,j,k,l")  - Gamma_aa_2("i,j,k,l");
//    std::cout << std::endl << "Diff 2: " << std::endl
//                           << "norm = "  << static_cast<double>(diff_2("i,j,k,l").norm())
//                           << "  max = " << static_cast<double>(diff_2("i,j,k,l").max())
//                           << "  min = " << static_cast<double>(diff_2("i,j,k,l").min()) << std::endl;

    TArray4 Gamma_IjKa_ab;
    Gamma_IjKa_ab("i,j,k,a") =  // - 1/4 tau^ca_ij lambda^k_c
                              - 0.25 * tau_ab("c,a,i,j") * L1("c,k")

                                // + 1/8 tau^cd_ij lambda^kl_cd t^a_l
                              + 0.25 * tau_ab("c,d,i,j") * L2("c,d,k,l") * T1("a,l")

                                // + 1/4 P-(ij) t^ac_il lambda^lk_cd t^d_j
                              + 0.25 * T2("a,c,l,i") * L2("c,d,k,l") * T1("d,j")
                                //
                              - 0.25 * (  (T2("a,c,j,l") - T2("a,c,l,j")) * L2("c,d,l,k")
                                        + T2("a,c,j,l") * (L2("c,d,l,k") - L2("c,d,k,l"))
                                       ) * T1("d,i")

                                // - 1/8 P-(ij) t^cd_il lambda^kl_cd t^a_j
                              - 0.125 * (  (T2("c,d,i,l") - T2("c,d,l,i")) * (L2("c,d,k,l") - L2("c,d,l,k"))
                                         + 2.0 * T2("c,d,i,l") * L2("c,d,k,l")
                                        ) * T1("a,j")
                                ;

    TArray4 Gamma_AKIJ_aa, Gamma_AkIj_ab;
    Gamma_AkIj_ab("a,k,i,j") =   // -1/4 t^c_k lambda^Ij_Ac
                               - 0.25 * T1("c,k") * L2("a,c,i,j");

    TArray4 Gamma_CiAb_ab;
    Gamma_CiAb_ab("c,i,a,b") =  //   1/4 tau^ab_ki lambda^k_c
                                0.25 * tau_ab("a,b,k,i") * L1("c,k")

                                // - 1/8 t^d_i lambda^kl_cd tau^ab_kl
                              - 0.25 * T1("d,i") * L2("c,d,k,l") * tau_ab("a,b,k,l")

                                // - 1/4 P-(ab) t^ad_il lambda^kl_cd t^b_k
                              - 0.25 * T2("a,d,l,i") * L2("c,d,l,k") * T1("b,k")
                                //
                              + 0.25 * (  (T2("b,d,i,l") - T2("b,d,l,i")) * L2("c,d,k,l")
                                        + T2("b,d,i,l") * (L2("c,d,k,l") - L2("c,d,l,k"))
                                       ) * T1("a,k")

                                // + 1/8 P-(ab) t^b_i t^ad_kl lambda^kl_cd
                              + 0.125 * T1("b,i") * (T2("a,d,k,l") - T2("d,a,k,l"))
                                       * (L2("c,d,k,l") - L2("c,d,l,k"))
                              + 0.25 * T1("b,i") * T2("a,d,k,l") * L2("c,d,k,l")
                              ;

    TArray4 Gamma_AbCi_ab;
    Gamma_AbCi_ab("a,b,c,i") = // 1/4 lambda^ki_ab t^c_k
                               0.25 * L2("a,b,k,i") * T1("c,k");

    TArray4 Gamma_aJiB_ab, Gamma_AjiB_ba;
    Gamma_aJiB_ab("a,j,i,b") =  // 1/4 lambda^ik_ac t^c_j t^cb_kj
                                0.25 * (L2("a,c,i,k") - L2("c,a,i,k")) * T2("c,b,k,j")
                              + 0.25 * L2("a,c,i,k") * (T2("c,b,k,j") - T2("b,c,k,j"))

                                // - 1/4 lambda^ik_ac t^c_j t^b_k
                              - 0.25 * L2("a,c,i,k") * T1("c,j") * T1("b,k")
                                // + 1/4 lambda^i_a t^b_j
                              + 0.25 * L1("a,i") * T1("b,j")
                              ;

    Gamma_AjiB_ba("a,j,i,b") =  // 1/4 lambda^ik_ac t^c_j t^cb_kj
                                0.25 * L2("c,a,i,k") * T2("b,c,k,j")

                                // - 1/4 lambda^ik_ac t^c_j t^b_k
                              + 0.25 * L2("a,c,k,i") * T1("c,j") * T1("b,k")
                                // + 1/4 lambda^i_a t^b_j
                              ;

    // compute CCSD density from amplitudes
    TArray2 Dij, Dab, Dia, Dai;
    compute_ccsd_1rdm_amp(T1, T2, L1, L2, Dij, Dab, Dia, Dai);

    TArray4d g_ijma = ijxy("<i j|g|m a>");
    TArray4d g_mabc = ijxy("<m a|g|b c>");
    TArray4d g_mabi = ijxy("<m a|g|b i>");
    TArray4d g_maib = ijxy("<m a|g|i b>");

    TArray4d g_imkl = ijxy("<i m|g|k l>");
    TArray4d g_miab = ijxy("<m i|g|a b>");
    TArray4d g_imbj = ijxy("<i m|g|b j>");
    TArray4d g_imjb = ijxy("<i m|g|j b>");

    TArray4d g_ijak = ijxy("<i j|g|a k>");
    TArray4d g_abci = ijxy("<a b|g|c i>");
    TArray4d g_aibj = ijxy("<a i|g|b j>");
    TArray4d g_iabj = ijxy("<i a|g|b j>");

    TArray4d g_ijab = ijxy("<i j|g|a b>");
    TArray4d g_abcd = ijxy("<a b|g|c d>");

    TArray2 Xam, Xai;
    Xam("a,m") =  (Dai("a,j") + Dia("j,a")) * _2("<j|F|m>")
                + (Dab("a,c") + Dab("c,a")) * _2("<c|F|m>")

                - (2.0 * g_mabc("m,c,a,d") - g_mabc("m,c,d,a"))
                  * (Dab("c,d") + Dab("d,c"))
                - (2.0 * g_imjb("l,m,k,a") - g_imbj("l,m,a,k"))
                  * (Dij("k,l") + Dij("l,k"))
                - (  2.0 * g_miab("m,j,a,b") - g_miab("m,j,b,a")
                   + 2.0 * g_mabi("m,b,a,j") - g_maib("m,b,j,a")
                  ) * (Dia("j,b") + Dai("b,j"))

                + 4.0 * (
                           //
                           // (Gamma(kl,ab) + Gamma(ab,kl)) g^kl_mb
                          (Gamma_IjAb_ab("k,l,a,b") + Gamma_AbIj_ab("a,b,k,l"))
                          * (g_ijma("k,l,m,b") * 2.0 - g_ijma("l,k,m,b"))

                           // (Gamma(cd,ab) + Gamma(ab, cd)) g^cd_mb
                         + (Gamma_AbCd_ab("c,d,a,b") + Gamma_AbCd_ab("a,b,c,d"))
                           * (g_mabc("m,b,c,d") * 2.0 - g_mabc("m,b,d,c"))

                           // 2 (Gamma(ck,ab) - Gamma(ba,ck)) g^ck_mb
                         + (  Gamma_CiAb_ab("c,k,a,b") - Gamma_CiAb_ab("c,k,b,a")
                            - Gamma_AbCi_ab("b,a,c,k") + Gamma_AbCi_ab("a,b,c,k"))
                           * (g_mabi("m,b,c,k") - g_maib("m,b,k,c"))
                         + (Gamma_CiAb_ab("c,k,a,b") + Gamma_AbCi_ab("a,b,c,k"))
                           * g_mabi("m,b,c,k")
                         - (- Gamma_CiAb_ab("c,k,b,a") - Gamma_AbCi_ab("b,a,c,k"))
                           * g_maib("m,b,k,c")

                           //
                           // (Gamma(lk, ja) + Gamma(aj, kl)) g^kl_mj
                         + (Gamma_IjKa_ab("l,k,j,a") + Gamma_AkIj_ab("a,j,k,l"))
                           * (g_imkl("j,m,l,k") * 2.0 - g_imkl("j,m,k,l"))

                           // (Gamma(cd, aj) + Gamma(aj, cd) g^mj_cd) g^cd_mj
                         + (Gamma_AbCi_ab("c,d,a,j") + Gamma_CiAb_ab("a,j,c,d"))
                           * (g_miab("m,j,c,d") * 2.0 - g_miab("m,j,d,c"))

                           // * (Gamma(ck, ja) + Gamma(aj,kc)) g^kc_mj
                         + (  Gamma_aJiB_ab("c,k,j,a") + Gamma_AjiB_ba("c,k,j,a")
                            + Gamma_aJiB_ab("a,j,k,c") + Gamma_AjiB_ba("a,j,k,c"))
                           * (g_imbj("j,m,c,k") - g_imjb("j,m,k,c"))
                         + (Gamma_aJiB_ab("c,k,j,a") + Gamma_aJiB_ab("a,j,k,c"))
                           * g_imbj("j,m,c,k")
                         - (Gamma_AjiB_ba("c,k,j,a") + Gamma_AjiB_ba("a,j,k,c"))
                           * g_imjb("j,m,k,c")
                        )
                ;

    Xai("a,i") = - _2("<a|F|b>") * (Dai("b,i") + Dia("i,b"))
                 - _2("<a|F|l>") * (Dij("l,i") + Dij("i,l"))

                 - 4.0 * (
                           //
                           // (Gamma(kl,ij) + Gamma(ij, kl)) g^kl_aj
                           (Gamma_IjKl_ab("k,l,i,j") + Gamma_IjKl_ab("i,j,k,l"))
                           * (g_ijak("k,l,a,j") * 2.0 - g_ijak("l,k,a,j"))

                           // (Gamma(cd,ij) + Gamma(ij, cd)) g^cd_aj
                         + (Gamma_AbIj_ab("c,d,i,j") + Gamma_IjAb_ab("i,j,c,d"))
                           * (g_abci("c,d,a,j") * 2.0 - g_abci("d,c,a,j"))

                           // (- Gamma(ck,ij) + Gamma(ij,kc)) g^kc_aj
                         + (- (Gamma_AkIj_ab("c,k,i,j") - Gamma_AkIj_ab("c,k,j,i"))
                            + Gamma_IjKa_ab("i,j,k,c") - Gamma_IjKa_ab("j,i,k,c"))
                           * (g_iabj("j,a,c,k") - g_aibj("a,j,c,k"))
                           // - Gamma_aKIj_ba("c,k,i,j") = Gamma_AkIj_ab("c,k,j,i")
                         + ( Gamma_AkIj_ab("c,k,j,i") + Gamma_IjKa_ab("i,j,k,c"))
                           * g_iabj("j,a,c,k")
                         - ( - Gamma_AkIj_ab("c,k,i,j") - Gamma_IjKa_ab("j,i,k,c"))
                           * g_aibj("a,j,c,k")

                           //
                           // (Gamma(kl,ib) - Gamma(bi,kl)) g^kl_ab
                         + (Gamma_IjKa_ab("k,l,i,b") + Gamma_AkIj_ab("b,i,l,k"))
                           * (g_ijab("k,l,a,b") * 2.0 - g_ijab("l,k,a,b"))

                           // (- Gamma(cd,bi) + Gamma(bi,dc)) g^cd_ab
                         + (Gamma_AbCi_ab("d,c,b,i") + Gamma_CiAb_ab("b,i,d,c"))
                           * (g_abcd("c,d,a,b") * 2.0 - g_abcd("d,c,a,b"))

                           // (Gamma(ck,ib) + Gamma(bi,kc)) g^ck_ab
                         + (  Gamma_aJiB_ab("c,k,i,b") + Gamma_AjiB_ba("c,k,i,b")
                            + Gamma_aJiB_ab("b,i,k,c") + Gamma_AjiB_ba("b,i,k,c"))
                           * (g_abci("a,b,c,k") - g_abci("b,a,c,k"))
                         + (Gamma_aJiB_ab("c,k,i,b") + Gamma_aJiB_ab("b,i,k,c"))
                           * g_abci("a,b,c,k")
                         - (Gamma_AjiB_ba("c,k,i,b") + Gamma_AjiB_ba("b,i,k,c"))
                           * g_abci("b,a,c,k")

                        )
                 ;
    Xam_tot = XaiAddToXam(Xam, Xai);

    // *** test code: print out Iam & Ima ***
//    TArray4d g_iabc = ijxy("<i a|g|b c>");
//    TArray4d g_ikbj = ijxy("<i k|g|b j>");
//    TArray4d g_ikjb = ijxy("<i k|g|j b>");
//    TArray4d g_kabi = ijxy("<k a|g|b i>");
//    TArray4d g_kaib = ijxy("<k a|g|i b>");
//    TArray4d g_ijkl = ijxy("<i j|g|k l>");
//
//    TArray2 Iam;
//    Iam("a,m") =   2.0 * (
//                           //
//                           // (Gamma(kl,ab) + Gamma(ab,kl)) g^kl_mb
//                          (Gamma_IjAb_ab("k,l,a,b") + Gamma_AbIj_ab("a,b,k,l"))
//                          * (g_ijma("k,l,m,b") * 2.0 - g_ijma("l,k,m,b"))
//
//                           // (Gamma(cd,ab) + Gamma(ab, cd)) g^cd_mb
//                         + (Gamma_AbCd_ab("c,d,a,b") + Gamma_AbCd_ab("a,b,c,d"))
//                           * (g_mabc("m,b,c,d") * 2.0 - g_mabc("m,b,d,c"))
//
//                           // 2 (Gamma(ck,ab) - Gamma(ba,ck)) g^ck_mb
//                         + (  Gamma_CiAb_ab("c,k,a,b") - Gamma_CiAb_ab("c,k,b,a")
//                            - Gamma_AbCi_ab("b,a,c,k") + Gamma_AbCi_ab("a,b,c,k"))
//                           * (g_mabi("m,b,c,k") - g_maib("m,b,k,c"))
//                         + (Gamma_CiAb_ab("c,k,a,b") + Gamma_AbCi_ab("a,b,c,k"))
//                           * g_mabi("m,b,c,k")
//                         - (- Gamma_CiAb_ab("c,k,b,a") - Gamma_AbCi_ab("b,a,c,k"))
//                           * g_maib("m,b,k,c")
//
//                           //
//                           // (Gamma(lk, ja) + Gamma(aj, kl)) g^kl_mj
//                         + (Gamma_IjKa_ab("l,k,j,a") + Gamma_AkIj_ab("a,j,k,l"))
//                           * (g_imkl("j,m,l,k") * 2.0 - g_imkl("j,m,k,l"))
//
//                           // (Gamma(cd, aj) + Gamma(aj, cd) g^mj_cd) g^cd_mj
//                         + (Gamma_AbCi_ab("c,d,a,j") + Gamma_CiAb_ab("a,j,c,d"))
//                           * (g_miab("m,j,c,d") * 2.0 - g_miab("m,j,d,c"))
//
//                           // * (Gamma(ck, ja) + Gamma(aj,kc)) g^kc_mj
//                         + (  Gamma_aJiB_ab("c,k,j,a") + Gamma_AjiB_ba("c,k,j,a")
//                            + Gamma_aJiB_ab("a,j,k,c") + Gamma_AjiB_ba("a,j,k,c"))
//                           * (g_imbj("j,m,c,k") - g_imjb("j,m,k,c"))
//                         + (Gamma_aJiB_ab("c,k,j,a") + Gamma_aJiB_ab("a,j,k,c"))
//                           * g_imbj("j,m,c,k")
//                         - (Gamma_AjiB_ba("c,k,j,a") + Gamma_AjiB_ba("a,j,k,c"))
//                           * g_imjb("j,m,k,c")
//                        )
//                ;
//
//    TArray2 Ima, Iia;
//    Ima("a,m") = - (2.0 * g_mabc("m,c,a,d") - g_mabc("m,c,d,a"))
//                   * (Dab("c,d") + Dab("d,c")) * 0.5
//                 - (2.0 * g_imjb("l,m,k,a") - g_imbj("l,m,a,k"))
//                   * (Dij("k,l") + Dij("l,k")) * 0.5
//                 ;
//    Iia("a,i") = - 2.0 * (
//                           //
//                           // (Gamma(kl,ij) + Gamma(ij, kl)) g^kl_aj
//                           (Gamma_IjKl_ab("k,l,i,j") + Gamma_IjKl_ab("i,j,k,l"))
//                           * (g_ijak("k,l,a,j") * 2.0 - g_ijak("l,k,a,j"))
//
//                           // (Gamma(cd,ij) + Gamma(ij, cd)) g^cd_aj
//                         + (Gamma_AbIj_ab("c,d,i,j") + Gamma_IjAb_ab("i,j,c,d"))
//                           * (g_abci("c,d,a,j") * 2.0 - g_abci("d,c,a,j"))
//
//                           // (- Gamma(ck,ij) + Gamma(ij,kc)) g^kc_aj
//                         + (- (Gamma_AkIj_ab("c,k,i,j") - Gamma_AkIj_ab("c,k,j,i"))
//                            + Gamma_IjKa_ab("i,j,k,c") - Gamma_IjKa_ab("j,i,k,c"))
//                           * (g_iabj("j,a,c,k") - g_aibj("a,j,c,k"))
//                           // - Gamma_aKIj_ba("c,k,i,j") = Gamma_AkIj_ab("c,k,j,i")
//                         + ( Gamma_AkIj_ab("c,k,j,i") + Gamma_IjKa_ab("i,j,k,c"))
//                           * g_iabj("j,a,c,k")
//                         - ( - Gamma_AkIj_ab("c,k,i,j") - Gamma_IjKa_ab("j,i,k,c"))
//                           * g_aibj("a,j,c,k")
//
//                           //
//                           // (Gamma(kl,ib) - Gamma(bi,kl)) g^kl_ab
//                         + (Gamma_IjKa_ab("k,l,i,b") + Gamma_AkIj_ab("b,i,l,k"))
//                           * (g_ijab("k,l,a,b") * 2.0 - g_ijab("l,k,a,b"))
//
//                           // (- Gamma(cd,bi) + Gamma(bi,dc)) g^cd_ab
//                         + (Gamma_AbCi_ab("d,c,b,i") + Gamma_CiAb_ab("b,i,d,c"))
//                           * (g_abcd("c,d,a,b") * 2.0 - g_abcd("d,c,a,b"))
//
//                           // (Gamma(ck,ib) + Gamma(bi,kc)) g^ck_ab
//                         + (  Gamma_aJiB_ab("c,k,i,b") + Gamma_AjiB_ba("c,k,i,b")
//                            + Gamma_aJiB_ab("b,i,k,c") + Gamma_AjiB_ba("b,i,k,c"))
//                           * (g_abci("a,b,c,k") - g_abci("b,a,c,k"))
//                         + (Gamma_aJiB_ab("c,k,i,b") + Gamma_aJiB_ab("b,i,k,c"))
//                           * g_abci("a,b,c,k")
//                         - (Gamma_AjiB_ba("c,k,i,b") + Gamma_AjiB_ba("b,i,k,c"))
//                           * g_abci("b,a,c,k")
//
//                        )
//                 ;
//    TArray2 I_ima = XaiAddToXam(Ima, Iia);
//
//    std::cout << "Ii/m a : " << std::endl << I_ima << std::endl;
//    std::cout << "Iam : " << std::endl << Iam << std::endl;

    // frozen-core contribution
    const std::size_t nocc = Xam.trange().elements().extent()[1];
    const std::size_t naocc = Xai.trange().elements().extent()[1];
    if (nocc != naocc) {
      TArray4d g_ijipk = ijxy("<i j|g|i' k>");
      TArray4d g_abipi = ijxy("<a b|g|i' i>");
      TArray4d g_ipibj = ijxy("<i' i|g|b j>");
      TArray4d g_iipbj = ijxy("<i i'|g|b j>");
      TArray4d g_ijipb = ijxy("<i j|g|i' b>");
      TArray4d g_abipd = ijxy("<a b|g|i' d>");
      TArray4d g_aipbi = ijxy("<a i'|g|b i>");

      Xiip("i,i'") = - 4.0 * (
                             //
                             // (Gamma(kl,ij) + Gamma(ij, kl)) g^kl_i'j
                             (Gamma_IjKl_ab("k,l,i,j") + Gamma_IjKl_ab("i,j,k,l"))
                             * (g_ijipk("k,l,i',j") * 2.0 - g_ijipk("l,k,i',j"))

                             // (Gamma(cd,ij) + Gamma(ij, cd)) g^cd_i'j
                           + (Gamma_AbIj_ab("c,d,i,j") + Gamma_IjAb_ab("i,j,c,d"))
                             * (g_abipi("c,d,i',j") * 2.0 - g_abipi("d,c,i',j"))

                             // (- Gamma(ck,ij) + Gamma(ij,kc)) g^kc_i'j
                           + (- (Gamma_AkIj_ab("c,k,i,j") - Gamma_AkIj_ab("c,k,j,i"))
                              + Gamma_IjKa_ab("i,j,k,c") - Gamma_IjKa_ab("j,i,k,c"))
                             * (g_iipbj("j,i',c,k") - g_ipibj("i',j,c,k"))
                           + ( Gamma_AkIj_ab("c,k,j,i") + Gamma_IjKa_ab("i,j,k,c"))
                             * g_iipbj("j,i',c,k")
                           - ( - Gamma_AkIj_ab("c,k,i,j") - Gamma_IjKa_ab("j,i,k,c"))
                             * g_ipibj("i',j,c,k")

                             //
                             // (Gamma(kl,ib) - Gamma(bi,kl)) g^kl_i'b
                           + (Gamma_IjKa_ab("k,l,i,b") + Gamma_AkIj_ab("b,i,l,k"))
                             * (g_ijipb("k,l,i',b") * 2.0 - g_ijipb("l,k,i',b"))

                             // (- Gamma(cd,bi) + Gamma(bi,dc)) g^cd_i'b
                           + (Gamma_AbCi_ab("d,c,b,i") + Gamma_CiAb_ab("b,i,d,c"))
                             * (g_abipd("c,d,i',b") * 2.0 - g_abipd("d,c,i',b"))

                             // (Gamma(ck,ib) + Gamma(bi,kc)) g^ck_i'b
                           + (  Gamma_aJiB_ab("c,k,i,b") + Gamma_AjiB_ba("c,k,i,b")
                              + Gamma_aJiB_ab("b,i,k,c") + Gamma_AjiB_ba("b,i,k,c"))
                             * (g_abipi("c,b,i',k") - g_aipbi("b,i',c,k"))
                           + (Gamma_aJiB_ab("c,k,i,b") + Gamma_aJiB_ab("b,i,k,c"))
                             * g_abipi("c,b,i',k")
                           - (Gamma_AjiB_ba("c,k,i,b") + Gamma_AjiB_ba("b,i,k,c"))
                             * g_aipbi("b,i',c,k")
                          )
                   ;

    }

  }

  template <typename T>
  void SingleReference_R12Intermediates<T>::compute_multipole() {

    ExEnv::out0() << std::endl << indent
                  << "Compute dipole and quadrupole moments" << std::endl;

    bool compute_dipole = true;
    bool compute_quadrupole = true;
    bool compute_EFG = false;

    #define INCLUDE_CABS_Singles_CONTRI 1
    #define INCLUDE_MP2_CONTRI 1
    #define INCLUDE_MP2F12C_CONTRI 1
    #define INCLUDE_F12_CONTRI 1

    #define INCLUDE_CCSD_CONTRI 1
    #define INCLUDE_CCSDF12C_CONTRI 1

    double conv_target = 1e-10;

    // Obtain property integrals which is needed in multiple places

    // dipole integrals
    TArray2 mu_z_mn, mu_z_am;
    if (compute_dipole) {
      mu_z_mn = xy("<m|mu_z|n>");
      mu_z_am = xy("<a|mu_z|m>");
    }

    // quadrupole integrals
    TArray2 q_xx_mn, q_yy_mn, q_zz_mn;
    TArray2 Qxx_mn, Qyy_mn, Qzz_mn, Qxz_mn, Qxy_mn, Qyz_mn;
    TArray2 q_xx_am, q_yy_am, q_zz_am;
    TArray2 Qxx_am, Qyy_am, Qzz_am, Qxy_am, Qxz_am, Qyz_am;
    if (compute_quadrupole) {
      q_xx_mn = xy("<m|q_xx|n>");
      q_yy_mn = xy("<m|q_yy|n>");
      q_zz_mn = xy("<m|q_zz|n>");
      Qxx_mn("m,n") = q_xx_mn("m,n") - (q_zz_mn("m,n") + q_yy_mn("m,n")) * 0.5;
      Qyy_mn("m,n") = q_yy_mn("m,n") - (q_zz_mn("m,n") + q_xx_mn("m,n")) * 0.5;
      Qzz_mn("m,n") = q_zz_mn("m,n") - (q_xx_mn("m,n") + q_yy_mn("m,n")) * 0.5;
      Qxz_mn("m,n") = _2("<m|q_xz|n>") * 1.5;
      Qxy_mn("m,n") = _2("<m|q_xy|n>") * 1.5;
      Qyz_mn("m,n") = _2("<m|q_yz|n>") * 1.5;

      q_xx_am = xy("<a|q_xx|m>");
      q_yy_am = xy("<a|q_yy|m>");
      q_zz_am = xy("<a|q_zz|m>");
      Qxx_am("a,m") = q_xx_am("a,m") - (q_zz_am("a,m") + q_yy_am("a,m")) * 0.5;
      Qyy_am("a,m") = q_yy_am("a,m") - (q_zz_am("a,m") + q_xx_am("a,m")) * 0.5;
      Qzz_am("a,m") = q_zz_am("a,m") - (q_xx_am("a,m") + q_yy_am("a,m")) * 0.5;
      Qxy_am("a,m") = _2("<a|q_xy|m>") * 1.5;
      Qxz_am("a,m") = _2("<a|q_xz|m>") * 1.5;
      Qyz_am("a,m") = _2("<a|q_yz|m>") * 1.5;
    }

    // electric field gradient integrals
    TArray2 v_xx_mn, v_yy_mn, v_zz_mn, v_xy_mn, v_xz_mn, v_yz_mn;
    TArray2 v_xx_am, v_yy_am, v_zz_am, v_xy_am, v_xz_am, v_yz_am;
    if (compute_EFG) {
      v_xx_mn = xy("<m|ddphi_xx|n>");
      v_yy_mn = xy("<m|ddphi_yy|n>");
      v_zz_mn = xy("<m|ddphi_zz|n>");
      v_xy_mn = xy("<m|ddphi_xy|n>");
      v_xz_mn = xy("<m|ddphi_xz|n>");
      v_yz_mn = xy("<m|ddphi_yz|n>");

      v_xx_am = xy("<a|ddphi_xx|m>");
      v_yy_am = xy("<a|ddphi_yy|m>");
      v_zz_am = xy("<a|ddphi_zz|m>");
      v_xy_am = xy("<a|ddphi_xy|m>");
      v_xz_am = xy("<a|ddphi_xz|m>");
      v_yz_am = xy("<a|ddphi_yz|m>");
    }

    // compute HF contribution
    Ref<Molecule> mol = r12world_->basis()->molecule();
    const int natom = mol->natom();

    TArray2 Imn = xy("<m|I|n>");

    // HF eletric dipole
    double mu_z_n = 0.0, mu_z_scf = 0.0;
    if (compute_dipole) {
      // Nuclear contribution to dipole
      for(int a = 0; a < natom; ++a) {
        const double x = mol->r(a, 0);
        const double y = mol->r(a, 1);
        const double z = mol->r(a, 2);
        const double Z = mol->Z(a);
        mu_z_n += Z * z;
      }

      // SCF contribution to electronic electric dipole
      mu_z_scf = dot(mu_z_mn("m,n"), Imn("m,n"));
    }

    // HF eletric quadrupole
    double q_xx_n = 0.0, q_yy_n = 0.0, q_zz_n = 0.0;
    double q_xy_n = 0.0, q_xz_n = 0.0, q_yz_n = 0.0;
    double q_xx_scf = 0.0, q_yy_scf = 0.0, q_zz_scf = 0.0,
           q_xy_scf = 0.0, q_xz_scf = 0.0, q_yz_scf = 0.0;
    if (compute_quadrupole) {
      // Nuclear electric quadrupole (traceless)
      for(int a = 0; a < natom; ++a) {
        const double x = mol->r(a, 0);
        const double y = mol->r(a, 1);
        const double z = mol->r(a, 2);
        const double Z = mol->Z(a);

        q_xx_n += Z * (x * x - (z * z + y * y) * 0.5); // traceless form of the quadrupole
        q_yy_n += Z * (y * y - (z * z + x * x) * 0.5); // traceless form of the quadrupole
        q_zz_n += Z * (z * z - (x * x + y * y) * 0.5); // traceless form of the quadrupole
        q_xy_n += Z * (x * y * 1.5); // traceless form of the quadrupole
        q_xz_n += Z * (x * z * 1.5); // traceless form of the quadrupole
        q_yz_n += Z * (y * z * 1.5); // traceless form of the quadrupole
      }

      // SCF contribution to electronic electric quadrupoles
      q_xx_scf = dot(Qxx_mn("m,n"), Imn("m,n"));
      q_yy_scf = dot(Qyy_mn("m,n"), Imn("m,n"));
      q_zz_scf = dot(Qzz_mn("m,n"), Imn("m,n"));
      q_xy_scf = dot(Qxy_mn("m,n"), Imn("m,n"));
      q_xz_scf = dot(Qxz_mn("m,n"), Imn("m,n"));
      q_yz_scf = dot(Qyz_mn("m,n"), Imn("m,n"));
    }

    // HF electric field gradient near the 1st nucleus
    if (compute_EFG) {
      // Nuclear contribution to electric field gradient near the 1st nucleus
      double v_xx_n = 0.0, v_yy_n = 0.0, v_zz_n = 0.0;
      double v_xy_n = 0.0, v_xz_n = 0.0, v_yz_n = 0.0;
      for(int i = 1; i < natom; i++) {
        const double x = mol->r(0, 0) - mol->r(i, 0);
        const double y = mol->r(0, 1) - mol->r(i, 1);
        const double z = mol->r(0, 2) - mol->r(i, 2);
        const double r2 = x * x + y * y + z * z;
        const double r = pow(r2, 0.5);
        const double Z = mol->Z(i);

        v_xx_n += - Z * (3.0 * x * x - r2) / (r * r2 * r2);
        v_yy_n += - Z * (3.0 * y * y - r2) / (r * r2 * r2);
        v_zz_n += - Z * (3.0 * z * z - r2) / (r * r2 * r2);

        v_xy_n += - Z * (3.0 * x * y) / (r * r2 * r2);
        v_xz_n += - Z * (3.0 * x * z) / (r * r2 * r2);
        v_yz_n += - Z * (3.0 * y * z) / (r * r2 * r2);
      }

      // Electric field gradient near the 1st nucleus
      const double v_xx_scf = dot(v_xx_mn("m,n"), Imn("m,n"));
      const double v_yy_scf = dot(v_yy_mn("m,n"), Imn("m,n"));
      const double v_zz_scf = dot(v_zz_mn("m,n"), Imn("m,n"));
      const double v_xy_scf = dot(v_xy_mn("m,n"), Imn("m,n"));
      const double v_xz_scf = dot(v_xz_mn("m,n"), Imn("m,n"));
      const double v_yz_scf = dot(v_yz_mn("m,n"), Imn("m,n"));

      std::cout << std::endl << indent
                << "electric gradient (HF=SCF+N)" << std::endl
                << indent << "v_xx (HF) = " << scprintf("%12.10f", v_xx_n + v_xx_scf * 2.0)
                << indent << "  v_yy (HF) = " << scprintf("%12.10f", v_yy_n + v_yy_scf * 2.0)
                << indent << "  v_zz (HF) = " << scprintf("%12.10f", v_zz_n + v_zz_scf * 2.0)
                << std::endl
                << indent << "v_xy (HF) = " << scprintf("%12.10f", v_xy_n + v_xy_scf * 2.0)
                << indent << "  v_xz (HF) = " << scprintf("%12.10f", v_xz_n + v_xz_scf * 2.0)
                << indent << "  v_yz (HF) = " << scprintf("%12.10f", v_yz_n + v_yz_scf * 2.0)
                << std::endl;
    }

    // compute integrals needed for orbital relaxation
    // i.e. solve Abnam Dbn = Xam (close-shell formula)
    TArray4d g_mnab = ijxy("<m n|g|a b>");
    TArray4 A_bnam;
    A_bnam("b,n,a,m") = - _4("<b n|g|a m>") - g_mnab("m,n,b,a") + 4.0 * g_mnab("n,m,b,a")
                        + _2("<b|F|a>") * Imn("m,n") - _2("<a|I|b>") * _2("<m|F|n>");

    // Make preconditioner: Delta_am = 1 / (<a|F|a> - <m|F|m>) for
    // solving k_bn A_bnam = X_am
    TArray2 mFmn, mFab;
    mFmn("m,n") = - _2("<m|F|n>");
    mFab("a,b") = - _2("<a|F|b>");
    typedef detail::diag_precond2<double> pceval_type; //!< evaluator of preconditioner
    pceval_type Delta_am_gen(TA::array_to_eigen(mFab), TA::array_to_eigen(mFmn));

    typedef TA::Array<T, 2, LazyTensor<T, 2, pceval_type > > TArray2d;
    TArray2 Iam = xy("<a|I|m>");
    TArray2d Delta_am(Iam.get_world(), Iam.trange());

    // construct local tiles
    for(auto t = Delta_am.trange().tiles().begin();
        t != Delta_am.trange().tiles().end(); ++t)
      if (Delta_am.is_local(*t)) {
        std::array<std::size_t, 2> index;
        std::copy(t->begin(), t->end(), index.begin());
        madness::Future < typename TArray2d::value_type >
          tile((LazyTensor<T, 2, pceval_type >(&Delta_am, index, &Delta_am_gen)
              ));

        // Insert the tile into the array
        Delta_am.set(*t, tile);
      }
    TArray2 preconditioner;
    preconditioner("a,m") = Delta_am("a,m");

    detail::Orbital_relaxation_Abjai<double> Orbital_relaxation_Abnam(A_bnam);
    TA::ConjugateGradientSolver<TiledArray::Array<T,2>,
                                detail::Orbital_relaxation_Abjai<double> > cg_solver2;

    // compute CABS singles contribution
    double mu_z_e2 = 0.0, mu_z_e2or = 0.0;
    double q_xx_e2 = 0.0, q_yy_e2 = 0.0, q_zz_e2 = 0.0;
    double q_xx_e2or = 0.0, q_yy_e2or = 0.0, q_zz_e2or = 0.0;
#if INCLUDE_CABS_Singles_CONTRI
    {
    ExEnv::out0() << std::endl << indent
                  << "Compute CABS_singles contributions" << std::endl;
    TArray2 TmA = xy("<m|T1|A'>");
    TArray2 Tma = xy("<m|T1|a>");

    // density from CABS Singles contribution
    // D^m_n =  t^m_A' * t^A'_n
    // D^A'_B' = t^A'_m * t^m_B'
    // D^A'_m = t^A'_m
    ExEnv::out0() << indent
                  << "Compute CABS_singles D" << std::endl;
    TArray2 D_e2_mn, D_e2_AB;
    D_e2_mn("m,n") = TmA("m,A'") * TmA("n,A'");
    D_e2_AB("A',B'") = TmA("m,A'") * TmA("m,B'");

    // CABS singles orbital response
    ExEnv::out0() << indent
                  << "Compute CABS_singles Xam" << std::endl;
    TArray2 Xam_E2 = Xam_CabsSingles(TmA, Tma);
    TArray2 Dbn_E2(Xam_E2.get_world(), Xam_E2.trange());
    // solve k_bn A_bnam = X_am
    auto resnorm_E2 = cg_solver2(Orbital_relaxation_Abnam,
                                 Xam_E2,
                                 Dbn_E2,
                                 preconditioner,
                                 conv_target);

    // compute CABS singles contribution to dipole
    if (compute_dipole) {
      TArray2 mu_z_AB = xy("<A'|mu_z|B'>");
      TArray2 mu_z_mA = xy("<m|mu_z|A'>");

      mu_z_e2 = - dot(mu_z_mn("m,n"), D_e2_mn("m,n"))
                + dot(mu_z_AB("A',B'"), D_e2_AB("A',B'"))
                + dot(mu_z_mA("m,A'"), TmA("m,A'")) * 2.0;

      // CABS singles orbital response contribution to dipole
      mu_z_e2or = dot(mu_z_am("a,m"), Dbn_E2("a,m"));
    }

    // CABS singles contribution to quadrupoles
    if (compute_quadrupole) {
      TArray2 q_xx_AB = xy("<A'|q_xx|B'>");
      TArray2 q_yy_AB = xy("<A'|q_yy|B'>");
      TArray2 q_zz_AB = xy("<A'|q_zz|B'>");
      TArray2 Qxx_AB, Qyy_AB, Qzz_AB;
      Qxx_AB("A',B'") = q_xx_AB("A',B'") - (q_zz_AB("A',B'") + q_yy_AB("A',B'")) * 0.5;
      Qyy_AB("A',B'") = q_yy_AB("A',B'") - (q_xx_AB("A',B'") + q_zz_AB("A',B'")) * 0.5;
      Qzz_AB("A',B'") = q_zz_AB("A',B'") - (q_xx_AB("A',B'") + q_yy_AB("A',B'")) * 0.5;

      TArray2 q_xx_mA = xy("<m|q_xx|A'>");
      TArray2 q_yy_mA = xy("<m|q_yy|A'>");
      TArray2 q_zz_mA = xy("<m|q_zz|A'>");
      TArray2 Qxx_mA, Qyy_mA, Qzz_mA;
      Qxx_mA("m,A'") = q_xx_mA("m,A'") - (q_zz_mA("m,A'")  + q_yy_mA("m,A'") ) * 0.5;
      Qyy_mA("m,A'") = q_yy_mA("m,A'") - (q_xx_mA("m,A'")  + q_zz_mA("m,A'") ) * 0.5;
      Qzz_mA("m,A'") = q_zz_mA("m,A'") - (q_xx_mA("m,A'")  + q_yy_mA("m,A'") ) * 0.5;

      q_xx_e2 = - dot(Qxx_mn("m,n"), D_e2_mn("m,n"))
                + dot(Qxx_AB("A',B'"), D_e2_AB("A',B'"))
                + dot(Qxx_mA("m,A'"), TmA("m,A'")) * 2.0;
      q_yy_e2 = - dot(Qyy_mn("m,n"), D_e2_mn("m,n"))
                + dot(Qyy_AB("A',B'"), D_e2_AB("A',B'"))
                + dot(Qyy_mA("m,A'"), TmA("m,A'")) * 2.0;
      q_zz_e2 = - dot(Qzz_mn("m,n"), D_e2_mn("m,n"))
                + dot(Qzz_AB("A',B'"), D_e2_AB("A',B'"))
                + dot(Qzz_mA("m,A'"), TmA("m,A'")) * 2.0;

      q_xx_e2or = dot(Qxx_am("a,m"), Dbn_E2("a,m"));
      q_yy_e2or = dot(Qyy_am("a,m"), Dbn_E2("a,m"));
      q_zz_e2or = dot(Qzz_am("a,m"), Dbn_E2("a,m"));
    }

    // CABS singles contribution to
    // electric field gradient near the 1st nucleus
    if (compute_EFG) {
      TArray2 v_xx_mA = xy("<m|ddphi_xx|A'>");
      TArray2 v_yy_mA = xy("<m|ddphi_yy|A'>");
      TArray2 v_zz_mA = xy("<m|ddphi_zz|A'>");
      TArray2 v_xy_mA = xy("<m|ddphi_xy|A'>");
      TArray2 v_xz_mA = xy("<m|ddphi_xz|A'>");
      TArray2 v_yz_mA = xy("<m|ddphi_yz|A'>");

      TArray2 v_xx_AB = xy("<A'|ddphi_xx|B'>");
      TArray2 v_yy_AB = xy("<A'|ddphi_yy|B'>");
      TArray2 v_zz_AB = xy("<A'|ddphi_zz|B'>");
      TArray2 v_xy_AB = xy("<A'|ddphi_xy|B'>");
      TArray2 v_xz_AB = xy("<A'|ddphi_xz|B'>");
      TArray2 v_yz_AB = xy("<A'|ddphi_yz|B'>");

      const double v_xx_e2 = - dot(v_xx_mn("m,n"), D_e2_mn("m,n"))
                             + dot(v_xx_AB("A',B'"), D_e2_AB("A',B'"))
                             + dot(v_xx_mA("m,A'"), TmA("m,A'")) * 2.0;
      const double v_yy_e2 = - dot(v_yy_mn("m,n"), D_e2_mn("m,n"))
                             + dot(v_yy_AB("A',B'"), D_e2_AB("A',B'"))
                             + dot(v_yy_mA("m,A'"), TmA("m,A'")) * 2.0;
      const double v_zz_e2 = - dot(v_zz_mn("m,n"), D_e2_mn("m,n"))
                             + dot(v_zz_AB("A',B'"), D_e2_AB("A',B'"))
                             + dot(v_zz_mA("m,A'"), TmA("m,A'")) * 2.0;
      const double v_xy_e2 = - dot(v_xy_mn("m,n"), D_e2_mn("m,n"))
                             + dot(v_xy_AB("A',B'"), D_e2_AB("A',B'"))
                             + dot(v_xy_mA("m,A'"), TmA("m,A'")) * 2.0;
      const double v_xz_e2 = - dot(v_xz_mn("m,n"), D_e2_mn("m,n"))
                             + dot(v_xz_AB("A',B'"), D_e2_AB("A',B'"))
                             + dot(v_xz_mA("m,A'"), TmA("m,A'")) * 2.0;
      const double v_yz_e2 = - dot(v_yz_mn("m,n"), D_e2_mn("m,n"))
                             + dot(v_yz_AB("A',B'"), D_e2_AB("A',B'"))
                             + dot(v_yz_mA("m,A'"), TmA("m,A'")) * 2.0;

      const double v_xx_e2or = dot(v_xx_am("a,m"), Dbn_E2("a,m"));
      const double v_yy_e2or = dot(v_yy_am("a,m"), Dbn_E2("a,m"));
      const double v_zz_e2or = dot(v_zz_am("a,m"), Dbn_E2("a,m"));
      const double v_xy_e2or = dot(v_xy_am("a,m"), Dbn_E2("a,m"));
      const double v_xz_e2or = dot(v_xz_am("a,m"), Dbn_E2("a,m"));
      const double v_yz_e2or = dot(v_yz_am("a,m"), Dbn_E2("a,m"));

      std::cout << std::endl << indent
                << "electric gradient (E2)" << std::endl
                << indent << "v_xx (E2) = " << scprintf("%12.10f", v_xx_e2 * 2.0)
                << indent << "  v_yy (E2) = " << scprintf("%12.10f", v_yy_e2 * 2.0)
                << indent << "  v_zz (E2) = " << scprintf("%12.10f", v_zz_e2 * 2.0)
                << std::endl
                << indent << "v_xy (E2) = " << scprintf("%12.10f", v_xy_e2 * 2.0)
                << indent << "  v_xz (E2) = " << scprintf("%12.10f", v_xz_e2 * 2.0)
                << indent << "  v_yz (E2) = " << scprintf("%12.10f", v_yz_e2 * 2.0)
                << std::endl;
      std::cout << std::endl << indent
                << "electric gradient (E2 orbital response)" << std::endl
                << indent  << "v_xx (E2 or) = " << scprintf("%12.10f", v_xx_e2or * 2.0)
                << indent << "  v_yy (E2 or) = " << scprintf("%12.10f", v_yy_e2or * 2.0)
                << indent << "  v_zz (E2 or) = " << scprintf("%12.10f", v_zz_e2or * 2.0)
                << std::endl
                << indent << "v_xy (E2 or) = " << scprintf("%12.10f", v_xy_e2or * 2.0)
                << indent << "  v_xz (E2 or) = " << scprintf("%12.10f", v_xz_e2or * 2.0)
                << indent << "  v_yz (E2 or) = " << scprintf("%12.10f", v_yz_e2or * 2.0)
                << std::endl << std::endl;
    }

    }
    world_.gop.fence(); // clear memory
#endif

    // determine if it is frozen-core
    TArray2 Fij = xy("<i|F|j>");
    const std::size_t nocc = mFmn.trange().elements().extent()[0];
    const std::size_t naocc = Fij.trange().elements().extent()[0];

    // compute integrals needed for MP2 and following calculations of properties
    // dipole integrals
    TArray2 mu_z_ij, mu_z_ab, mu_z_ijp;
    if (compute_dipole) {
      mu_z_ij = xy("<i|mu_z|j>");
      mu_z_ab = xy("<a|mu_z|b>");

      if (nocc != naocc) {
        mu_z_ijp = xy("<i|mu_z|j'>");
      }
    }

    // quadrupole integrals
    TArray2 q_xx_ij, q_yy_ij, q_zz_ij, Qxx_ij, Qyy_ij, Qzz_ij;
    TArray2 q_xx_ab, q_yy_ab, q_zz_ab, Qxx_ab, Qyy_ab, Qzz_ab;
    // integrals for frozen-core
    TArray2 Qxx_ijp, Qyy_ijp, Qzz_ijp;
    if (compute_quadrupole){
      q_xx_ij = xy("<i|q_xx|j>");
      q_yy_ij = xy("<i|q_yy|j>");
      q_zz_ij = xy("<i|q_zz|j>");
      Qxx_ij("i,j") = q_xx_ij("i,j") - (q_zz_ij("i,j") + q_yy_ij("i,j")) * 0.5;
      Qyy_ij("i,j") = q_yy_ij("i,j") - (q_zz_ij("i,j") + q_xx_ij("i,j")) * 0.5;
      Qzz_ij("i,j") = q_zz_ij("i,j") - (q_xx_ij("i,j") + q_yy_ij("i,j")) * 0.5;

      q_xx_ab = xy("<a|q_xx|b>");
      q_yy_ab = xy("<a|q_yy|b>");
      q_zz_ab = xy("<a|q_zz|b>");
      Qxx_ab("a,b") = q_xx_ab("a,b") - (q_zz_ab("a,b") + q_yy_ab("a,b")) * 0.5;
      Qyy_ab("a,b") = q_yy_ab("a,b") - (q_zz_ab("a,b") + q_xx_ab("a,b")) * 0.5;
      Qzz_ab("a,b") = q_zz_ab("a,b") - (q_xx_ab("a,b") + q_yy_ab("a,b")) * 0.5;

      if (nocc != naocc) {
        TArray2 q_xx_ijp = xy("<i|q_xx|j'>");
        TArray2 q_yy_ijp = xy("<i|q_yy|j'>");
        TArray2 q_zz_ijp = xy("<i|q_zz|j'>");
        Qxx_ijp("i,j'") = q_xx_ijp("i,j'") - (q_zz_ijp("i,j'") + q_yy_ijp("i,j'")) * 0.5;
        Qyy_ijp("i,j'") = q_yy_ijp("i,j'") - (q_zz_ijp("i,j'") + q_xx_ijp("i,j'")) * 0.5;
        Qzz_ijp("i,j'") = q_zz_ijp("i,j'") - (q_xx_ijp("i,j'") + q_yy_ijp("i,j'")) * 0.5;
      }
    }

    // electric field gradient integrals
    TArray2 v_xx_ij, v_yy_ij, v_zz_ij, v_xy_ij, v_xz_ij, v_yz_ij;
    TArray2 v_xx_ab, v_yy_ab, v_zz_ab, v_xy_ab, v_xz_ab, v_yz_ab;
    if (compute_EFG) {
      v_xx_ij = xy("<i|ddphi_xx|j>");
      v_yy_ij = xy("<i|ddphi_yy|j>");
      v_zz_ij = xy("<i|ddphi_zz|j>");
      v_xy_ij = xy("<i|ddphi_xy|j>");
      v_xz_ij = xy("<i|ddphi_xz|j>");
      v_yz_ij = xy("<i|ddphi_yz|j>");

      v_xx_ab = xy("<a|ddphi_xx|b>");
      v_yy_ab = xy("<a|ddphi_yy|b>");
      v_zz_ab = xy("<a|ddphi_zz|b>");
      v_xy_ab = xy("<a|ddphi_xy|b>");
      v_xz_ab = xy("<a|ddphi_xz|b>");
      v_yz_ab = xy("<a|ddphi_yz|b>");
    }

    // compute terms needed for frozen-core
    TArray2 Delta_ijp_F;
    TArray4 A_ijpam;
    if (nocc != naocc) {
      // compute 1 / (Fi'i' - Fjj)
      TArray2 Fipjp = xy("<i'|F|j'>");
      pceval_type Delta_ijp_gen(TA::array_to_eigen(Fij), TA::array_to_eigen(Fipjp));

      TArray2 Fijp = xy("<i|F|j'>");
      TArray2d Delta_ijp(Fijp.get_world(), Fijp.trange());
      // construct local tiles
      for(auto t = Delta_ijp.trange().tiles().begin();
          t != Delta_ijp.trange().tiles().end(); ++t)
        if (Delta_ijp.is_local(*t)) {
          std::array<std::size_t, 2> index;
          std::copy(t->begin(), t->end(), index.begin());
          madness::Future < typename TArray2d::value_type >
            tile((LazyTensor<T, 2, pceval_type >(&Delta_ijp, index, &Delta_ijp_gen)
                ));

          // Insert the tile into the array
          Delta_ijp.set(*t, tile);
        }
      Delta_ijp_F("i,j'") = Delta_ijp("i,j'");

      A_ijpam("i,j',a,m") = 4.0 * _4("<i a|g|j' m>") - _4("<i a|g|m j'>") - _4("<i m|g|a j'>");
    }

    // compute MP2 and its orbital response contribution
    // compute Delta_ijab = - 1 / (- <i|F|i> - <j|F|j> + <a|F|a> + <b|F|b>)
    // which is needed for MP2 amplitudes & MP2 F12 coupling
    TArray2 mFij;
    mFij("i,j") = - _2("<i|F|j>");
    TArray4d g_ijab = ijxy("<i j|g|a b>");

    typedef detail::diag_precond4<double> pc4eval_type;
    typedef TA::Array<T, 4, LazyTensor<T, 4, pc4eval_type > > TArray4dLazy;
    TArray4dLazy Delta_ijab(g_ijab.get_world(), g_ijab.trange());
    pc4eval_type Delta_ijab_gen(TA::array_to_eigen(mFij), TA::array_to_eigen(mFij),
                                TA::array_to_eigen(mFab),TA::array_to_eigen(mFab));
    // construct local tiles
    for(auto t = Delta_ijab.trange().tiles().begin();
        t != Delta_ijab.trange().tiles().end(); ++t)
      if (Delta_ijab.is_local(*t)) {
        std::array<std::size_t, 4> index;
        std::copy(t->begin(), t->end(), index.begin());
        madness::Future < typename TArray4dLazy::value_type >
          tile((LazyTensor<T, 4, pc4eval_type >(&Delta_ijab, index, &Delta_ijab_gen)
              ));

        // Insert the tile into the array
        Delta_ijab.set(*t, tile);
      }

    // MP2 amplitues which are needed for MP2 and coupling contibutions
    TArray4 T2_ijab;
    T2_ijab("i,j,a,b") = g_ijab("i,j,a,b") * Delta_ijab("i,j,a,b");

    double mu_z_mp2 = 0.0, mu_z_mp2or = 0.0;
    double q_xx_mp2 = 0.0, q_yy_mp2 = 0.0, q_zz_mp2 = 0.0;
    double q_xx_mp2or = 0.0, q_yy_mp2or = 0.0, q_zz_mp2or = 0.0;
#if INCLUDE_MP2_CONTRI
    {
    ExEnv::out0() << std::endl << indent
                  << "Compute MP2 contributions" << std::endl;
    // MP2 density
    TArray2 D_mp2_ij, D_mp2_ab;
    D_mp2_ij("i,j") = (2.0 * T2_ijab("i,k,a,b") - T2_ijab("k,i,a,b"))
                      * T2_ijab("j,k,a,b");
    D_mp2_ab("a,b") = (2.0 * T2_ijab("i,j,a,c") - T2_ijab("i,j,c,a"))
                      * T2_ijab("i,j,b,c");

    // MP2 orbital response
    TArray2 X_mp2_nfzc = Xam_mp2(T2_ijab, D_mp2_ij, D_mp2_ab);

    TArray2 X_mp2;
    // frozen-core contribution
    TArray2 Dijp_or_mp2;
    if (nocc != naocc) {
      ExEnv::out0() << std::endl << indent
                    << "Include frozen-core MP2 contributions" << std::endl;
      Dijp_or_mp2("i,j'") = - ( _4("<j' l|g|c d>")
                                * (2.0 * T2_ijab("i,l,c,d") - T2_ijab("l,i,c,d"))
                              ) * Delta_ijp_F("i,j'") * 2.0;
      TArray2 Xam_Dijp_or;
      Xam_Dijp_or("a,m") = A_ijpam("i,j',a,m") * Dijp_or_mp2("i,j'");
      X_mp2("a,m") = X_mp2_nfzc("a,m") - Xam_Dijp_or("a,m");

      // test code: print out Ima in frozen-core case
//      TArray2 Iia, Ima;
//      Iia("a,i") = _4("<a l|g|c d>")
//                   * (2.0 * T2_ijab("i,l,c,d") - T2_ijab("l,i,c,d"));
//      Ima("a,m") =   (2.0 * _4("<a c|g|m d>") - _4("<a c|g|d m>"))
//                     * D_mp2_ab("d,c")
//                   - (2.0 * _4("<a l|g|m k>") - _4("<a l|g|k m>"))
//                     * D_mp2_ij("k,l")
//                   ;
//      TArray2 Iima = XaiAddToXam(Ima, Iia);
//      TArray2 Iima_fzc;
//      Iima_fzc("a,m") = Iima("a,m") + Xam_Dijp_or("a,m") * 0.5;
//      ExEnv::out0() << std::endl << "FZC Ima : " << std::endl << Iima_fzc << std::endl;
    } else {
      X_mp2("a,m") = X_mp2_nfzc("a,m");
    }
     // *** end of frozen-core ***

    TArray2 Dbn_mp2(X_mp2.get_world(), X_mp2.trange());
    auto resnorm_mp2 = cg_solver2(Orbital_relaxation_Abnam,
                                  X_mp2,
                                  Dbn_mp2,
                                  preconditioner,
                                  conv_target);

    // MP2 dipole
    if (compute_dipole) {
      // MP2 density contribution to dipole
      mu_z_mp2 = - dot(mu_z_ij("i,j"), D_mp2_ij("i,j"))
                 + dot(mu_z_ab("a,b"), D_mp2_ab("a,b"));

      // MP2 orbital response contribution to dipole
      mu_z_mp2or = dot(mu_z_am("a,m"), Dbn_mp2("a,m"));

      if (nocc != naocc) {
        mu_z_mp2or += dot(mu_z_ijp("i,j'"), Dijp_or_mp2("i,j'"));
      }
    }

    // MP2 contribution to quadrupoles
    if (compute_quadrupole) {
      q_xx_mp2 = - dot(Qxx_ij("i,j"), D_mp2_ij("i,j"))
                 + dot(Qxx_ab("a,b"), D_mp2_ab("a,b"));
      q_yy_mp2 = - dot(Qyy_ij("i,j"), D_mp2_ij("i,j"))
                 + dot(Qyy_ab("a,b"), D_mp2_ab("a,b"));
      q_zz_mp2 = - dot(Qzz_ij("i,j"), D_mp2_ij("i,j"))
                 + dot(Qzz_ab("a,b"), D_mp2_ab("a,b"));
//      // x2, y2, z2 components to quadrupole moment
//      q_xx_mp2 = - dot(q_xx_ij("i,j"), D_mp2_ij("i,j"))
//                 + dot(q_xx_ab("a,b"), D_mp2_ab("a,b"));
//      q_yy_mp2 = - dot(q_yy_ij("i,j"), D_mp2_ij("i,j"))
//                 + dot(q_yy_ab("a,b"), D_mp2_ab("a,b"));
//      q_zz_mp2 = - dot(q_zz_ij("i,j"), D_mp2_ij("i,j"))
//                 + dot(q_zz_ab("a,b"), D_mp2_ab("a,b"));

      q_xx_mp2or = dot(Qxx_am("a,m"), Dbn_mp2("a,m"));
      q_yy_mp2or = dot(Qyy_am("a,m"), Dbn_mp2("a,m"));
      q_zz_mp2or = dot(Qzz_am("a,m"), Dbn_mp2("a,m"));
//      // x2, y2, z2 components to quadrupole moment
//      q_xx_mp2or = dot(q_xx_am("a,m"), Dbn_mp2("a,m"));
//      q_yy_mp2or = dot(q_yy_am("a,m"), Dbn_mp2("a,m"));
//      q_zz_mp2or = dot(q_zz_am("a,m"), Dbn_mp2("a,m"));

      if (nocc != naocc) {
          q_xx_mp2or += dot(Qxx_ijp("i,j'"), Dijp_or_mp2("i,j'"));
          q_yy_mp2or += dot(Qyy_ijp("i,j'"), Dijp_or_mp2("i,j'"));
          q_zz_mp2or += dot(Qzz_ijp("i,j'"), Dijp_or_mp2("i,j'"));
      }
    }

    // Electric field gradient near the 1st nucleus
    if (compute_EFG) {
      const double v_xx_mp2 = - dot(v_xx_ij("i,j"), D_mp2_ij("i,j"))
                              + dot(v_xx_ab("a,b"), D_mp2_ab("a,b"));
      const double v_yy_mp2 = - dot(v_yy_ij("i,j"), D_mp2_ij("i,j"))
                              + dot(v_yy_ab("a,b"), D_mp2_ab("a,b"));
      const double v_zz_mp2 = - dot(v_zz_ij("i,j"), D_mp2_ij("i,j"))
                              + dot(v_zz_ab("a,b"), D_mp2_ab("a,b"));
      const double v_xy_mp2 = - dot(v_xy_ij("i,j"), D_mp2_ij("i,j"))
                              + dot(v_xy_ab("a,b"), D_mp2_ab("a,b"));
      const double v_xz_mp2 = - dot(v_xz_ij("i,j"), D_mp2_ij("i,j"))
                              + dot(v_xz_ab("a,b"), D_mp2_ab("a,b"));
      const double v_yz_mp2 = - dot(v_yz_ij("i,j"), D_mp2_ij("i,j"))
                              + dot(v_yz_ab("a,b"), D_mp2_ab("a,b"));

      const double v_xx_mp2or = dot(v_xx_am("a,m"), Dbn_mp2("a,m"));
      const double v_yy_mp2or = dot(v_yy_am("a,m"), Dbn_mp2("a,m"));
      const double v_zz_mp2or = dot(v_zz_am("a,m"), Dbn_mp2("a,m"));
      const double v_xy_mp2or = dot(v_xy_am("a,m"), Dbn_mp2("a,m"));
      const double v_xz_mp2or = dot(v_xz_am("a,m"), Dbn_mp2("a,m"));
      const double v_yz_mp2or = dot(v_yz_am("a,m"), Dbn_mp2("a,m"));

      std::cout << std::endl << indent
                << "electric gradient (MP2)" << std::endl
                << indent << "v_xx (MP2) = " << scprintf("%12.10f", - v_xx_mp2 * 2.0)
                << indent << "  v_yy (MP2) = " << scprintf("%12.10f", - v_yy_mp2 * 2.0)
                << indent << "  v_zz (MP2) = " << scprintf("%12.10f", - v_zz_mp2 * 2.0)
                << std::endl
                << indent << "v_xy (MP2) = " << scprintf("%12.10f", - v_xy_mp2 * 2.0)
                << indent << "  v_xz (MP2) = " << scprintf("%12.10f", - v_xz_mp2 * 2.0)
                << indent << "  v_yz (MP2) = " << scprintf("%12.10f", - v_yz_mp2 * 2.0)
                << std::endl;
      std::cout << std::endl << indent
                << "electric gradient (MP2 orbital response)" << std::endl
                << indent << "v_xx (MP2 or) = " << scprintf("%12.10f", - v_xx_mp2or * 2.0)
                << indent << "  v_yy (MP2 or) = " << scprintf("%12.10f", - v_yy_mp2or * 2.0)
                << indent << "  v_zz (MP2 or) = " << scprintf("%12.10f", - v_zz_mp2or * 2.0)
                << std::endl
                << indent << "v_xy (MP2 or) = " << scprintf("%12.10f", - v_xy_mp2or * 2.0)
                << indent << "  v_xz (MP2 or) = " << scprintf("%12.10f", - v_xz_mp2or * 2.0)
                << indent << "  v_yz (MP2 or) = " << scprintf("%12.10f", - v_yz_mp2or * 2.0)
                << std::endl << std::endl;
    }

    }
    world_.gop.fence();
#endif

    // compute integrals needed for F12 related calculations of properties
    TArray2 mu_z_apb;
    if (compute_dipole) {
      mu_z_apb = xy("<a'|mu_z|b>");
    }

    TArray2 q_xx_apb, q_yy_apb, q_zz_apb;
    TArray2 Qxx_apb, Qyy_apb, Qzz_apb, Qxy_apb, Qxz_apb, Qyz_apb;
    if (compute_quadrupole) {
      q_xx_apb = xy("<a'|q_xx|b>");
      q_yy_apb = xy("<a'|q_yy|b>");
      q_zz_apb = xy("<a'|q_zz|b>");

      Qxx_apb("a',b") = q_xx_apb("a',b") - (q_zz_apb("a',b") + q_yy_apb("a',b")) * 0.5;
      Qyy_apb("a',b") = q_yy_apb("a',b") - (q_zz_apb("a',b") + q_xx_apb("a',b")) * 0.5 ;
      Qzz_apb("a',b") = q_zz_apb("a',b") - (q_xx_apb("a',b") + q_yy_apb("a',b")) * 0.5;
      Qxy_apb("a',b") = _2("<a'|q_xy|b>") * 1.5;
      Qxz_apb("a',b") = _2("<a'|q_xz|b>") * 1.5;
      Qyz_apb("a',b") = _2("<a'|q_yz|b>") * 1.5;
    }

    TArray2 v_xx_apb, v_yy_apb, v_zz_apb, v_xy_apb, v_xz_apb, v_yz_apb;
    if (compute_EFG) {
      v_xx_apb = xy("<a'|ddphi_xx|b>");
      v_yy_apb = xy("<a'|ddphi_yy|b>");
      v_zz_apb = xy("<a'|ddphi_zz|b>");
      v_xy_apb = xy("<a'|ddphi_xy|b>");
      v_xz_apb = xy("<a'|ddphi_xz|b>");
      v_yz_apb = xy("<a'|ddphi_yz|b>");
    }

    // singlet and triplet coefficients for F12 and coupling terms
    const double C_0 = 1.0 / 2.0;
    const double C_1 = 1.0 / 4.0;
    // compute coefficients needed in the F12 and coupling calculations
    const double R_C1 = (0.5 * C_0 + 1.5 * C_1);
    const double R_C2 = (0.5 * C_0 - 1.5 * C_1);
    const double RR_C1 = 0.5 * C_0 * C_0 + 1.5 * C_1 * C_1;
    const double RR_C2 = 0.5 * C_0 * C_0 - 1.5 * C_1 * C_1;

    // MP2 F12 coupling part
    double mu_z_mp2f12C = 0.0, mu_z_mp2f12Cor = 0.0;
    double q_xx_mp2f12C = 0.0, q_yy_mp2f12C = 0.0, q_zz_mp2f12C = 0.0;
    double q_xx_mp2f12Cor = 0.0, q_yy_mp2f12Cor = 0.0, q_zz_mp2f12Cor = 0.0;
#if INCLUDE_MP2F12C_CONTRI
    {
    ExEnv::out0() << std::endl << indent
                  << "Compute MP2 F12 coupling contributions" << std::endl;
    TArray4 C_ijab, A_ijab;
    C_ijab("i,j,a,b") = _4("<i j|r|a_F(a') b>") + _4("<i j|r|a b_F(a')>");
    A_ijab("i,j,a,b") = C_ijab("i,j,a,b") * Delta_ijab("i,j,a,b");

    // Coupling and F12 part of MP2F12 density
    ExEnv::out0() << indent
                  << "Compute density from coupling contribution" << std::endl;
    TArray2 D_mp2f12_ij, D_mp2f12_ab;
    D_mp2f12_ij("i,j") =  (R_C1 * A_ijab("i,k,a,b") + R_C2 * A_ijab("k,i,a,b"))
                          * T2_ijab("j,k,a,b")
                        // A_jk^ab T_ab^ik
                        + (R_C1 * A_ijab("j,k,a,b") + R_C2 * A_ijab("k,j,a,b"))
                          * T2_ijab("i,k,a,b")
                        // A^ik_ab A^ab_jk
                        + (RR_C1 * A_ijab("i,k,a,b") + RR_C2 * A_ijab("k,i,a,b"))
                          * A_ijab("j,k,a,b");
    D_mp2f12_ab("a,b") =  (R_C1 * T2_ijab("i,j,a,c") + R_C2 * T2_ijab("i,j,c,a"))
                          * A_ijab("i,j,b,c")
                        // A^ac_ij T^ij_bc
                        + (R_C1 * T2_ijab("i,j,b,c") + R_C2 * T2_ijab("i,j,c,b"))
                          * A_ijab("i,j,a,c")
                        //  A_ac^ij A_ij^bc
                        + (RR_C1 * A_ijab("i,j,a,c") + RR_C2 * A_ijab("i,j,c,a"))
                          * A_ijab("i,j,b,c");

    // 1/2 R^kl_a'c \tilde{T}^bc_kl
    TArray2 RT_apb;
    RT_apb("a',b") = _4("<a' c|r|k l>")
                     * ( R_C1 * T2_ijab("k,l,b,c") + R_C2 * T2_ijab("k,l,c,b")
                       + RR_C1 * A_ijab("k,l,b,c") + RR_C2 * A_ijab("k,l,c,b")
                       );

    // MP2 F12 coupling contribution to orbital response
    ExEnv::out0() << indent
                  << "Compute orbital response from coupling contribution" << std::endl;
    TArray2 Xmp2f12_contri_nfzc = Xam_Cmp2f12(C_0, C_1,T2_ijab, A_ijab,
                                              D_mp2f12_ij, D_mp2f12_ab, RT_apb);

    TArray2 Xmp2f12_contri;
    // frozen-core contribution
    TArray2 Diip_Cmp2;
    if (nocc != naocc) {
      ExEnv::out0() << std::endl << indent
                    << "Include frozen-core contributions for MP2 F12 coupling" << std::endl;
      TArray2 Xiip_CT2_mp2;
      Xiip_CT2_mp2("i,i'") =  // F^a'_b R^i'l_a'c \tilde{T}^bc_il
                              _4("<i' l|r|b_F(a') c>")
                              * (  R_C1 * T2_ijab("i,l,b,c") + R_C2 * T2_ijab("l,i,b,c")
                                 + RR_C1 * A_ijab("i,l,b,c") + RR_C2 * A_ijab("l,i,b,c")
                                )
                             +  _4("<i' l|r|c b_F(a')>")
                              * (  R_C1 * T2_ijab("i,l,c,b") + R_C2 * T2_ijab("l,i,c,b")
                                 + RR_C1 * A_ijab("i,l,c,b") + RR_C2 * A_ijab("l,i,c,b")
                                )
                                // g^i'l_cd A^cd_il
                             + _4("<i' l|g|c d>")
                             * (R_C1 * A_ijab("i,l,c,d") + R_C2 * A_ijab("l,i,c,d"))
                            ; // 2.0 is not multiplied as 4.0 is multiplied in the end

      Diip_Cmp2("i,i'") = - Xiip_CT2_mp2("i,i'") * Delta_ijp_F("i,i'");

      Xmp2f12_contri("a,m") = Xmp2f12_contri_nfzc("a,m") - A_ijpam("i,j',a,m") * Diip_Cmp2("i,j'");

    } else {
      Xmp2f12_contri("a,m") = Xmp2f12_contri_nfzc("a,m");
    }

    TArray2 Dbn_mp2f12(Xmp2f12_contri.get_world(), Xmp2f12_contri.trange());
    auto resnorm_mp2f12 = cg_solver2(Orbital_relaxation_Abnam,
                                     Xmp2f12_contri,
                                     Dbn_mp2f12,
                                     preconditioner,
                                     conv_target);

    // MP2-F12 coupling contribution to dipole
    if (compute_dipole) {
      // MP2 F12 coupling density contribution to dipole
      mu_z_mp2f12C = - dot(mu_z_ij("i,j"), D_mp2f12_ij("i,j"))
                     + dot(mu_z_ab("a,b"), D_mp2f12_ab("a,b"))
                     + dot(mu_z_apb("a',b"), RT_apb("a',b")) * 2.0;

      // MP2-F12 coupling orbital response contribution to dipole
      mu_z_mp2f12Cor = dot(mu_z_am("a,m"), Dbn_mp2f12("a,m"));

      if (nocc != naocc) {
        mu_z_mp2f12Cor +=  dot(mu_z_ijp("i,j'"), Diip_Cmp2("i,j'"));
      }
    }

    // MP2-F12 coupling contribution to quadrupoles
    if (compute_quadrupole) {
      q_xx_mp2f12C = - dot(Qxx_ij("i,j"), D_mp2f12_ij("i,j"))
                     + dot(Qxx_ab("a,b"), D_mp2f12_ab("a,b"))
                     + dot(Qxx_apb("a',b"), RT_apb("a',b")) * 2.0;
      q_yy_mp2f12C = - dot(Qyy_ij("i,j"), D_mp2f12_ij("i,j"))
                     + dot(Qyy_ab("a,b"), D_mp2f12_ab("a,b"))
                     + dot(Qyy_apb("a',b"), RT_apb("a',b")) * 2.0;
      q_zz_mp2f12C = - dot(Qzz_ij("i,j"), D_mp2f12_ij("i,j"))
                     + dot(Qzz_ab("a,b"), D_mp2f12_ab("a,b"))
                     + dot(Qzz_apb("a',b"), RT_apb("a',b")) * 2.0;
//      // x2, y2, z2 components to quadrupole moment
//      q_xx_mp2f12C = - dot(q_xx_ij("i,j"), D_mp2f12_ij("i,j"))
//                     + dot(q_xx_ab("a,b"), D_mp2f12_ab("a,b"))
//                     + dot(q_xx_apb("a',b"), RT_apb("a',b")) * 2.0;
//      q_yy_mp2f12C = - dot(q_yy_ij("i,j"), D_mp2f12_ij("i,j"))
//                     + dot(q_yy_ab("a,b"), D_mp2f12_ab("a,b"))
//                     + dot(q_yy_apb("a',b"), RT_apb("a',b")) * 2.0;
//      q_zz_mp2f12C = - dot(q_zz_ij("i,j"), D_mp2f12_ij("i,j"))
//                     + dot(q_zz_ab("a,b"), D_mp2f12_ab("a,b"))
//                     + dot(q_zz_apb("a',b"), RT_apb("a',b")) * 2.0;

      q_xx_mp2f12Cor = dot(Qxx_am("a,m"), Dbn_mp2f12("a,m"));
      q_yy_mp2f12Cor = dot(Qyy_am("a,m"), Dbn_mp2f12("a,m"));
      q_zz_mp2f12Cor = dot(Qzz_am("a,m"), Dbn_mp2f12("a,m"));
//      // x2, y2, z2 components to quadrupole moment
//      q_xx_mp2f12Cor = dot(q_xx_am("a,m"), Dbn_mp2f12("a,m"));
//      q_yy_mp2f12Cor = dot(q_yy_am("a,m"), Dbn_mp2f12("a,m"));
//      q_zz_mp2f12Cor = dot(q_zz_am("a,m"), Dbn_mp2f12("a,m"));

      if (nocc != naocc) {
          q_xx_mp2f12Cor += dot(Qxx_ijp("i,j'"), Diip_Cmp2("i,j'"));
          q_yy_mp2f12Cor += dot(Qyy_ijp("i,j'"), Diip_Cmp2("i,j'"));
          q_zz_mp2f12Cor += dot(Qzz_ijp("i,j'"), Diip_Cmp2("i,j'"));
      }
    }

    // MP2-F12 coupling contribution to electric field gradient on 1st nucleus
    if (compute_EFG) {
      const double v_xx_mp2f12C = - dot(v_xx_ij("i,j"), D_mp2f12_ij("i,j"))
                                  + dot(v_xx_ab("a,b"), D_mp2f12_ab("a,b"))
                                  + dot(v_xx_apb("a',b"), RT_apb("a',b")) * 2.0;
      const double v_yy_mp2f12C = - dot(v_yy_ij("i,j"), D_mp2f12_ij("i,j"))
                                  + dot(v_yy_ab("a,b"), D_mp2f12_ab("a,b"))
                                  + dot(v_yy_apb("a',b"), RT_apb("a',b")) * 2.0;
      const double v_zz_mp2f12C = - dot(v_zz_ij("i,j"), D_mp2f12_ij("i,j"))
                                  + dot(v_zz_ab("a,b"), D_mp2f12_ab("a,b"))
                                  + dot(v_zz_apb("a',b"), RT_apb("a',b")) * 2.0;
      const double v_xy_mp2f12C = - dot(v_xy_ij("i,j"), D_mp2f12_ij("i,j"))
                                  + dot(v_xy_ab("a,b"), D_mp2f12_ab("a,b"))
                                  + dot(v_xy_apb("a',b"), RT_apb("a',b")) * 2.0;
      const double v_xz_mp2f12C = - dot(v_xz_ij("i,j"), D_mp2f12_ij("i,j"))
                                  + dot(v_xz_ab("a,b"), D_mp2f12_ab("a,b"))
                                  + dot(v_xz_apb("a',b"), RT_apb("a',b")) * 2.0;
      const double v_yz_mp2f12C = - dot(v_yz_ij("i,j"), D_mp2f12_ij("i,j"))
                                  + dot(v_yz_ab("a,b"), D_mp2f12_ab("a,b"))
                                  + dot(v_yz_apb("a',b"), RT_apb("a',b")) * 2.0;

      const double v_xx_mp2f12Cor = dot(v_xx_am("a,m"), Dbn_mp2f12("a,m"));
      const double v_yy_mp2f12Cor = dot(v_yy_am("a,m"), Dbn_mp2f12("a,m"));
      const double v_zz_mp2f12Cor = dot(v_zz_am("a,m"), Dbn_mp2f12("a,m"));
      const double v_xy_mp2f12Cor = dot(v_xy_am("a,m"), Dbn_mp2f12("a,m"));
      const double v_xz_mp2f12Cor = dot(v_xz_am("a,m"), Dbn_mp2f12("a,m"));
      const double v_yz_mp2f12Cor = dot(v_yz_am("a,m"), Dbn_mp2f12("a,m"));

      std::cout << std::endl << indent
                << "electric gradient (MP2F12 coupling)" << std::endl
                << indent << "v_xx (MP2F12 C) = " << scprintf("%12.10f", v_xx_mp2f12C * 2.0)
                << indent << "  v_yy (MP2F12 C) = " << scprintf("%12.10f", v_yy_mp2f12C * 2.0)
                << indent << "  v_zz (MP2F12 C) = " << scprintf("%12.10f", v_zz_mp2f12C * 2.0)
                << std::endl
                << indent << "v_xy (MP2F12 C) = " << scprintf("%12.10f", v_xy_mp2f12C * 2.0)
                << indent << "  v_xz (MP2F12 C) = " << scprintf("%12.10f", v_xz_mp2f12C * 2.0)
                << indent << "  v_yz (MP2F12 C) = " << scprintf("%12.10f", v_yz_mp2f12C * 2.0)
                << std::endl;
      std::cout << std::endl << indent
                << "electric gradient (MP2F12 coupling orbital response)" << std::endl
                << indent << "v_xx (MP2F12 C or) = " << scprintf("%12.10f", v_xx_mp2f12Cor * 2.0)
                << indent << "  v_yy (MP2F12 C or) = " << scprintf("%12.10f", v_yy_mp2f12Cor * 2.0)
                << indent << "  v_zz (MP2F12 C or) = " << scprintf("%12.10f", v_zz_mp2f12Cor * 2.0)
                << std::endl
                << indent << "v_xy (MP2F12 C or) = " << scprintf("%12.10f", v_xy_mp2f12Cor * 2.0)
                << indent << "  v_xz (MP2F12 C or) = " << scprintf("%12.10f", v_xz_mp2f12Cor * 2.0)
                << indent << "  v_yz (MP2F12 C or) = " << scprintf("%12.10f", v_yz_mp2f12Cor * 2.0)
                << std::endl << std::endl;
    }

    }
    world_.gop.fence();
#endif

    // F12 contribution
    double mu_z_f12 = 0.0, mu_z_f12or = 0.0;
    double q_xx_f12 = 0.0, q_yy_f12 = 0.0, q_zz_f12 = 0.0;
    double q_xx_f12or = 0.0, q_yy_f12or = 0.0, q_zz_f12or = 0.0;
#if INCLUDE_F12_CONTRI
    {
    ExEnv::out0() << std::endl << indent
                  << "Compute F12 contributions" << std::endl;

    ExEnv::out0() << std::endl << indent << "Compute D_f12" << std::endl;
    // F12 density from X and B terms
    TArray2 D_f12_ij, D_f12_ab, D_f12_apbp, D_f12_apb;
    compute_Df12_XB(C_0, C_1, D_f12_ij, D_f12_ab, D_f12_apbp, D_f12_apb);


    // F12 orbital response contribution

    ExEnv::out0() << std::endl << indent << "Compute X and B density contribution to Xam" << std::endl;
    // X and B density contribution to Xam
    TArray2 gdf12_am = Xam_Df12_XB(C_0, C_1, D_f12_ij,D_f12_ab,
                                   D_f12_apbp, D_f12_apb);

    // Xam from F12 (V, X, and B) contributions
    ExEnv::out0() << std::endl << indent << "Compute V Xam" << std::endl;
    TArray2 Xam_Vcontri = Xam_V(C_0, C_1);
    ExEnv::out0() << std::endl << indent << "Compute X Xam" << std::endl;
    TArray2 Xam_Xcontri = Xam_X(C_0, C_1);
    ExEnv::out0() << std::endl << indent << "Compute B Xam" << std::endl;
    TArray2 Xam_Bcontri = Xam_B(C_0, C_1);

    // F12 contribution to orbital response
    TArray2 Xam_f12_nfzc;
    Xam_f12_nfzc("a,m") = 2.0 * (  Xam_Vcontri("a,m")
                                 - Xam_Xcontri("a,m")
                                 + Xam_Bcontri("a,m")
                                 + gdf12_am("a,m")  // F12 density terms
                                );
    TArray2 Xam_f12;
    // frozen-core contribution of F12 part: Xii'
    TArray2 Diip_f12;
    if (nocc != naocc) {
      ExEnv::out0() << std::endl << indent
                      << "Include frozen-core F12 contributions" << std::endl;
      TArray2 Xiip_f12 = Xiip_VBX(C_0, C_1);
      Diip_f12("i,i'") = - Xiip_f12("i,i'") * Delta_ijp_F("i,i'");

      Xam_f12("a,m") = Xam_f12_nfzc("a,m") - A_ijpam("i,j',a,m") * Diip_f12("i,j'");
    } else {
      Xam_f12("a,m") = Xam_f12_nfzc("a,m");
    }

    TArray2 Dbn_f12(Xam_f12.get_world(), Xam_f12.trange());
    auto resnorm_f12 = cg_solver2(Orbital_relaxation_Abnam,
                                  Xam_f12,
                                  Dbn_f12,
                                  preconditioner,
                                  conv_target);

    // print out each contribution in F12 contribution
#if 0
    TArray2 gdf12_X_am;
    gdf12_X_am("a,m") = // X related contribution
                        (2.0 * _4("<a k|g|m l>") - _4("<a k|g|l m>"))
                        * D_f12_ij("k,l");
    TArray2 gdf12_B_am;
    gdf12_B_am("a,m") = // B related contribution
                       - (2.0 * g_abmc("a,b,m,c") - g_abmc("b,a,m,c"))
                         * D_f12_ab("b,c")
                       - (2.0 * _4("<a b'|g|m c'>") - _4("<a b'|g|c' m>"))
                         * D_f12_apbp("b',c'")
                       - (2.0 * _4("<a b'|g|m c>") - _4("<a b'|g|c m>"))
                         * D_f12_apb("b',c")
                       - (2.0 * _4("<a b|g|m c'>") - _4("<a b|g|c' m>"))
                         * D_f12_apb("c',b");

    TArray2 Xam_V_nfzc;
    Xam_V_nfzc("a,m") = 2.0 * Xam_Vcontri("a,m");
    TArray2 Xam_X_nfzc;
    Xam_X_nfzc("a,m") = 2.0 * (- Xam_Xcontri("a,m") + gdf12_X_am("a,m"));
    TArray2 Xam_B_nfzc;
    Xam_B_nfzc("a,m") = 2.0 * (Xam_Bcontri("a,m")+ gdf12_B_am("a,m"));

    TArray2 Xam_V, Xam_X, Xam_B;
    // frozen-core contribution of F12 part: Xii'
    TArray2 D_iip_V, D_iip_X, D_iip_B;
    if (nocc != naocc) {
      ExEnv::out0() << std::endl << indent
                    << "Include frozen-core V, X, and B contributions" << std::endl;
      TArray2 Xiip_V;
      const char* ip = "i'";
      const char* i = "i";
      TArray2 V_ipi = VRk_Sk(ip,i,C_0, C_1);
      TArray2 V_iip = VRk_Sk(i,ip,C_0, C_1);
      Xiip_V("i,i'") = // 2 * (1/2 R^i'k_A'C' g^A'C'_ik + 1/2 R^ik_A'C' g^A'C'_i'k)
                       2.0 * (V_ipi("i',i") + V_iip("i,i'"));
      D_iip_V("i,i'") = - Xiip_V("i,i'") * Delta_ijp_F("i,i'");
      Xam_V("a,m") = Xam_V_nfzc("a,m") - A_ijpam("i,j',a,m") * D_iip_V("i,j'");

      TArray2 Xiip_X;
      TArray2 F_ij = xy("<i|F|j>");
      TArray4d r2_ipk_jl = ijxy("<i' k|r2|j l>");
      TArray4d r_ipk_pq = ijxy("<i' k|r|p q>");
      TArray4d r_pqkl = ijxy("<p q|r|k l>");
      TArray4d r_ipk_apn = ijxy("<i' k|r|a' n>");
      TArray4d r_apn_kl = ijxy("<a' n|r|k l>");
      TArray4d r_kip_apn = ijxy("<k i'|r|a' n>");
      Xiip_X("i,i'") = // 2 * (F^i_i + F^k_k) 1/2 R^i'k_A'B' R_ik^A'B'
                      - 2.0 * (  (RR_C1 * r2_ipk_jl("i',k,j,l") + RR_C2 * r2_ipk_jl("i',k,l,j"))
                                 * _2("<l|I|k>") * F_ij("j,i")
                               + (RR_C1 * r2_ipk_jl("i',k,i,l") + RR_C2 * r2_ipk_jl("i',k,l,i"))
                                 * F_ij("l,k")

                               - (RR_C1 * r_ipk_pq("i',k,p,q") + RR_C2 * r_ipk_pq("i',k,q,p"))
                                 * ( r_pqkl("p,q,j,k") * F_ij("j,i")
                                   + r_pqkl("p,q,i,l") * F_ij("l,k"))

                               - (RR_C1 * r_ipk_apn("i',k,a',n") + RR_C2 * r_kip_apn("k,i',a',n"))
                                 * ( r_apn_kl("a',n,j,k") * F_ij("j,i")
                                   + r_apn_kl("a',n,i,l") * F_ij("l,k"))

                               - (RR_C1 * r_kip_apn("k,i',a',n") + RR_C2 * r_ipk_apn("i',k,a',n"))
                                 * ( r_apn_kl("a',n,k,j") * F_ij("j,i")
                                   + r_apn_kl("a',n,l,i") * F_ij("l,k"))
                              );
      D_iip_X("i,i'") = - Xiip_X("i,i'") * Delta_ijp_F("i,i'");
      Xam_X("a,m") = Xam_X_nfzc("a,m") - A_ijpam("i,j',a,m") * D_iip_X("i,j'");

      TArray2 Xiip_B;
      TArray2 B_ipi = BPk_Qk(ip,i,C_0,C_1);
      Xiip_B("i,i'") = // 2 R^i'k_A'B' F^A'_C' R^ik_C'B'
                       2.0 * B_ipi("i',i");
      D_iip_B("i,i'") = - Xiip_B("i,i'") * Delta_ijp_F("i,i'");
      Xam_B("a,m") = Xam_B_nfzc("a,m") - A_ijpam("i,j',a,m") * D_iip_B("i,j'");
    } else {
      Xam_V("a,m") = Xam_V_nfzc("a,m");
      Xam_X("a,m") = Xam_X_nfzc("a,m");
      Xam_B("a,m") = Xam_B_nfzc("a,m");
    }

    TArray2 Dbn_V(Xam_V.get_world(), Xam_V.trange());
    auto resnorm_V = cg_solver2(Orbital_relaxation_Abnam,
                                Xam_V,
                                Dbn_V,
                                preconditioner,
                                conv_target);
    TArray2 Dbn_X(Xam_X.get_world(), Xam_X.trange());
    auto resnorm_X = cg_solver2(Orbital_relaxation_Abnam,
                                Xam_X,
                                Dbn_X,
                                preconditioner,
                                conv_target);
    TArray2 Dbn_B(Xam_B.get_world(), Xam_B.trange());
    auto resnorm_B = cg_solver2(Orbital_relaxation_Abnam,
                                Xam_B,
                                Dbn_B,
                                preconditioner,
                                conv_target);

    if (compute_dipole) {
      double mu_z_X, mu_z_B, mu_z_Vor, mu_z_Xor, mu_z_Bor;
      // F12 density contribution to dipole
      mu_z_X = // X related contribution
               - dot(mu_z_ij("i,j"), D_f12_ij("i,j"));
      TArray2 mu_z_apbp = xy("<a'|mu_z|b'>");
      mu_z_B = // B related contribution
                dot(mu_z_ab("a,b"), D_f12_ab("a,b"))
              + dot(mu_z_apbp("a',b'"), D_f12_apbp("a',b'"))
              + dot(mu_z_apb("a',b"), D_f12_apb("a',b")) * 2.0;

      mu_z_Vor = dot(mu_z_am("a,m"), Dbn_V("a,m"));
      mu_z_Xor = dot(mu_z_am("a,m"), Dbn_X("a,m"));
      mu_z_Bor = dot(mu_z_am("a,m"), Dbn_B("a,m"));
      if (nocc != naocc) {
        mu_z_Vor += dot(mu_z_ijp("i,j'"), D_iip_V("i,j'"));
        mu_z_Xor += dot(mu_z_ijp("i,j'"), D_iip_X("i,j'"));
        mu_z_Bor += dot(mu_z_ijp("i,j'"), D_iip_B("i,j'"));
      }

      std::cout << std::endl << indent
                << "mu_z (V orbital response) = " << scprintf("%15.12f", - mu_z_Vor * 2.0)
                << std::endl << indent
                << "mu_z (X) = " << scprintf("%15.12f",- mu_z_X * 2.0)
                << std::endl << indent
                << "mu_z (X orbital response) = " << scprintf("%15.12f", - mu_z_Xor * 2.0)
                << std::endl << indent
                << "mu_z (B) = " << scprintf("%15.12f",- mu_z_B * 2.0)
                << std::endl << indent
                << "mu_z (B orbital response) = " << scprintf("%15.12f", - mu_z_Bor * 2.0)
                << std::endl;
    }
#endif

    if (compute_dipole) {
      TArray2 mu_z_apbp = xy("<a'|mu_z|b'>");
      // F12 density contribution to dipole
      mu_z_f12 = // X related contribution
                 - dot(mu_z_ij("i,j"), D_f12_ij("i,j"))
                 // B related contribution
                 + dot(mu_z_ab("a,b"), D_f12_ab("a,b"))
                 + dot(mu_z_apbp("a',b'"), D_f12_apbp("a',b'"))
                 + dot(mu_z_apb("a',b"), D_f12_apb("a',b")) * 2.0
                 ;
      mu_z_f12or = dot(mu_z_am("a,m"), Dbn_f12("a,m"));
      if (nocc != naocc) {
        mu_z_f12or += dot(mu_z_ijp("i,j'"), Diip_f12("i,j'"));
      }
    }

    // F12 contribution to quadrupoles
    if (compute_quadrupole) {
      TArray2 q_xx_apbp = xy("<a'|q_xx|b'>");
      TArray2 q_yy_apbp = xy("<a'|q_yy|b'>");
      TArray2 q_zz_apbp = xy("<a'|q_zz|b'>");

      TArray2 Qxx_apbp, Qyy_apbp, Qzz_apbp;
      Qxx_apbp("a',b'") = q_xx_apbp("a',b'") - (q_zz_apbp("a',b'") + q_yy_apbp("a',b'")) * 0.5;
      Qyy_apbp("a',b'") = q_yy_apbp("a',b'") - (q_zz_apbp("a',b'") + q_xx_apbp("a',b'")) * 0.5;
      Qzz_apbp("a',b'") = q_zz_apbp("a',b'") - (q_xx_apbp("a',b'") + q_yy_apbp("a',b'")) * 0.5;

      q_xx_f12 = - dot(Qxx_ij("i,j"), D_f12_ij("i,j"))
                 + dot(Qxx_ab("a,b"), D_f12_ab("a,b"))
                 + dot(Qxx_apbp("a',b'"), D_f12_apbp("a',b'"))
                 + dot(Qxx_apb("a',b"), D_f12_apb("a',b")) * 2.0
                 ;
      q_yy_f12 = - dot(Qyy_ij("i,j"), D_f12_ij("i,j"))
                 + dot(Qyy_ab("a,b"), D_f12_ab("a,b"))
                 + dot(Qyy_apbp("a',b'"), D_f12_apbp("a',b'"))
                 + dot(Qyy_apb("a',b"), D_f12_apb("a',b")) * 2.0
                 ;
      q_zz_f12 = - dot(Qzz_ij("i,j"), D_f12_ij("i,j"))
                 + dot(Qzz_ab("a,b"), D_f12_ab("a,b"))
                 + dot(Qzz_apbp("a',b'"), D_f12_apbp("a',b'"))
                 + dot(Qzz_apb("a',b"), D_f12_apb("a',b")) * 2.0
                 ;
//      // x2, y2, z2 components to quadrupole moment
//      q_xx_f12 = - dot(q_xx_ij("i,j"), D_f12_ij("i,j"))
//                 + dot(q_xx_ab("a,b"), D_f12_ab("a,b"))
//                 + dot(q_xx_apbp("a',b'"), D_f12_apbp("a',b'"))
//                 + dot(q_xx_apb("a',b"), D_f12_apb("a',b")) * 2.0;
//      q_yy_f12 = - dot(q_yy_ij("i,j"), D_f12_ij("i,j"))
//                 + dot(q_yy_ab("a,b"), D_f12_ab("a,b"))
//                 + dot(q_yy_apbp("a',b'"), D_f12_apbp("a',b'"))
//                 + dot(q_yy_apb("a',b"), D_f12_apb("a',b")) * 2.0;
//      q_zz_f12 = - dot(q_zz_ij("i,j"), D_f12_ij("i,j"))
//                 + dot(q_zz_ab("a,b"), D_f12_ab("a,b"))
//                 + dot(q_zz_apbp("a',b'"), D_f12_apbp("a',b'"))
//                 + dot(q_zz_apb("a',b"), D_f12_apb("a',b")) * 2.0;

      q_xx_f12or = dot(Qxx_am("a,m"), Dbn_f12("a,m"));
      q_yy_f12or = dot(Qyy_am("a,m"), Dbn_f12("a,m"));
      q_zz_f12or = dot(Qzz_am("a,m"), Dbn_f12("a,m"));
//      // x2, y2, z2 components to quadrupole moment
//      q_xx_f12or = dot(q_xx_am("a,m"), Dbn_f12("a,m"));
//      q_yy_f12or = dot(q_yy_am("a,m"), Dbn_f12("a,m"));
//      q_zz_f12or = dot(q_zz_am("a,m"), Dbn_f12("a,m"));

      if (nocc != naocc) {
          q_xx_f12or += dot(Qxx_ijp("i,j'"), Diip_f12("i,j'"));
          q_yy_f12or += dot(Qyy_ijp("i,j'"), Diip_f12("i,j'"));
          q_zz_f12or += dot(Qzz_ijp("i,j'"), Diip_f12("i,j'"));
      }
    }

    if (compute_EFG) {
      TArray2 v_xx_apbp = xy("<a'|ddphi_xx|b'>");
      TArray2 v_yy_apbp = xy("<a'|ddphi_yy|b'>");
      TArray2 v_zz_apbp = xy("<a'|ddphi_zz|b'>");
      TArray2 v_xy_apbp = xy("<a'|ddphi_xy|b'>");
      TArray2 v_xz_apbp = xy("<a'|ddphi_xz|b'>");
      TArray2 v_yz_apbp = xy("<a'|ddphi_yz|b'>");

      // F12 contribution to electric gradient
      const double v_xx_f12 = - dot(v_xx_ij("i,j"), D_f12_ij("i,j"))
                              + dot(v_xx_ab("a,b"), D_f12_ab("a,b"))
                              + dot(v_xx_apbp("a',b'"), D_f12_apbp("a',b'"))
                              + dot(v_xx_apb("a',b"), D_f12_apb("a',b")) * 2.0;
      const double v_yy_f12 = - dot(v_yy_ij("i,j"), D_f12_ij("i,j"))
                              + dot(v_yy_ab("a,b"), D_f12_ab("a,b"))
                              + dot(v_yy_apbp("a',b'"), D_f12_apbp("a',b'"))
                              + dot(v_yy_apb("a',b"), D_f12_apb("a',b")) * 2.0;
      const double v_zz_f12 = - dot(v_zz_ij("i,j"), D_f12_ij("i,j"))
                              + dot(v_zz_ab("a,b"), D_f12_ab("a,b"))
                              + dot(v_zz_apbp("a',b'"), D_f12_apbp("a',b'"))
                              + dot(v_zz_apb("a',b"), D_f12_apb("a',b")) * 2.0;
      const double v_xy_f12 = - dot(v_xy_ij("i,j"), D_f12_ij("i,j"))
                              + dot(v_xy_ab("a,b"), D_f12_ab("a,b"))
                              + dot(v_xy_apbp("a',b'"), D_f12_apbp("a',b'"))
                              + dot(v_xy_apb("a',b"), D_f12_apb("a',b")) * 2.0;
      const double v_xz_f12 = - dot(v_xz_ij("i,j"), D_f12_ij("i,j"))
                              + dot(v_xz_ab("a,b"), D_f12_ab("a,b"))
                              + dot(v_xz_apbp("a',b'"), D_f12_apbp("a',b'"))
                              + dot(v_xz_apb("a',b"), D_f12_apb("a',b")) * 2.0;
      const double v_yz_f12 = - dot(v_yz_ij("i,j"), D_f12_ij("i,j"))
                              + dot(v_yz_ab("a,b"), D_f12_ab("a,b"))
                              + dot(v_yz_apbp("a',b'"), D_f12_apbp("a',b'"))
                              + dot(v_yz_apb("a',b"), D_f12_apb("a',b")) * 2.0;

      const double v_xx_f12or = dot(v_xx_am("a,m"), Dbn_f12("a,m"));
      const double v_yy_f12or = dot(v_yy_am("a,m"), Dbn_f12("a,m"));
      const double v_zz_f12or = dot(v_zz_am("a,m"), Dbn_f12("a,m"));
      const double v_xy_f12or = dot(v_xy_am("a,m"), Dbn_f12("a,m"));
      const double v_xz_f12or = dot(v_xz_am("a,m"), Dbn_f12("a,m"));
      const double v_yz_f12or = dot(v_yz_am("a,m"), Dbn_f12("a,m"));

      std::cout << indent << "electric gradient (F12)" << std::endl
                << indent << "v_xx (F12) = " << scprintf("%12.10f", v_xx_f12 * 2.0)
                << indent << "  v_yy (F12) = " << scprintf("%12.10f", v_yy_f12 * 2.0)
                << indent << "  v_zz (F12) = " << scprintf("%12.10f", v_zz_f12 * 2.0)
                << std::endl
                << indent << "v_xy (F12) = " << scprintf("%12.10f", v_xy_f12 * 2.0)
                << indent << "  v_xz (F12) = " << scprintf("%12.10f", v_xz_f12 * 2.0)
                << indent << "  v_yz (F12) = " << scprintf("%12.10f", v_yz_f12 * 2.0)
                << std::endl;
      std::cout << indent << "electric gradient (F12 orbital response)" << std::endl
                << indent << "v_xx (F12 or) = " << scprintf("%12.10f", v_xx_f12or * 2.0)
                << indent << "  v_yy (F12 or) = " << scprintf("%12.10f", v_yy_f12or * 2.0)
                << indent << "  v_zz (F12 or) = " << scprintf("%12.10f", v_zz_f12or * 2.0)
                << std::endl
                << indent << "v_xy (F12 or) = " << scprintf("%12.10f", v_xy_f12or * 2.0)
                << indent << "  v_xz (F12 or) = " << scprintf("%12.10f", v_xz_f12or * 2.0)
                << indent << "  v_yz (F12 or) = " << scprintf("%12.10f", v_yz_f12or * 2.0)
                << std::endl << std::endl;
    }
    }
#endif

    ExEnv::out0() << std::endl << indent << "Start CCSD computation" << std::endl;
    double mu_z_ccsd = 0.0, mu_z_ccsdor = 0.0;
    double q_xx_ccsd = 0.0, q_yy_ccsd = 0.0, q_zz_ccsd = 0.0;
    double q_xx_ccsdor = 0.0, q_yy_ccsdor = 0.0, q_zz_ccsdor = 0.0;

    TArray2 mu_z_ai;
    TArray2 Qxx_ai, Qyy_ai, Qzz_ai, Qxy_ai, Qxz_ai, Qyz_ai;
#if 0
    // compute CC2 amplitudes using Schaefer's formula
    TArray2 T1_cc2, L1_cc2;
    TArray4 T2_cc2, L2_cc2;

    ExEnv::out0() << std::endl << indent << "Compute CC2 T amplitudes " << std::endl;
    compute_T_cc2(T1_cc2,T2_cc2);

    ExEnv::out0() << std::endl << indent << "Compute CC2 lambda amplitudes " << std::endl;
    compute_lambda_cc2(T1_cc2, T2_cc2, L1_cc2, L2_cc2);

    // compute CC2 density from amplitudes
    TArray2 Dij_cc2, Dab_cc2, Dia_cc2, Dai_cc2;
    ExEnv::out0() << indent << "Compute CC2 density from amplitudes" << std::endl;
    compute_cc2_1rdm_amp(T1_cc2, T2_cc2, L1_cc2, L2_cc2,
                        Dij_cc2, Dab_cc2, Dia_cc2, Dai_cc2);
#endif

#if 0
    // compute CC2 lambda amplitudes using Gauss's formula
    TArray2 L1_cc2_2;
    TArray4 L2_cc2_2;
    ExEnv::out0() << indent << "Using Gauss's formula" << std::endl;
    compute_lambda_cc2_2(T1_cc2, T2_cc2, L1_cc2_2, L2_cc2_2);
#endif

    ExEnv::out0() << indent << "Compute CCSD T amplitudes " << std::endl;
    TArray2 T1;
    TArray4 T2;
    compute_T_ccsd(T1, T2, "CCSD");

    ExEnv::out0() << std::endl << indent << "Compute CCSD L amplitudes " << std::endl;
    TArray2 L1;
    TArray4 L2;
    compute_lambda_ccsd(T1, T2, L1, L2, "CCSD"); // do not include F12 contributions

#if INCLUDE_CCSD_CONTRI
    {
    // compute CCSD density from amplitudes
    TArray2 Dij_ccsd, Dab_ccsd, Dia_ccsd, Dai_ccsd;
    ExEnv::out0() << std::endl << indent
                  << "Compute CCSD density from amplitudes" << std::endl;
    compute_ccsd_1rdm_amp(T1, T2, L1, L2,
                          Dij_ccsd, Dab_ccsd, Dia_ccsd, Dai_ccsd);

    ExEnv::out0() << std::endl << indent << "Compute CCSD Xam and Xai" << std::endl;
    TArray2 Xam_ccsd_nfzc, Xiip_ccsd;
    compute_Xam_ccsd(T1, T2, L1, L2, Xam_ccsd_nfzc, Xiip_ccsd);

    TArray2 Diip_or_ccsd, Xam_ccsd;
    if (nocc != naocc) {
      ExEnv::out0() << std::endl << indent
                      << "Include frozen-core CCSD contributions" << std::endl;
      Diip_or_ccsd("i,i'") = Xiip_ccsd("i,i'") * Delta_ijp_F("i,i'");

      TArray2 Xam_Diip_ccsd;
      Xam_Diip_ccsd("a,m") = A_ijpam("i,i',a,m") * Diip_or_ccsd("i,i'");
      Xam_ccsd("a,m") = Xam_ccsd_nfzc("a,m") - Xam_Diip_ccsd("a,m");

      // test code: print out frozen-core contribution to Ima
//      TArray2 I_ma_fzc;
//      I_ma_fzc("a,m") = - Xam_Diip_ccsd("a,m") * 0.5;
//      std::cout << "Ima frozen-core part: " << std::endl << I_ma_fzc << std::endl;
    } else {
      Xam_ccsd("a,m") = Xam_ccsd_nfzc("a,m");
    }

    // solve the Z-vector equation for CCSD
    TArray2 Dbn_ccsd(Xam_ccsd.get_world(), Xam_ccsd.trange());
    auto resnorm_ccsd = cg_solver2(Orbital_relaxation_Abnam,
                                   Xam_ccsd,
                                   Dbn_ccsd,
                                   preconditioner,
                                   conv_target);

    if (compute_dipole) {
      mu_z_ai = xy("<a|mu_z|i>");
      mu_z_ccsd =  dot(mu_z_ij("i,j"), Dij_ccsd("i,j"))
                 + dot(mu_z_ab("a,b"), Dab_ccsd("a,b"))
                 + dot(mu_z_ai("a,i"), Dai_ccsd("a,i"))
                 + dot(mu_z_ai("a,i"), Dia_ccsd("i,a"))
                 ;

      mu_z_ccsdor = dot(mu_z_am("a,m"), Dbn_ccsd("a,m"));
      if (nocc != naocc) {
        mu_z_ccsdor += dot(mu_z_ijp("i,i'"), Diip_or_ccsd("i,i'"));
      }

//      std::cout << std::endl << indent
//                  << "mu_z (CCSD) = " << scprintf("%15.12f", - mu_z_ccsd * 2.0)
//                  << std::endl << indent
//                  << "mu_z (CCSD orbital response) = " << scprintf("%15.12f", - mu_z_ccsdor * 2.0)
//                  << std::endl;
    }

    if (compute_quadrupole) {
      TArray2 q_xx_ai = xy("<a|q_xx|i>");
      TArray2 q_yy_ai = xy("<a|q_yy|i>");
      TArray2 q_zz_ai = xy("<a|q_zz|i>");
      Qxx_ai("a,i") = q_xx_ai("a,i") - (q_zz_ai("a,i") + q_yy_ai("a,i")) * 0.5;
      Qyy_ai("a,i") = q_yy_ai("a,i") - (q_zz_ai("a,i") + q_xx_ai("a,i")) * 0.5;
      Qzz_ai("a,i") = q_zz_ai("a,i") - (q_xx_ai("a,i") + q_yy_ai("a,i")) * 0.5;

      q_xx_ccsd =  dot(Qxx_ij("i,j"), Dij_ccsd("i,j"))
                 + dot(Qxx_ab("a,b"), Dab_ccsd("a,b"))
                 + dot(Qxx_ai("a,i"), Dai_ccsd("a,i"))
                 + dot(Qxx_ai("a,i"), Dia_ccsd("i,a"))
                 ;
      q_yy_ccsd =  dot(Qyy_ij("i,j"), Dij_ccsd("i,j"))
                 + dot(Qyy_ab("a,b"), Dab_ccsd("a,b"))
                 + dot(Qyy_ai("a,i"), Dai_ccsd("a,i"))
                 + dot(Qyy_ai("a,i"), Dia_ccsd("i,a"))
                 ;
      q_zz_ccsd =  dot(Qzz_ij("i,j"), Dij_ccsd("i,j"))
                 + dot(Qzz_ab("a,b"), Dab_ccsd("a,b"))
                 + dot(Qzz_ai("a,i"), Dai_ccsd("a,i"))
                 + dot(Qzz_ai("a,i"), Dia_ccsd("i,a"))
                 ;

      q_xx_ccsdor = dot(Qxx_am("a,m"), Dbn_ccsd("a,m"));
      q_yy_ccsdor = dot(Qyy_am("a,m"), Dbn_ccsd("a,m"));
      q_zz_ccsdor = dot(Qzz_am("a,m"), Dbn_ccsd("a,m"));

      if (nocc != naocc) {
        q_xx_ccsdor += dot(Qxx_ijp("i,j'"), Diip_or_ccsd("i,j'"));
        q_yy_ccsdor += dot(Qyy_ijp("i,j'"), Diip_or_ccsd("i,j'"));
        q_zz_ccsdor += dot(Qzz_ijp("i,j'"), Diip_or_ccsd("i,j'"));
      }
    }
    }
#endif

    // CC F12 coupling contribution to Xam
    double mu_z_Cccsdf12 = 0.0, mu_z_Cccsdf12or = 0.0;
    double q_xx_Cccsdf12 = 0.0, q_yy_Cccsdf12 = 0.0, q_zz_Cccsdf12 = 0.0;
    double q_xx_Cccsdf12or = 0.0, q_yy_Cccsdf12or = 0.0, q_zz_Cccsdf12or = 0.0;
#if INCLUDE_CCSDF12C_CONTRI
    {
    // resolve CCSD lambda amplitudes, density, and Xam in the framework of CCSD(2)_F12
    ExEnv::out0() << std::endl << indent << "Compute CCSD(2)_F12 L amplitudes " << std::endl;
    TArray2 L1_f12;
    TArray4 L2_f12;
    compute_lambda_ccsd(T1, T2, L1_f12, L2_f12, "F12");

    // compute CCSD density from amplitudes
    TArray2 Dij_ccsd_f12, Dab_ccsd_f12, Dia_ccsd_f12, Dai_ccsd_f12;
    ExEnv::out0() << std::endl << indent
                  << "Compute CCSD density from amplitudes in CCSD(2)_F12" << std::endl;
    compute_ccsd_1rdm_amp(T1, T2, L1_f12, L2_f12,
                          Dij_ccsd_f12, Dab_ccsd_f12, Dia_ccsd_f12, Dai_ccsd_f12);

    ExEnv::out0() << std::endl << indent
                  << "Compute CCSD Xam and Xai in CCSD(2)_F12" << std::endl;
    TArray2 Xam_ccsd_f12_nfzc, Xiip_ccsd_f12;
    compute_Xam_ccsd(T1, T2, L1_f12, L2_f12, Xam_ccsd_f12_nfzc, Xiip_ccsd_f12);
    //compute_Xam_ccsd(T1, T2, L1, L2_f12, Xam_ccsd_f12_nfzc, Xiip_ccsd_f12);

    TArray2 Diip_ccsd_f12, Xam_ccsd_f12;
    if (nocc != naocc) {
      Diip_ccsd_f12("i,i'") = Xiip_ccsd_f12("i,i'") * Delta_ijp_F("i,i'");
      Xam_ccsd_f12("a,m") = Xam_ccsd_f12_nfzc("a,m") - A_ijpam("i,i',a,m") * Diip_ccsd_f12("i,i'");
    } else {
      Xam_ccsd_f12("a,m") = Xam_ccsd_f12_nfzc("a,m");
    }

    TArray2 Dbn_ccsd_f12(Xam_ccsd_f12.get_world(), Xam_ccsd_f12.trange());
    auto resnorm_ccsd_f12 = cg_solver2(Orbital_relaxation_Abnam,
                                       Xam_ccsd_f12,
                                       Dbn_ccsd_f12,
                                       preconditioner,
                                       1e-12);

    ExEnv::out0() << std::endl << indent
                  << "Compute CCSD F12 coupling contributions" << std::endl;
    // CC CT2
    // 1/2 R^kl_a'c T2^bc_kl
    TArray2 RT2_aPb;
    RT2_aPb("a',b") = _4("<a' c|r|k l>")
                      * (R_C1 * T2("b,c,k,l") + R_C2 * T2("c,b,k,l"));

    TArray2 Xam_CT2_cc = Xam_CT2_ccsd(C_0, C_1, T2, RT2_aPb);
    // VT1 & VT2 coupling contribution to F12 Xam
    TArray2 Xam_VT_cc = Xam_VT_ccsd(C_0, C_1, T1, T2);
    TArray2 Xam_CVT_cc_nfzc;
    Xam_CVT_cc_nfzc("a,m") = Xam_CT2_cc("a,m") + Xam_VT_cc("a,m");

    // frozen-core contribution
    TArray2 Diip_CVT_cc, Xam_CVT_cc;
    if (nocc != naocc) {
      ExEnv::out0() << std::endl << indent
                      << "Include frozen-core contributions for CCSD F12 coupling" << std::endl;
      TArray2 Xiip_couling = Xiip_CVT(C_0, C_1, T1, T2);
      Diip_CVT_cc("i,i'") = - Xiip_couling("i,i'") * Delta_ijp_F("i,i'");

      TArray2 Xam_Diip_CVT;
      Xam_Diip_CVT("a,m") = - A_ijpam("i,j',a,m") * Diip_CVT_cc("i,j'");

      Xam_CVT_cc("a,m") = Xam_CVT_cc_nfzc("a,m") + Xam_Diip_CVT("a,m");
    } else {
      Xam_CVT_cc("a,m") = Xam_CVT_cc_nfzc("a,m");
    }

    TArray2 Dbn_CVT_cc(Xam_CVT_cc.get_world(), Xam_CVT_cc.trange());
    auto resnorm_CVT_cc = cg_solver2(Orbital_relaxation_Abnam,
                                    Xam_CVT_cc,
                                    Dbn_CVT_cc,
                                    preconditioner,
                                    conv_target);

    if (compute_dipole) {
      double mu_z_ccsd_f12 =  dot(mu_z_ij("i,j"), Dij_ccsd_f12("i,j"))
                            + dot(mu_z_ab("a,b"), Dab_ccsd_f12("a,b"))
                            + dot(mu_z_ai("a,i"), Dai_ccsd_f12("a,i"))
                            + dot(mu_z_ai("a,i"), Dia_ccsd_f12("i,a"))
                            // F12 coupling density contribution
                            + dot(mu_z_apb("a',b"), RT2_aPb("a',b")) * 2.0
                            ;

      double mu_z_ccsdor_f12 =  dot(mu_z_am("a,m"), Dbn_ccsd_f12("a,m"))
                              + dot(mu_z_am("a,m"), Dbn_CVT_cc("a,m"));

      if (nocc != naocc) {
        mu_z_ccsdor_f12 +=  dot(mu_z_ijp("i,j'"), Diip_ccsd_f12("i,j'"))
                          + dot(mu_z_ijp("i,j'"), Diip_CVT_cc("i,j'"));
      }

      mu_z_Cccsdf12 = mu_z_ccsd_f12 - mu_z_ccsd;
      mu_z_Cccsdf12or = mu_z_ccsdor_f12 - mu_z_ccsdor;
    }

    if (compute_quadrupole) {
      double q_xx_ccsd_f12 =  dot(Qxx_ij("i,j"), Dij_ccsd_f12("i,j"))
                            + dot(Qxx_ab("a,b"), Dab_ccsd_f12("a,b"))
                            + dot(Qxx_ai("a,i"), Dai_ccsd_f12("a,i"))
                            + dot(Qxx_ai("a,i"), Dia_ccsd_f12("i,a"))
                            + dot(Qxx_apb("a',b"), RT2_aPb("a',b")) * 2.0
                            ;
      double q_yy_ccsd_f12 =  dot(Qyy_ij("i,j"), Dij_ccsd_f12("i,j"))
                            + dot(Qyy_ab("a,b"), Dab_ccsd_f12("a,b"))
                            + dot(Qyy_ai("a,i"), Dai_ccsd_f12("a,i"))
                            + dot(Qyy_ai("a,i"), Dia_ccsd_f12("i,a"))
                            + dot(Qyy_apb("a',b"), RT2_aPb("a',b")) * 2.0
                            ;
      double q_zz_ccsd_f12 =  dot(Qzz_ij("i,j"), Dij_ccsd_f12("i,j"))
                            + dot(Qzz_ab("a,b"), Dab_ccsd_f12("a,b"))
                            + dot(Qzz_ai("a,i"), Dai_ccsd_f12("a,i"))
                            + dot(Qzz_ai("a,i"), Dia_ccsd_f12("i,a"))
                            + dot(Qzz_apb("a',b"), RT2_aPb("a',b")) * 2.0
                            ;

      double q_xx_ccsdor_f12 =  dot(Qxx_am("a,m"), Dbn_ccsd_f12("a,m"))
                              + dot(Qxx_am("a,m"), Dbn_CVT_cc("a,m"));
      double q_yy_ccsdor_f12 =  dot(Qyy_am("a,m"), Dbn_ccsd_f12("a,m"))
                              + dot(Qyy_am("a,m"), Dbn_CVT_cc("a,m"));
      double q_zz_ccsdor_f12 =  dot(Qzz_am("a,m"), Dbn_ccsd_f12("a,m"))
                              + dot(Qzz_am("a,m"), Dbn_CVT_cc("a,m"));

      if (nocc != naocc) {
          q_xx_ccsdor_f12 +=  dot(Qxx_ijp("i,j'"), Diip_ccsd_f12("i,j'"))
                            + dot(Qxx_ijp("i,j'"), Diip_CVT_cc("i,j'"));
          q_yy_ccsdor_f12 +=  dot(Qyy_ijp("i,j'"), Diip_ccsd_f12("i,j'"))
                            + dot(Qyy_ijp("i,j'"), Diip_CVT_cc("i,j'"));
          q_zz_ccsdor_f12 +=  dot(Qzz_ijp("i,j'"), Diip_ccsd_f12("i,j'"))
                            + dot(Qzz_ijp("i,j'"), Diip_CVT_cc("i,j'"));
      }

      q_xx_Cccsdf12 = q_xx_ccsd_f12 - q_xx_ccsd;
      q_yy_Cccsdf12 = q_yy_ccsd_f12 - q_yy_ccsd;
      q_zz_Cccsdf12 = q_zz_ccsd_f12 - q_zz_ccsd;
      q_xx_Cccsdf12or = q_xx_ccsdor_f12 - q_xx_ccsdor;
      q_yy_Cccsdf12or = q_yy_ccsdor_f12 - q_yy_ccsdor;
      q_zz_Cccsdf12or = q_zz_ccsdor_f12 - q_zz_ccsdor;
    }

    }
#endif

    // print results
    if (compute_dipole) {
      std::cout << std::endl << indent << "Dipole moment "
                << std::endl << indent
                << "mu_z (HF=SCF+N) = " << scprintf("%15.12f", mu_z_n - mu_z_scf * 2.0)
                                                //electron charge = -1, hence the minus
                << std::endl << std::endl << indent
                << "mu_z (E2) = " << scprintf("%15.12f", - mu_z_e2 * 2.0)
                << std::endl << indent
                << "mu_z (E2 orbital response) = " << scprintf("%15.12f", - mu_z_e2or * 2.0)
                << std::endl << std::endl << indent
                << "mu_z (MP2) = " << scprintf("%15.12f", - mu_z_mp2 * 2.0)
                << std::endl << indent
                << "mu_z (MP2 orbital response) = "<< scprintf("%15.12f", - mu_z_mp2or * 2.0)
                << std::endl << std::endl << indent
                << "mu_z (MP2F12 coupling) = " << scprintf("%15.12f", - mu_z_mp2f12C * 2.0)
                << std::endl << indent
                << "mu_z (MP2-F12 coupling orbital response) = " << scprintf("%15.12f", - mu_z_mp2f12Cor * 4.0)
                << std::endl << std::endl << indent
                << "mu_z (F12) = " << scprintf("%15.12f",- mu_z_f12 * 2.0)
                << std::endl << indent
                << "mu_z (F12 orbital response) = " << scprintf("%15.12f", - mu_z_f12or * 2.0)
                << std::endl << std::endl << indent
                << "mu_z (CCSD) = " << scprintf("%15.12f", - mu_z_ccsd * 2.0)
                << std::endl << indent
                << "mu_z (CCSD orbital response) = " << scprintf("%15.12f", - mu_z_ccsdor * 2.0)
                << std::endl << std::endl << indent
                << "mu_z (CCSDF12 coupling) = " << scprintf("%15.12f", - mu_z_Cccsdf12 * 2.0)
                << std::endl << indent
                << "mu_z (CCSDF12 coupling orbital response) = " << scprintf("%15.12f", - mu_z_Cccsdf12or * 2.0)
                << std::endl;
    }

    if (compute_quadrupole) {
      std::cout << std::endl << indent
                << "Traceless quadrupole moment"
                << std::endl << indent
                << "q_xx (HF) = "   << scprintf("%12.10f", q_xx_n - q_xx_scf * 2.0)
                << "  q_yy (HF) = " << scprintf("%12.10f", q_yy_n - q_yy_scf * 2.0)
                << "  q_zz (HF) = " << scprintf("%12.10f", q_zz_n - q_zz_scf * 2.0)
                << std::endl << indent
                << "q_xy (HF) = "   << scprintf("%12.10f", q_xy_n - q_xy_scf * 2.0)
                << "  q_xz (HF) = " << scprintf("%12.10f", q_xz_n - q_xz_scf * 2.0)
                << "  q_yz (HF) = " << scprintf("%12.10f", q_yz_n - q_yz_scf * 2.0)
                << std::endl;

      std::cout << std::endl << indent
                << "q_xx (E2) = "   << scprintf("%12.10f", - q_xx_e2 * 2.0)
                << "  q_yy (E2) = " << scprintf("%12.10f", - q_yy_e2 * 2.0)
                << "  q_zz (E2) = " << scprintf("%12.10f", - q_zz_e2 * 2.0)
                << std::endl << indent
                << "q_xx (E2 or) = "   << scprintf("%12.10f", - q_xx_e2or * 2.0)
                << "  q_yy (E2 or) = " << scprintf("%12.10f", - q_yy_e2or * 2.0)
                << "  q_zz (E2 or) = " << scprintf("%12.10f", - q_zz_e2or * 2.0)
                << std::endl;

      std::cout << std::endl << indent
                << "q_xx (MP2) = "   << scprintf("%12.10f", - q_xx_mp2 * 2.0)
                << "  q_yy (MP2) = " << scprintf("%12.10f", - q_yy_mp2 * 2.0)
                << "  q_zz (MP2) = " << scprintf("%12.10f", - q_zz_mp2 * 2.0)
                << std::endl << indent
                << "q_xx (MP2 or) = "   << scprintf("%12.10f", - q_xx_mp2or * 2.0)
                << "  q_yy (MP2 or) = " << scprintf("%12.10f", - q_yy_mp2or * 2.0)
                << "  q_zz (MP2 or) = " << scprintf("%12.10f", - q_zz_mp2or * 2.0)
                << std::endl;

      std::cout << std::endl << indent
                << "q_xx (MP2F12 C) = "   << scprintf("%12.10f", - q_xx_mp2f12C * 2.0)
                << "  q_yy (MP2F12 C) = " << scprintf("%12.10f", - q_yy_mp2f12C * 2.0)
                << "  q_zz (MP2F12 C) = " << scprintf("%12.10f", - q_zz_mp2f12C * 2.0)
                << std::endl << indent
                << "q_xx (MP2F12 C or) = "   << scprintf("%12.10f", - q_xx_mp2f12Cor * 2.0)
                << "  q_yy (MP2F12 C or) = " << scprintf("%12.10f", - q_yy_mp2f12Cor * 2.0)
                << "  q_zz (MP2F12 C or) = " << scprintf("%12.10f", - q_zz_mp2f12Cor * 2.0)
                << std::endl;

      std::cout << std::endl << indent
                << "q_xx (F12) = "   << scprintf("%12.10f", - q_xx_f12 * 2.0)
                << "  q_yy (F12) = " << scprintf("%12.10f", - q_yy_f12 * 2.0)
                << "  q_zz (F12) = " << scprintf("%12.10f", - q_zz_f12 * 2.0)
                << std::endl << indent
                << "q_xx (F12 or) = "   << scprintf("%12.10f", - q_xx_f12or * 2.0)
                << "  q_yy (F12 or) = " << scprintf("%12.10f", - q_yy_f12or * 2.0)
                << "  q_zz (F12 or) = " << scprintf("%12.10f", - q_zz_f12or * 2.0)
                << std::endl;

      std::cout << std::endl << indent
                << "q_xx (CCSD) = "   << scprintf("%12.10f", - q_xx_ccsd * 2.0)
                << "  q_yy (CCSD) = " << scprintf("%12.10f", - q_yy_ccsd * 2.0)
                << "  q_zz (CCSD) = " << scprintf("%12.10f", - q_zz_ccsd * 2.0)
                << std::endl << indent
                << "q_xx (CCSD or) = "   << scprintf("%12.10f", - q_xx_ccsdor * 2.0)
                << "  q_yy (CCSD or) = " << scprintf("%12.10f", - q_yy_ccsdor * 2.0)
                << "  q_zz (CCSD or) = " << scprintf("%12.10f", - q_zz_ccsdor * 2.0)
                << std::endl;

      std::cout << std::endl << indent
                << "q_xx (CCSDF12 C) = "   << scprintf("%12.10f", - q_xx_Cccsdf12 * 2.0)
                << "  q_yy (CCSDF12 C) = " << scprintf("%12.10f", - q_yy_Cccsdf12 * 2.0)
                << "  q_zz (CCSDF12 C) = " << scprintf("%12.10f", - q_zz_Cccsdf12 * 2.0)
                << std::endl << indent
                << "q_xx (CCSDF12 C or) = "   << scprintf("%12.10f", - q_xx_Cccsdf12or * 2.0)
                << "  q_yy (CCSDF12 C or) = " << scprintf("%12.10f", - q_yy_Cccsdf12or * 2.0)
                << "  q_zz (CCSDF12 C or) = " << scprintf("%12.10f", - q_zz_Cccsdf12or * 2.0)
                << std::endl << std::endl;
    }
  }

  // F12b method
  // Xam contribution of CT2 part resulted from CCSD F12 coupling
  template <typename T>
  typename SingleReference_R12Intermediates<T>::TArray2
  SingleReference_R12Intermediates<T>::Xam_CT2L2_f12b(const double C_0, const double C_1,
                                                      const TArray4& T2, const TArray2& RT2_aPb,
                                                      const TArray4& L2, const TArray2& RL2_aPb) {
    const double R_C1 = (0.5 * C_0 + 1.5 * C_1);
    const double R_C2 = (0.5 * C_0 - 1.5 * C_1);

    TArray4d g_abmaP = ijxy("<a b|g|m a'>");
    TArray2 Xam_CT2, Xai_CT2;
    Xam_CT2("a,m") = // + 1/2 F^a'_m R^kl_a'b (T2^ab_kl + L2^ab_kl)
                     // + 1/2 F^a'_b R^kl_ma' (T2^ab_kl + L2^ab_kl)
                     (  R_C1 * T2("a,b,k,l") + R_C2 * T2("b,a,k,l")
                      + R_C1 * L2("a,b,k,l") + R_C2 * L2("b,a,k,l")
                     ) * (_4("<k l|r|m_F(a') b>") + _4("<k l|r|m b_F(a')>"))

                      // - g^aa'_mb 1/2 R^kl_a'c (T2^bc_kl + L2^bc_kl)
                      // - g^ab_ma' 1/2 R^kl_a'c (T2^bc_kl + L2^bc_kl)
                    - (  2.0 * _4("<a a'|g|m b>") - _4("<a a'|g|b m>")
                       + 2.0 * g_abmaP("a,b,m,a'") - g_abmaP("b,a,m,a'")
                      ) * (RT2_aPb("a',b") + RL2_aPb("a',b")
                          )
                    ;
    Xai_CT2("a,i") =  // - F^a'_b R^al_a'c (T2^bc_il + L2^bc_il)
                    - _4("<a l|r|b_F(a') c>")
                      * (  R_C1 * T2("b,c,i,l") + R_C2 * T2("b,c,l,i")
                         + R_C1 * L2("b,c,i,l") + R_C2 * L2("b,c,l,i")
                        )
                    - _4("<a l|r|c b_F(a')>")
                      * (  R_C1 * T2("c,b,i,l") + R_C2 * T2("c,b,l,i")
                         + R_C1 * L2("c,b,i,l") + R_C2 * L2("c,b,l,i")
                        )
                    ;

    return XaiAddToXam(Xam_CT2, Xai_CT2);
  }

  // F12b method
  // Xam contribution of VT(T1&T2) part resulted from CCSD F12 coupling
  template <typename T>
  typename SingleReference_R12Intermediates<T>::TArray2
  SingleReference_R12Intermediates<T>::Xam_VTL_f12b(const double C_0, const double C_1,
                                                   const TArray2& T1, const TArray4& T2,
                                                   const TArray2& L1, const TArray4& L2) {

    const double R_C1 = (0.5 * C_0 + 1.5 * C_1);
    const double R_C2 = (0.5 * C_0 - 1.5 * C_1);

    const char* i = "i";
    const char* j = "j";
    const char* k = "k";
    const char* l = "l";
    const char* m = "m";

    const char* a = "a";
    const char* c = "c";
    const char* d = "d";

    TArray4 V_alic = VPq_Rs(a,l,i,c, C_0, C_1);
    TArray4 V_ilac = VPq_Rs(i,l,a,c,C_0, C_1);
    TArray2 V_ac = VRk_Sk(a,c,C_0, C_1);
    TArray2 V_jm = VRk_Sk(j,m,C_0, C_1);
    TArray4d r_klmaP = ijxy("<k l|r|m a'>");
    TArray4d r_abPkl = ijxy("<a b'|r|k l>");
    TArray2 Xam_VT1L1, Xai_VT1L1;
    Xam_VT1L1("a,m") =  // + 1/2 R^jk_A'B' g^A'B'_mk (t1^a_j + l1^a_j)
                        (T1("a,j") + L1("a,j")) * V_jm("j,m")

                      + (  // + R^kl_mb' g^ab'_kc (t1^c_l + l1^c_l)
                           _4("<a b'|g|k c>")
                           * (R_C1 * r_klmaP("k,l,m,b'") + R_C2 * r_klmaP("l,k,m,b'"))
                         + _4("<a b'|g|c k>")
                           * (R_C1 * r_klmaP("l,k,m,b'") + R_C2 * r_klmaP("k,l,m,b'"))
                         // + R^kl_ab' g^mb'_kc (t1^c_l + l1^c_l)
                         + (R_C1 * r_abPkl("a,b',k,l") + R_C2 * r_abPkl("a,b',l,k"))
                           * _4("<k c|g|m b'>")
                         + (R_C1 * r_abPkl("a,b',l,k") + R_C2 * r_abPkl("a,b',k,l"))
                           * _4("<c k|g|m b'>")
                        ) * (T1("c,l") + L1("c,l"))
                      ;
    Xai_VT1L1("a,i") =  // - 1/2 R^al_A'B' g^A'B'_ic (t1^c_l + l1^c_l)
                        // - 1/2 R^il_A'B' g^A'B'_ac (t1^c_l + l1^c_l)
                      - (V_alic("a,l,i,c") + V_ilac("i,l,a,c"))
                        * (T1("c,l") + L1("c,l"))

                        // - 1/2 R^ak_A'B' g^A'B'_ck (t1^c_i + l1^c_i)
                      - V_ac("a,c") * (T1("c,i") + L1("c,i"))
                      ;
    TArray2 Xam_VT1L1_tot = XaiAddToXam(Xam_VT1L1, Xai_VT1L1);

    TArray4 V_alcd = VPq_Rs(a,l,c,d,C_0, C_1);
    TArray4 V_klmd = VPq_Rs(k,l,m,d,C_0, C_1);
    TArray2 Xam_VT2L2, Xai_VT2L2;
    Xam_VT2L2("a,m") =  //   1/4 R^kl_A'B' g^A'B'_md (T^ad_kl + L^ad_kl)
                        (T2("a,d,k,l") + L2("a,d,k,l")) * V_klmd("k,l,m,d")

                      + (  // + 1/4 R^kl_ab' g^mb'_cd (T^cd_kl + L^cd_kl)
                           (R_C1 * r_abPkl("a,b',k,l") + R_C2 * r_abPkl("a,b',l,k"))
                           * _4("<c d|g|m b'>")
                           // + 1/4 R^kl_mb' g^ab'_cd (T^cd_kl + L^cd_kl)
                         + _4("<a b'|g|c d>")
                           * (R_C1 * r_klmaP("k,l,m,b'") + R_C2 * r_klmaP("l,k,m,b'"))
                        ) * (T2("c,d,k,l") + L2("c,d,k,l"))
                      ;
    Xai_VT2L2("a,i") =  // - 1/4 R^al_A'B' g^A'B'_cd (T^cd_il + L^cd_il)
                      - V_alcd("a,l,c,d") * (T2("c,d,i,l") + L2("c,d,i,l"))
                      ;
    TArray2 Xam_VT2L2_tot = XaiAddToXam(Xam_VT2L2, Xai_VT2L2);

    TArray2 Xam_VTL_tot;
    Xam_VTL_tot("a,m") = Xam_VT1L1_tot("a,m") + Xam_VT2L2_tot("a,m");
    return Xam_VTL_tot;
  }

  // F12b method
  // Xam contribution of VT1T1 part resulted from CCSD F12 coupling
  // VT1T1: 1/2 [V^dagger]^cd_kl t^k_c t^l_d
  //      (1/4 R^AB_kl g^cd_AB t^k_c t^l_d)
  template <typename T>
  typename SingleReference_R12Intermediates<T>::TArray2
  SingleReference_R12Intermediates<T>::Xam_VT1T1_f12b(const double C_0, const double C_1,
                                                      const TArray2& T1) {
    const char* k = "k";
    const char* l = "l";
    const char* m = "m";

    const char* a = "a";
    const char* c = "c";
    const char* d = "d";

    const double R_C1 = (0.5 * C_0 + 1.5 * C_1);
    const double R_C2 = (0.5 * C_0 - 1.5 * C_1);

    TArray4 Vklmd = VPq_Rs(k, l, m, d, C_0, C_1); // V^Rk_Sk = 1/2 bar{R}^Rk_A'B' bar{g}^A'B'_Sk
    TArray4 Valcd = VPq_Rs(a, l, c, d, C_0, C_1);

    // test for VPQ_RS
    //TArray4 V_KLMD = VPQ_RS(k, l, m, d);

    TArray4d rmbp_kl = ijxy("<m b'|r|k l>");
    TArray4d rabp_kl = ijxy("<a b'|r|k l>");

    TArray2 Xam_VT1T1;
    Xam_VT1T1("a,m") =  // 1/2 R^AB_kl g^md_AB t^k_a t^l_d
                        Vklmd("k,l,m,d")
                        //(R_C1 * V_KLMD("k,l,m,d") + R_C2 * V_KLMD("l,k,m,d"))
                        * T1("a,k") * T1("d,l")

                        // 1/2 R^mb'_kl g^cd_ab' t^k_c t^l_d
                        // 1/2 R^ab'_kl g^cd_mb' t^k_c t^l_d
                       + (  (R_C1 * rmbp_kl("m,b',k,l") + R_C2 * rmbp_kl("m,b',l,k"))
                            * _4("<a b'|g|c d>")
                          + (R_C1 * rabp_kl("a,b',k,l") + R_C2 * rabp_kl("a,b',l,k"))
                             * _4("<m b'|g|c d>")
                         ) * T1("c,k") * T1("d,l")
                       ;
    TArray2 Xai_VT1T1;
    Xai_VT1T1("a,i") =  // - 1/2 R^AB_al g^cd_AB t^i_c t^l_d
                       - Valcd("a,l,c,d") * T1("c,i") * T1("d,l")
                       ;

    return XaiAddToXam(Xam_VT1T1, Xai_VT1T1);
  }

  template <typename T>
  typename SingleReference_R12Intermediates<T>::TArray2
  SingleReference_R12Intermediates<T>::Xam_VT1T1_f12b_test(const double C_0, const double C_1,
                                                           const TArray4& tauT1_ab) {

    const double R_C1 = (0.5 * C_0 + 1.5 * C_1);
    const double R_C2 = (0.5 * C_0 - 1.5 * C_1);

    const char* k = "k";
    const char* l = "l";
    const char* m = "m";

    const char* a = "a";
    const char* c = "c";
    const char* d = "d";

    TArray4d r_klmaP = ijxy("<k l|r|m a'>");
    TArray4d r_abPkl = ijxy("<a b'|r|k l>");


    TArray4 V_alcd = VPq_Rs(a,l,c,d,C_0, C_1);
    TArray4 V_klmd = VPq_Rs(k,l,m,d,C_0, C_1);
    TArray2 Xam, Xai;
    Xam("a,m") =  //
                  tauT1_ab("a,d,k,l") * V_klmd("k,l,m,d")

                + (  //
                     (R_C1 * r_abPkl("a,b',k,l") + R_C2 * r_abPkl("a,b',l,k"))
                     * _4("<c d|g|m b'>")
                     //
                   + _4("<a b'|g|c d>")
                     * (R_C1 * r_klmaP("k,l,m,b'") + R_C2 * r_klmaP("l,k,m,b'"))
                  ) * tauT1_ab("c,d,k,l")
                ;
    Xai("a,i") =  //
                - V_alcd("a,l,c,d") * tauT1_ab("c,d,i,l")
                ;

    return XaiAddToXam(Xam, Xai);
  }

  // F12b method
  // Xam contribution of VL2T1 part resulted from CCSD F12 coupling
  // VL2T1: 1/2 L^cd_kl V^cd_jc t^j_d
  //     (- 1/4 L^cd_kl R^kl_AB g^AB_cj t^j_d)
  template <typename T>
  typename SingleReference_R12Intermediates<T>::TArray2
  SingleReference_R12Intermediates<T>::Xam_VL2T1_f12b(const double C_0, const double C_1,
                                                      const TArray2& T1, const TArray4& L2) {

    const char* k = "k";
    const char* l = "l";
    const char* j = "j";
    const char* m = "m";
    const char* a = "a";
    const char* c = "c";

    const double R_C1 = (0.5 * C_0 + 1.5 * C_1);
    const double R_C2 = (0.5 * C_0 - 1.5 * C_1);

    TArray4d rmbp_kl = ijxy("<m b'|r|k l>");
    TArray4d rabp_kl = ijxy("<a b'|r|k l>");

    TArray4 V_ALCJ = VPQ_RS(a, l, c, j);
    TArray4 V_LACJ = VPQ_RS(l, a, c, j);
    TArray4 V_KLCA = VPQ_RS(k, l, c, a);
    TArray4 V_KLMJ = VPQ_RS(k, l, m, j);

    TArray2 Xam_VL2T1;
    Xam_VL2T1("a,m") =  // - 1/2 l^ad_kl R^kl_AB g^AB_mj t^j_d
                      - (R_C1 * V_KLMJ("k,l,m,j") + R_C2 * V_KLMJ("l,k,m,j"))
                        * L2("a,d,k,l") * T1("d,j")

                        // - 1/2 l^cd_kl R^kl_mb' g^ab'_cj t^j_d
                        // - 1/2 l^cd_kl R^kl_ab' g^mb'_cj t^j_d
                      - (  (R_C1 * rmbp_kl("m,b',k,l") + R_C2 * rmbp_kl("m,b',l,k"))
                           * _4("<a b'|g|c j>")
                         + (R_C1 * rmbp_kl("m,b',l,k") + R_C2 * rmbp_kl("m,b',k,l"))
                           * _4("<a b'|g|j c>")
                           //
                         + (R_C1 * rabp_kl("a,b',k,l") + R_C2 * rabp_kl("a,b',l,k"))
                           * _4("<m b'|g|c j>")
                         + (R_C1 * rabp_kl("a,b',l,k") + R_C2 * rabp_kl("a,b',k,l"))
                           * _4("<m b'|g|j c>")
                        ) * L2("c,d,k,l") * T1("d,j")
                       ;
    TArray2 Xai_VL2T1;
    Xai_VL2T1("a,i") =  // + 1/2 l^cd_il R^al_AB g^AB_cj t^j_d
                        (  R_C1 * (  V_ALCJ("a,l,c,j") * L2("c,d,i,l")
                                   + V_LACJ("l,a,c,j") * L2("c,d,l,i"))
                         + R_C2 * (  V_ALCJ("a,l,c,j") * L2("c,d,l,i")
                                   + V_LACJ("l,a,c,j") * L2("c,d,i,l"))
                        ) * T1("d,j")

                        // + 1/2 l^cd_kl R^kl_AB g^AB_ca t^i_d
                      + (R_C1 * V_KLCA("k,l,c,a") + R_C2 * V_KLCA("l,k,c,a"))
                        * L2("c,d,k,l") * T1("d,i")
                       ;

    return XaiAddToXam(Xam_VL2T1, Xai_VL2T1);

  }

  template <typename T>
  typename SingleReference_R12Intermediates<T>::TArray2
  SingleReference_R12Intermediates<T>::Xam_VL2T1_f12b_test(const double C_0, const double C_1,
                                                           const TArray2& T1, const TArray4& L2) {

    const double R_C1 = (0.5 * C_0 + 1.5 * C_1);
    const double R_C2 = (0.5 * C_0 - 1.5 * C_1);

    const char* j = "j";
    const char* k = "k";
    const char* l = "l";
    const char* m = "m";

    const char* a = "a";
    const char* c = "c";

    TArray4d r_klmaP = ijxy("<k l|r|m a'>");
    TArray4d r_abPkl = ijxy("<a b'|r|k l>");

    TArray4 V_alcj = VPq_Rs(a,l,c,j,C_0, C_1);
    TArray4 V_lacj = VPq_Rs(l,a,c,j,C_0, C_1);
    TArray4 V_klca = VPq_Rs(k,l,c,a,C_0, C_1);
    TArray4 V_klmj = VPq_Rs(k,l,m,j,C_0, C_1);
    TArray2 Xam, Xai;
    Xam("a,m") =  //
                - L2("a,d,k,l") * V_klmj("k,l,m,j") * T1("d,j")

                - (  //
                     (R_C1 * r_abPkl("a,b',k,l") + R_C2 * r_abPkl("a,b',l,k"))
                     * _4("<c j|g|m b'>")
                   + (R_C1 * r_abPkl("a,b',l,k") + R_C2 * r_abPkl("a,b',k,l"))
                     * _4("<j c|g|m b'>")
                     //
                   + _4("<a b'|g|c j>")
                     * (R_C1 * r_klmaP("k,l,m,b'") + R_C2 * r_klmaP("l,k,m,b'"))
                   + _4("<a b'|g|j c>")
                     * (R_C1 * r_klmaP("l,k,m,b'") + R_C2 * r_klmaP("k,l,m,b'"))
                  ) * L2("c,d,k,l") * T1("d,j")
                ;
    Xai("a,i") =  //
                  (V_alcj("a,l,c,j") * L2("c,d,i,l") + V_lacj("l,a,c,j") * L2("c,d,l,i"))
                  * T1("d,j")
                  //
                + L2("c,d,k,l") * V_klca("k,l,c,a") * T1("d,i")
                ;

    return XaiAddToXam(Xam, Xai);
  }

  // F12b method
  // Xam contribution of VL2T1T1 part resulted from CCSD F12 coupling
  template <typename T>
  typename SingleReference_R12Intermediates<T>::TArray2
  SingleReference_R12Intermediates<T>::Xam_VL2T1T1_f12b(const double C_0, const double C_1,
                                                        const TArray2& T1, const TArray4& L2) {

    const char* k = "k";
    const char* l = "l";
    const char* j = "j";
    const char* a = "a";

    const double R_C1 = (0.5 * C_0 + 1.5 * C_1);
    const double R_C2 = (0.5 * C_0 - 1.5 * C_1);

    TArray4d rmbp_kl = ijxy("<m b'|r|k l>");
    TArray4d rabp_kl = ijxy("<a b'|r|k l>");

    TArray4 V_AJKL = VPQ_RS(a, j, k, l);
    TArray4 V_KLAJ = VPQ_RS(k, l, a, j);

    // compute Xam from 1/4 L^cd_kl V^cd_k1l1 t^k1_c t^l1_d
    //                 (1/8 L^cd_kl R^kl_AB g^AB_k1l1 t^k1_c t^l1_d)
    TArray2 Xam_VL2T1T1;
    Xam_VL2T1T1("a,m") =  // 1/4 l^cd_kl R^kl_mb' g^ab'_k1l1 t^k1_c t^l1_d
                          // 1/4 l^cd_kl R^kl_ab' g^mb'_k1l1 t^k1_c t^l1_d
                         ( (R_C1 * rmbp_kl("m,b',k,l") + R_C2 * rmbp_kl("m,b',l,k"))
                            * _4("<a b'|g|i j>")
                            //
                          + (R_C1 * rabp_kl("a,b',k,l") + R_C2 * rabp_kl("a,b',l,k"))
                            * _4("<m b'|g|i j>")
                         ) * L2("c,d,k,l") * T1("c,i") * T1("d,j")
                         ;
    TArray2 Xai_VL2T1T1;
    Xai_VL2T1T1("a,i") =  // - 1/4 l^cd_il R^al_AB g^AB_k1l1 t^k1_c t^l1_d
                        - (R_C1 * V_AJKL("a,j,k,l") + R_C2 * V_AJKL("a,j,l,k"))
                          * L2("c,d,i,j") * T1("c,k") * T1("d,l")

                          // - 1/4 l^cd_kl R^kl_AB g^AB_al1 t^i_c t^l1_d
                        - (R_C1 * V_KLAJ("k,l,a,j") + R_C2 * V_KLAJ("l,k,a,j"))
                          * L2("c,d,k,l") * T1("c,i") * T1("d,j")
                        ;

    return XaiAddToXam(Xam_VL2T1T1, Xai_VL2T1T1);
  }

  template <typename T>
  typename SingleReference_R12Intermediates<T>::TArray2
  SingleReference_R12Intermediates<T>::Xam_VL2T1T1_f12b_test(const double C_0, const double C_1,
                                                             const TArray2& T1, const TArray4& L2) {

    const double R_C1 = (0.5 * C_0 + 1.5 * C_1);
    const double R_C2 = (0.5 * C_0 - 1.5 * C_1);

    const char* j = "j";
    const char* k = "k";
    const char* l = "l";
    const char* a = "a";

    TArray4d r_klmaP = ijxy("<k l|r|m a'>");
    TArray4d r_abPkl = ijxy("<a b'|r|k l>");

    TArray4 V_ajkl = VPq_Rs(a,j,k,l,C_0, C_1);
    TArray4 V_klaj = VPq_Rs(k,l,a,j,C_0, C_1);
    TArray2 Xam, Xai;
    Xam("a,m") =  (  //
                     (R_C1 * r_abPkl("a,b',k,l") + R_C2 * r_abPkl("a,b',l,k"))
                     * _4("<i j|g|m b'>")
                     //
                   + _4("<a b'|g|i j>")
                     * (R_C1 * r_klmaP("k,l,m,b'") + R_C2 * r_klmaP("l,k,m,b'"))
                  ) * L2("c,d,k,l") * T1("c,i") * T1("d,j")
                ;
    Xai("a,i") =  //
                - L2("c,d,k,l") * V_klaj("k,l,a,j") * T1("c,i") * T1("d,j")
                  //
                - V_ajkl("a,j,k,l") * T1("c,k") * T1("d,l") * L2("c,d,i,j")
                ;

    return XaiAddToXam(Xam, Xai);
  }

  // frozen-core Xii' contribution resulting from CT2, VT2, and VT1 in F12b
  template <typename T>
  typename SingleReference_R12Intermediates<T>::TArray2
  SingleReference_R12Intermediates<T>::Xiip_CVT_f12b(const double C_0, const double C_1,
                                                     const TArray2& T1, const TArray4& T2,
                                                     const TArray2& L1, const TArray4& L2){

    const double R_C1 = (0.5 * C_0 + 1.5 * C_1);
    const double R_C2 = (0.5 * C_0 - 1.5 * C_1);

    const char* ip = "i'";
    const char* i = "i";
    const char* l = "l";
    const char* c = "c";
    const char* d = "d";

    TArray2 Xiip_CTL;
    Xiip_CTL("i,i'") = // F^a'_b R^i'l_a'c (T2^bc_il + L2^bc_il)
                       _4("<i' l|r|b_F(a') c>")
                       * (  R_C1 * T2("b,c,i,l") + R_C2 * T2("b,c,l,i")
                          + R_C1 * L2("b,c,i,l") + R_C2 * L2("b,c,l,i"))
                      +  _4("<i' l|r|c b_F(a')>")
                       * (  R_C1 * T2("c,b,i,l") + R_C2 * T2("c,b,l,i")
                          + R_C1 * L2("c,b,i,l") + R_C2 * L2("c,b,l,i"))
                      ;

    TArray2 Xiip_VT2L2;
    TArray4 V_ipl_cd = VPq_Rs(ip,l,c,d,C_0, C_1);
    Xiip_VT2L2("i,i'") = // 1/4 R^i'l_A'B' g^A'B'_cd (T^cd_il +L^cd_il)
                         V_ipl_cd("i',l,c,d") * (T2("c,d,i,l") + L2("c,d,i,l"));

    TArray2 Xiip_VT1L1;
    TArray4 V_ipl_ic = VPq_Rs(ip,l,i,c, C_0, C_1);
    TArray4 V_il_ipc = VPq_Rs(i,l,ip,c,C_0, C_1);
    TArray2 V_ipc = VRk_Sk(ip,c,C_0, C_1);
    Xiip_VT1L1("i,i'") = // 1/2 R^i'l_A'B' g^A'B'_ic (T^c_l + L^c_l)
                         // 1/2 R^il_A'B' g^A'B'_i'c (T^c_l + L^c_l)
                         (V_ipl_ic("i',l,i,c") + V_il_ipc("i,l,i',c"))
                         * (T1("c,l") + L1("c,l"))

                          // 1/2 R^i'k_A'B' g^A'B'_ck (T^c_i + L^c_i)
                        + V_ipc("i',c") * (T1("c,i") + L1("c,i"))
                        ;

    TArray2 Xiip_C1;
    Xiip_C1("i,i'") = Xiip_CTL("i,i'") + Xiip_VT2L2("i,i'") + Xiip_VT1L1("i,i'");
    return Xiip_C1;
  }

  // frozen-core Xii' contribution resulting from VT1T1 in F12b
  template <typename T>
  typename SingleReference_R12Intermediates<T>::TArray2
  SingleReference_R12Intermediates<T>::Xiip_VT1T1_f12b(const double C_0, const double C_1,
                                                       const TArray2& T1, const TArray4& T2,
                                                       const TArray4& L2) {
    const double R_C1 = (0.5 * C_0 + 1.5 * C_1);
    const double R_C2 = (0.5 * C_0 - 1.5 * C_1);

    const char* ip = "i'";
    const char* i = "i";
    const char* j = "j";
    const char* k = "k";
    const char* l = "l";
    const char* c = "c";
    const char* d = "d";

    TArray4 V_ipl_cd = VPq_Rs(ip,l,c,d,C_0, C_1);
    TArray2 Xiip_VT1T1;
    Xiip_VT1T1("i,i'") = // 1/2 R^i'l_A'B' g^A'B'_cd T^i_c t^l_d
                         V_ipl_cd("i',l,c,d") * T1("c,i") * T1("d,l");


    TArray4 V_IPL_CJ = VPQ_RS(ip, l, c, j);
    TArray4 V_LIP_CJ = VPQ_RS(l, ip, c, j);
    TArray4 V_KL_CIP = VPQ_RS(k, l, c, ip);
    TArray2 Xai_VL2T1;
    Xai_VL2T1("i,i'") =  // - 1/2 l^cd_il R^i'l_AB g^AB_cj t^j_d
                         (  R_C1 * (  V_IPL_CJ("i',l,c,j") * L2("c,d,i,l")
                                    + V_LIP_CJ("l,i',c,j") * L2("c,d,l,i"))
                          + R_C2 * (  V_IPL_CJ("i',l,c,j") * L2("c,d,l,i")
                                    + V_LIP_CJ("l,i',c,j") * L2("c,d,i,l"))
                         ) * T1("d,j")

                         // - 1/2 l^cd_kl R^kl_AB g^AB_ci' t^i_d
                       + (R_C1 * V_KL_CIP("k,l,c,i'") + R_C2 * V_KL_CIP("l,k,c,i'"))
                         * L2("c,d,k,l") * T1("d,i")
                        ;

    TArray4 V_IPJ_KL = VPQ_RS(ip, j, k, l);
    TArray4 V_KL_IPJ = VPQ_RS(k, l, ip, j);
    TArray2 Xai_VL2T1T1;
    Xai_VL2T1T1("i,i'") =  //   1/4 l^cd_il R^i'l_AB g^AB_k1l1 t^k1_c t^l1_d
                           (R_C1 * V_IPJ_KL("i',j,k,l") + R_C2 * V_IPJ_KL("i',j,l,k"))
                           * L2("c,d,i,j") * T1("c,k") * T1("d,l")

                           // + 1/4 l^cd_kl R^kl_AB g^AB_i'l1 t^i_c t^l1_d
                         + (R_C1 * V_KL_IPJ("k,l,i',j") + R_C2 * V_KL_IPJ("l,k,i',j"))
                           * L2("c,d,k,l") * T1("c,i") * T1("d,j")
                          ;

    TArray2 Xiip_C2;
    Xiip_C2("i,i'") = Xiip_VT1T1("i,i'") + Xai_VL2T1("i,i'") + Xai_VL2T1T1("i,i'");
    return Xiip_C2;
  }

  template <typename T>
  void SingleReference_R12Intermediates<T>::compute_multipole_F12b_coupling() {

    bool compute_dipole = true;
    bool compute_quadrupole = true;

    double conv_target = 1e-10;

    // singlet and triplet coefficients for F12 and coupling terms
    const double C_0 = 1.0 / 2.0;
    const double C_1 = 1.0 / 4.0;
    const double RC1 = 0.5 * (C_0 + C_1);
    const double RC2 = 0.5 * (C_0 - C_1);
    // compute coefficients needed in the F12 and coupling calculations
    const double R_C1 = (0.5 * C_0 + 1.5 * C_1);
    const double R_C2 = (0.5 * C_0 - 1.5 * C_1);

    const char* i = "i";
    const char* j = "j";
    const char* k = "k";
    const char* l = "l";
    const char* m = "m";

    const char* a = "a";
    const char* b = "b";
    const char* c = "c";
    const char* d = "d";

    // compute integrals needed for orbital relaxation
    // i.e. solve Abnam Dbn = Xam (close-shell formula)
    TArray4d g_mnab = ijxy("<m n|g|a b>");
    TArray4 A_bnam;
    A_bnam("b,n,a,m") = - _4("<b n|g|a m>") - g_mnab("m,n,b,a") + 4.0 * g_mnab("n,m,b,a")
                        + _2("<b|F|a>") * _2("<m|I|n>") - _2("<a|I|b>") * _2("<m|F|n>");

    // Make preconditioner: Delta_am = 1 / (<a|F|a> - <m|F|m>) for
    // solving k_bn A_bnam = X_am
    TArray2 mFmn, mFab;
    mFmn("m,n") = - _2("<m|F|n>");
    mFab("a,b") = - _2("<a|F|b>");
    typedef detail::diag_precond2<double> pceval_type; //!< evaluator of preconditioner
    pceval_type Delta_am_gen(TA::array_to_eigen(mFab), TA::array_to_eigen(mFmn));

    typedef TA::Array<T, 2, LazyTensor<T, 2, pceval_type > > TArray2d;
    TArray2 Iam = xy("<a|I|m>");
    TArray2d Delta_am(Iam.get_world(), Iam.trange());

    // construct local tiles
    for(auto t = Delta_am.trange().tiles().begin();
        t != Delta_am.trange().tiles().end(); ++t)
      if (Delta_am.is_local(*t)) {
        std::array<std::size_t, 2> index;
        std::copy(t->begin(), t->end(), index.begin());
        madness::Future < typename TArray2d::value_type >
          tile((LazyTensor<T, 2, pceval_type >(&Delta_am, index, &Delta_am_gen)
              ));

        // Insert the tile into the array
        Delta_am.set(*t, tile);
      }
    TArray2 preconditioner;
    preconditioner("a,m") = Delta_am("a,m");

    detail::Orbital_relaxation_Abjai<double> Orbital_relaxation_Abnam(A_bnam);
    TA::ConjugateGradientSolver<TiledArray::Array<T,2>,
                                detail::Orbital_relaxation_Abjai<double> > cg_solver2;


    // determine if it is frozen-core
    TArray2 Fij = xy("<i|F|j>");
    const std::size_t nocc = mFmn.trange().elements().extent()[0];
    const std::size_t naocc = Fij.trange().elements().extent()[0];

    // compute terms needed for frozen-core
    TArray2 Delta_ijp_F;
    TArray4 A_ijpam;
    if (nocc != naocc) {
      // compute 1 / (Fi'i' - Fjj)
      TArray2 Fipjp = xy("<i'|F|j'>");
      pceval_type Delta_ijp_gen(TA::array_to_eigen(Fij), TA::array_to_eigen(Fipjp));

      TArray2 Fijp = xy("<i|F|j'>");
      TArray2d Delta_ijp(Fijp.get_world(), Fijp.trange());
      // construct local tiles
      for(auto t = Delta_ijp.trange().tiles().begin();
          t != Delta_ijp.trange().tiles().end(); ++t)
        if (Delta_ijp.is_local(*t)) {
          std::array<std::size_t, 2> index;
          std::copy(t->begin(), t->end(), index.begin());
          madness::Future < typename TArray2d::value_type >
            tile((LazyTensor<T, 2, pceval_type >(&Delta_ijp, index, &Delta_ijp_gen)
                ));

          // Insert the tile into the array
          Delta_ijp.set(*t, tile);
        }
      Delta_ijp_F("i,j'") = Delta_ijp("i,j'");

      A_ijpam("i,j',a,m") = 4.0 * _4("<i a|g|j' m>") - _4("<i a|g|m j'>") - _4("<i m|g|a j'>");
    }

    TArray2 mu_z_ij, mu_z_ijp, mu_z_ab, mu_z_ai, mu_z_am;
    if (compute_dipole) {
      mu_z_ij = xy("<i|mu_z|j>");
      mu_z_ab = xy("<a|mu_z|b>");
      mu_z_ai = xy("<a|mu_z|i>");
      mu_z_am = xy("<a|mu_z|m>");

      if (nocc != naocc)
        mu_z_ijp = xy("<i|mu_z|j'>");
    }

    // quadrupole integrals
    TArray2 Qxx_ij, Qyy_ij, Qzz_ij, Qxx_ab, Qyy_ab, Qzz_ab,
            Qxx_ai, Qyy_ai, Qzz_ai, Qxx_am, Qyy_am, Qzz_am,
            Qxx_apb, Qyy_apb, Qzz_apb,
            // integrals for frozen-core
            Qxx_ijp, Qyy_ijp, Qzz_ijp;

    if (compute_quadrupole){

      TArray2 q_xx_ij = xy("<i|q_xx|j>");
      TArray2 q_yy_ij = xy("<i|q_yy|j>");
      TArray2 q_zz_ij = xy("<i|q_zz|j>");
      Qxx_ij("i,j") = q_xx_ij("i,j") - (q_zz_ij("i,j") + q_yy_ij("i,j")) * 0.5;
      Qyy_ij("i,j") = q_yy_ij("i,j") - (q_zz_ij("i,j") + q_xx_ij("i,j")) * 0.5;
      Qzz_ij("i,j") = q_zz_ij("i,j") - (q_xx_ij("i,j") + q_yy_ij("i,j")) * 0.5;

      TArray2 q_xx_ab = xy("<a|q_xx|b>");
      TArray2 q_yy_ab = xy("<a|q_yy|b>");
      TArray2 q_zz_ab = xy("<a|q_zz|b>");
      Qxx_ab("a,b") = q_xx_ab("a,b") - (q_zz_ab("a,b") + q_yy_ab("a,b")) * 0.5;
      Qyy_ab("a,b") = q_yy_ab("a,b") - (q_zz_ab("a,b") + q_xx_ab("a,b")) * 0.5;
      Qzz_ab("a,b") = q_zz_ab("a,b") - (q_xx_ab("a,b") + q_yy_ab("a,b")) * 0.5;

      TArray2 q_xx_ai = xy("<a|q_xx|i>");
      TArray2 q_yy_ai = xy("<a|q_yy|i>");
      TArray2 q_zz_ai = xy("<a|q_zz|i>");
      Qxx_ai("a,i") = q_xx_ai("a,i") - (q_zz_ai("a,i") + q_yy_ai("a,i")) * 0.5;
      Qyy_ai("a,i") = q_yy_ai("a,i") - (q_zz_ai("a,i") + q_xx_ai("a,i")) * 0.5;
      Qzz_ai("a,i") = q_zz_ai("a,i") - (q_xx_ai("a,i") + q_yy_ai("a,i")) * 0.5;

      TArray2 q_xx_am = xy("<a|q_xx|m>");
      TArray2 q_yy_am = xy("<a|q_yy|m>");
      TArray2 q_zz_am = xy("<a|q_zz|m>");
      Qxx_am("a,m") = q_xx_am("a,m") - (q_zz_am("a,m") + q_yy_am("a,m")) * 0.5;
      Qyy_am("a,m") = q_yy_am("a,m") - (q_zz_am("a,m") + q_xx_am("a,m")) * 0.5;
      Qzz_am("a,m") = q_zz_am("a,m") - (q_xx_am("a,m") + q_yy_am("a,m")) * 0.5;

      TArray2 q_xx_apb = xy("<a'|q_xx|b>");
      TArray2 q_yy_apb = xy("<a'|q_yy|b>");
      TArray2 q_zz_apb = xy("<a'|q_zz|b>");
      Qxx_apb("a',b") = q_xx_apb("a',b") - (q_zz_apb("a',b") + q_yy_apb("a',b")) * 0.5;
      Qyy_apb("a',b") = q_yy_apb("a',b") - (q_zz_apb("a',b") + q_xx_apb("a',b")) * 0.5 ;
      Qzz_apb("a',b") = q_zz_apb("a',b") - (q_xx_apb("a',b") + q_yy_apb("a',b")) * 0.5;

      if (nocc != naocc) {
        TArray2 q_xx_ijp = xy("<i|q_xx|j'>");
        TArray2 q_yy_ijp = xy("<i|q_yy|j'>");
        TArray2 q_zz_ijp = xy("<i|q_zz|j'>");
        Qxx_ijp("i,j'") = q_xx_ijp("i,j'") - (q_zz_ijp("i,j'") + q_yy_ijp("i,j'")) * 0.5;
        Qyy_ijp("i,j'") = q_yy_ijp("i,j'") - (q_zz_ijp("i,j'") + q_xx_ijp("i,j'")) * 0.5;
        Qzz_ijp("i,j'") = q_zz_ijp("i,j'") - (q_xx_ijp("i,j'") + q_yy_ijp("i,j'")) * 0.5;
      }
    }

    double muz_ccsd = 0.0, muz_ccsdor = 0.0;
    double q_xx_ccsd = 0.0, q_yy_ccsd = 0.0, q_zz_ccsd = 0.0;
    double q_xx_ccsdor = 0.0, q_yy_ccsdor = 0.0, q_zz_ccsdor = 0.0;
#if 1 // CCSD
    ExEnv::out0() << indent << "Compute CCSD T amplitudes " << std::endl;
    TArray2 T1_ccsd;
    TArray4 T2_ccsd;
    compute_T_ccsd(T1_ccsd, T2_ccsd, "CCSD");

    ExEnv::out0() << indent << "Compute CCSD L amplitudes " << std::endl;
    TArray2 L1_ccsd;
    TArray4 L2_ccsd;
    compute_lambda_ccsd(T1_ccsd, T2_ccsd, L1_ccsd, L2_ccsd, "CCSD");

    // compute CCSD density from amplitudes
    TArray2 Dij_ccsd, Dab_ccsd, Dia_ccsd, Dai_ccsd;
    ExEnv::out0() << std::endl << indent
                  << "Compute CCSD density from amplitudes" << std::endl;
    compute_ccsd_1rdm_amp(T1_ccsd, T2_ccsd, L1_ccsd, L2_ccsd,
                          Dij_ccsd, Dab_ccsd, Dia_ccsd, Dai_ccsd);

    ExEnv::out0() << std::endl << indent << "Compute CCSD Xam and Xai" << std::endl;
    TArray2 Xam_ccsd_nfzc, Xiip_ccsd;
    compute_Xam_ccsd(T1_ccsd, T2_ccsd, L1_ccsd, L2_ccsd,
                     Xam_ccsd_nfzc, Xiip_ccsd);

    TArray2 Diip_or_ccsd, Xam_ccsd;
    if (nocc != naocc) {
      ExEnv::out0() << std::endl << indent
                      << "Include frozen-core CCSD contributions" << std::endl;
      Diip_or_ccsd("i,i'") = Xiip_ccsd("i,i'") * Delta_ijp_F("i,i'");

      TArray2 Xam_Diip_ccsd;
      Xam_Diip_ccsd("a,m") = A_ijpam("i,i',a,m") * Diip_or_ccsd("i,i'");
      Xam_ccsd("a,m") = Xam_ccsd_nfzc("a,m") - Xam_Diip_ccsd("a,m");
    } else {
      Xam_ccsd("a,m") = Xam_ccsd_nfzc("a,m");
    }

    // solve the Z-vector equation for CCSD
    TArray2 Dbn_ccsd(Xam_ccsd.get_world(), Xam_ccsd.trange());
    auto resnorm_ccsd = cg_solver2(Orbital_relaxation_Abnam,
                                   Xam_ccsd,
                                   Dbn_ccsd,
                                   preconditioner,
                                   conv_target);

    if (compute_dipole) {
      muz_ccsd =  dot(mu_z_ij("i,j"), Dij_ccsd("i,j"))
                + dot(mu_z_ab("a,b"), Dab_ccsd("a,b"))
                + dot(mu_z_ai("a,i"), Dai_ccsd("a,i"))
                + dot(mu_z_ai("a,i"), Dia_ccsd("i,a"))
                ;

      muz_ccsdor = dot(mu_z_am("a,m"), Dbn_ccsd("a,m"));
      if (nocc != naocc) {
        muz_ccsdor += dot(mu_z_ijp("i,i'"), Diip_or_ccsd("i,i'"));
      }

      std::cout << std::endl << indent
                << "muz_CCSD = " << scprintf("%15.12f", - muz_ccsd * 2.0)
                << std::endl << indent
                << "muz_CCSDor = " << scprintf("%15.12f", - muz_ccsdor * 2.0)
                << std::endl;
    }

    if (compute_quadrupole) {
      q_xx_ccsd =  dot(Qxx_ij("i,j"), Dij_ccsd("i,j"))
                 + dot(Qxx_ab("a,b"), Dab_ccsd("a,b"))
                 + dot(Qxx_ai("a,i"), Dai_ccsd("a,i"))
                 + dot(Qxx_ai("a,i"), Dia_ccsd("i,a"))
                 ;
      q_yy_ccsd =  dot(Qyy_ij("i,j"), Dij_ccsd("i,j"))
                 + dot(Qyy_ab("a,b"), Dab_ccsd("a,b"))
                 + dot(Qyy_ai("a,i"), Dai_ccsd("a,i"))
                 + dot(Qyy_ai("a,i"), Dia_ccsd("i,a"))
                 ;
      q_zz_ccsd =  dot(Qzz_ij("i,j"), Dij_ccsd("i,j"))
                 + dot(Qzz_ab("a,b"), Dab_ccsd("a,b"))
                 + dot(Qzz_ai("a,i"), Dai_ccsd("a,i"))
                 + dot(Qzz_ai("a,i"), Dia_ccsd("i,a"))
                 ;

      q_xx_ccsdor = dot(Qxx_am("a,m"), Dbn_ccsd("a,m"));
      q_yy_ccsdor = dot(Qyy_am("a,m"), Dbn_ccsd("a,m"));
      q_zz_ccsdor = dot(Qzz_am("a,m"), Dbn_ccsd("a,m"));

      if (nocc != naocc) {
        q_xx_ccsdor += dot(Qxx_ijp("i,j'"), Diip_or_ccsd("i,j'"));
        q_yy_ccsdor += dot(Qyy_ijp("i,j'"), Diip_or_ccsd("i,j'"));
        q_zz_ccsdor += dot(Qzz_ijp("i,j'"), Diip_or_ccsd("i,j'"));
      }

      std::cout << std::endl << indent
                << "qxx_ccsd = " << scprintf("%12.10f", - q_xx_ccsd * 2.0)
                << "  qyy_ccsd = " << scprintf("%12.10f", - q_yy_ccsd * 2.0)
                << "  qzz_ccsd = " << scprintf("%12.10f", - q_zz_ccsd * 2.0)
                << std::endl << indent
                << "qxx_ccsdor = " << scprintf("%12.10f", - q_xx_ccsdor * 2.0)
                << "  qyy_ccsdor = " << scprintf("%12.10f", - q_yy_ccsdor * 2.0)
                << "  qzz_ccsdor = " << scprintf("%12.10f", - q_zz_ccsdor * 2.0)
                << std::endl;
    }
#endif

    // resolve CCSD T amplitudes in the framework of CCSD-F12b
    ExEnv::out0() << indent << "Compute CCSD-F12b T amplitudes " << std::endl;
    TArray2 T1_f12b;
    TArray4 T2_f12b;
    compute_T_ccsd(T1_f12b, T2_f12b, "F12b");

    // resolve CCSD lambda amplitudes in the framework of CCSD-F12b
    ExEnv::out0() << std::endl << indent << "Compute CCSD-F12b L amplitudes " << std::endl;
    TArray2 L1_f12b;
    TArray4 L2_f12b;
    compute_lambda_ccsd(T1_f12b, T2_f12b, L1_f12b, L2_f12b, "F12b");

    // compute energy from coupling in F12b
#if 1
    // CT2 energy
    // F^a'_a R^ij_a'b + F^a'_b R^ij_aa'
    TArray4 Cijab;
    TArray4d r_ijaFb = ijxy("<i j|r|a_F(c') b>");
    Cijab("i,j,a,b") =  R_C1 * r_ijaFb("i,j,a,b") + R_C2 * r_ijaFb("j,i,a,b")
                      + R_C1 * r_ijaFb("j,i,b,a") + R_C2 * r_ijaFb("i,j,b,a");
    TArray2 CT2;
    CT2("i,j") = Cijab("i,k,a,b") * T2_f12b("a,b,j,k");
    double E_CT2 = dot(CT2("i,j"), _2("<j|I|i>"));
    std::cout << std::endl
              << indent << "E_CT2 =  " << scprintf("%15.12f", E_CT2)
              << std::endl;

    // VT1 + VT2 energy
    TArray2 V_ia;
    V_ia = VRk_Sk(i,a, C_0, C_1);

    TArray4 Vijab = VPq_Rs(i, j, a, b, C_0, C_1);
    TArray2 VT2;
    VT2("i,j") = Vijab("i,k,a,b") * T2_f12b("a,b,j,k");

    double E_VT =  dot(V_ia("i,a"), T1_f12b("a,i")) * 2.0
                 + dot(VT2("i,j"), _2("<j|I|i>"));
    std::cout << indent << "E_VT =  " << scprintf("%15.12f", E_VT)
              << std::endl;

    // VT1T1 energy
    double E_VT1T1 = dot(Vijab("i,j,a,b"), T1_f12b("a,i") * T1_f12b("b,j"));
    std::cout << indent << "E_VT1T1 =  " << scprintf("%15.12f", E_VT1T1)
              << std::endl;

    std::cout << indent << "E (Coupling in F12b) =  "
              << scprintf("%15.12f", (E_CT2+E_VT+E_VT1T1))
              << std::endl;
#endif

    // CC F12b coupling contribution to Xam
    TArray2 mu_z_apb;
#if 1
    // compute CCSD-F12b density from amplitudes
    TArray2 Dij_f12b, Dab_f12b, Dia_f12b, Dai_f12b;
    ExEnv::out0() << std::endl << indent
                  << "Compute CCSD density from amplitudes in CCSD-F12b" << std::endl;
    compute_ccsd_1rdm_amp(T1_f12b, T2_f12b, L1_f12b, L2_f12b,
                          Dij_f12b, Dab_f12b, Dia_f12b, Dai_f12b);

    ExEnv::out0() << std::endl << indent
                  << "Compute CCSD Xam and Xai in CCSD-F12b" << std::endl;
    TArray2 Xam_f12b_nfzc, Xiip_f12b;
    compute_Xam_ccsd(T1_f12b, T2_f12b, L1_f12b, L2_f12b, Xam_f12b_nfzc, Xiip_f12b);

    TArray2 Diip_f12b, Xam_f12b;
    if (nocc != naocc) {
      Diip_f12b("i,i'") = Xiip_f12b("i,i'") * Delta_ijp_F("i,i'");
      Xam_f12b("a,m") = Xam_f12b_nfzc("a,m") - A_ijpam("i,i',a,m") * Diip_f12b("i,i'");
    } else {
      Xam_f12b("a,m") = Xam_f12b_nfzc("a,m");
    }

    TArray2 Dbn_f12b(Xam_f12b.get_world(), Xam_f12b.trange());
    auto resnorm_f12b = cg_solver2(Orbital_relaxation_Abnam,
                                   Xam_f12b,
                                   Dbn_f12b,
                                   preconditioner,
                                   conv_target);


    ExEnv::out0() << std::endl << indent
                  << "Compute CCSD F12b coupling contributions" << std::endl;
    // CC CT2
    TArray2 RT2_aPb, RL2_aPb;
    // 1/2 R^kl_a'c T2^bc_kl
    RT2_aPb("a',b") = _4("<a' c|r|k l>")
                      * (R_C1 * T2_f12b("b,c,k,l") + R_C2 * T2_f12b("c,b,k,l"));
    // 1/2 R^kl_a'c L2^bc_kl
    RL2_aPb("a',b") = _4("<a' c|r|k l>")
                      * (R_C1 * L2_f12b("b,c,k,l") + R_C2 * L2_f12b("c,b,k,l"));

    TArray2 Xam_CT2L2 = Xam_CT2L2_f12b(C_0, C_1, T2_f12b, RT2_aPb, L2_f12b, RL2_aPb);

    // VT1 & VT2 coupling contribution to F12 Xam
    TArray2 Xam_VTL = Xam_VTL_f12b(C_0, C_1, T1_f12b, T2_f12b, L1_f12b, L2_f12b);

    // VT1T1
    TArray2 X_VT1T1 = Xam_VT1T1_f12b(C_0, C_1, T1_f12b);
//    // *** test
//    TArray4 tauT1_ab;
//    tauT1_ab("a,b,i,j") = T1_f12b("a,i") * T1_f12b("b,j");
//    TArray2 X_VT1T1 = Xam_VT1T1_f12b_test(C_0, C_1, tauT1_ab);

    // VL2T1
    TArray2 X_VL2T1 = Xam_VL2T1_f12b(C_0, C_1, T1_f12b, L2_f12b);
//    // *** test
//    TArray2 X_VL2T1 = Xam_VL2T1_f12b_test(C_0, C_1, T1_f12b, L2_f12b);

    // VL2T1T1
    TArray2 X_VL2T1T1 = Xam_VL2T1T1_f12b(C_0, C_1, T1_f12b, L2_f12b);
//    // *** test
//    TArray2 X_VL2T1T1 = Xam_VL2T1T1_f12b_test(C_0, C_1, T1_f12b, L2_f12b);

    TArray2 Xam_C_nfzc;
    Xam_C_nfzc("a,m") =  Xam_CT2L2("a,m")
                       + Xam_VTL("a,m")
                       + X_VT1T1("a,m")
                       + X_VL2T1("a,m")
                       + X_VL2T1T1("a,m")
                       ;

    // frozen-core contribution
    TArray2 Diip_C, Xam_C;
    if (nocc != naocc) {
      ExEnv::out0() << std::endl << indent
                      << "Include frozen-core contributions for CCSD F12 coupling" << std::endl;
      TArray2 Xiip_couling1 = Xiip_CVT_f12b(C_0, C_1, T1_f12b, T2_f12b, L1_f12b, L2_f12b);
      TArray2 Xiip_couling2 = Xiip_VT1T1_f12b(C_0, C_1, T1_f12b, T2_f12b, L2_f12b);
      Diip_C("i,i'") = - (Xiip_couling1("i,i'") + Xiip_couling2("i,i'"))
                         * Delta_ijp_F("i,i'");

      TArray2 Xam_Diip_C;
      Xam_Diip_C("a,m") = - A_ijpam("i,j',a,m") * Diip_C("i,j'");

      Xam_C("a,m") = Xam_C_nfzc("a,m") + Xam_Diip_C("a,m");
    } else {
      Xam_C("a,m") = Xam_C_nfzc("a,m");
    }

    TArray2 Dbn_C(Xam_C.get_world(), Xam_C.trange());
    auto resnorm_C = cg_solver2(Orbital_relaxation_Abnam,
                                Xam_C,
                                Dbn_C,
                                preconditioner,
                                conv_target);

//    TArray2 Dbn_CT2L2(Xam_CT2L2.get_world(), Xam_CT2L2.trange());
//    auto resnorm_CT2L2 = cg_solver2(Orbital_relaxation_Abnam,
//                                    Xam_CT2L2,
//                                    Dbn_CT2L2,
//                                    preconditioner,
//                                    conv_target);
//
//    TArray2 Dbn_VTL(Xam_VTL.get_world(), Xam_VTL.trange());
//    auto resnorm_VTL = cg_solver2(Orbital_relaxation_Abnam,
//                                  Xam_VTL,
//                                  Dbn_VTL,
//                                  preconditioner,
//                                  conv_target);

    if (compute_dipole) {
      mu_z_apb = xy("<a'|mu_z|b>");
      double muz_f12b =  dot(mu_z_ij("i,j"), Dij_f12b("i,j"))
                       + dot(mu_z_ab("a,b"), Dab_f12b("a,b"))
                       + dot(mu_z_ai("a,i"), Dai_f12b("a,i"))
                       + dot(mu_z_ai("a,i"), Dia_f12b("i,a"))
                       // F12 coupling density contribution
                       + dot(mu_z_apb("a',b"), RT2_aPb("a',b"))
                       + dot(mu_z_apb("a',b"), RL2_aPb("a',b"))
                       ;

      double muz_f12bor =  dot(mu_z_am("a,m"), Dbn_f12b("a,m"))
                         + dot(mu_z_am("a,m"), Dbn_C("a,m"))
                         ;

      if (nocc != naocc) {
        muz_f12b +=  dot(mu_z_ijp("i,j'"), Diip_f12b("i,j'"))
                   + dot(mu_z_ijp("i,j'"), Diip_C("i,j'"));
      }

//      std::cout << std::endl << indent
//                << "muz_ccsd (F12b) = " << scprintf("%15.12f", - muz_f12b * 2.0)
//                << std::endl << indent
//                << "muz_ccsdor (F12b) = " << scprintf("%15.12f", - muz_f12bor * 2.0)
//                << std::endl;

      double muz_Cf12b =  muz_f12b - muz_ccsd;
      double muz_Cf12bor =  muz_f12bor - muz_ccsdor;

      std::cout << std::endl << indent
                << "mu_z (F12b C) = " << scprintf("%15.12f", - muz_Cf12b * 2.0)
                << std::endl << indent
                << "mu_z (F12b C or) = " << scprintf("%15.12f", - muz_Cf12bor * 2.0)
                << std::endl << std::endl;
    }

    if (compute_quadrupole) {
      double q_xx_f12b =  dot(Qxx_ij("i,j"), Dij_f12b("i,j"))
                        + dot(Qxx_ab("a,b"), Dab_f12b("a,b"))
                        + dot(Qxx_ai("a,i"), Dai_f12b("a,i"))
                        + dot(Qxx_ai("a,i"), Dia_f12b("i,a"))
                        + dot(Qxx_apb("a',b"), RT2_aPb("a',b"))
                        + dot(Qxx_apb("a',b"), RL2_aPb("a',b"))
                        ;
      double q_yy_f12b =  dot(Qyy_ij("i,j"), Dij_f12b("i,j"))
                        + dot(Qyy_ab("a,b"), Dab_f12b("a,b"))
                        + dot(Qyy_ai("a,i"), Dai_f12b("a,i"))
                        + dot(Qyy_ai("a,i"), Dia_f12b("i,a"))
                        + dot(Qyy_apb("a',b"), RT2_aPb("a',b"))
                        + dot(Qyy_apb("a',b"), RL2_aPb("a',b"))
                        ;
      double q_zz_f12b =  dot(Qzz_ij("i,j"), Dij_f12b("i,j"))
                        + dot(Qzz_ab("a,b"), Dab_f12b("a,b"))
                        + dot(Qzz_ai("a,i"), Dai_f12b("a,i"))
                        + dot(Qzz_ai("a,i"), Dia_f12b("i,a"))
                        + dot(Qzz_apb("a',b"), RT2_aPb("a',b"))
                        + dot(Qzz_apb("a',b"), RL2_aPb("a',b"))
                        ;

      double q_xx_f12bor =  dot(Qxx_am("a,m"), Dbn_f12b("a,m"))
                          + dot(Qxx_am("a,m"), Dbn_C("a,m"));
      double q_yy_f12bor =  dot(Qyy_am("a,m"), Dbn_f12b("a,m"))
                          + dot(Qyy_am("a,m"), Dbn_C("a,m"));
      double q_zz_f12bor =  dot(Qzz_am("a,m"), Dbn_f12b("a,m"))
                          + dot(Qzz_am("a,m"), Dbn_C("a,m"));

      if (nocc != naocc) {
        q_xx_f12bor +=  dot(Qxx_ijp("i,j'"), Diip_f12b("i,j'"))
                      + dot(Qxx_ijp("i,j'"), Diip_C("i,j'"));
        q_yy_f12bor +=  dot(Qyy_ijp("i,j'"), Diip_f12b("i,j'"))
                      + dot(Qyy_ijp("i,j'"), Diip_C("i,j'"));
        q_zz_f12bor +=  dot(Qzz_ijp("i,j'"), Diip_f12b("i,j'"))
                      + dot(Qzz_ijp("i,j'"), Diip_C("i,j'"));
      }

//      std::cout << std::endl << indent
//                << "q_xx (F12b) = "   << scprintf("%12.10f", - q_xx_f12b * 2.0)
//                << "  q_yy (F12b) = " << scprintf("%12.10f", - q_yy_f12b * 2.0)
//                << "  q_zz (F12b) = " << scprintf("%12.10f", - q_zz_f12b * 2.0)
//                << std::endl << indent
//                << "q_xx (F12b or) = "   << scprintf("%12.10f", - q_xx_f12bor * 2.0)
//                << "  q_yy (F12b or) = " << scprintf("%12.10f", - q_yy_f12bor * 2.0)
//                << "  q_zz (F12b or) = " << scprintf("%12.10f", - q_zz_f12bor * 2.0)
//                << std::endl << std::endl;

      double q_xx_f12b_C = q_xx_f12b - q_xx_ccsd;
      double q_yy_f12b_C = q_yy_f12b - q_yy_ccsd;
      double q_zz_f12b_C = q_zz_f12b - q_zz_ccsd;
      double q_xx_f12bor_C = q_xx_f12bor - q_xx_ccsdor;
      double q_yy_f12bor_C = q_yy_f12bor - q_yy_ccsdor;
      double q_zz_f12bor_C = q_zz_f12bor - q_zz_ccsdor;
      std::cout << std::endl << indent
                << "q_xx (F12b C) = "   << scprintf("%12.10f", - q_xx_f12b_C * 2.0)
                << "  q_yy (F12b C) = " << scprintf("%12.10f", - q_yy_f12b_C * 2.0)
                << "  q_zz (F12b C) = " << scprintf("%12.10f", - q_zz_f12b_C * 2.0)
                << std::endl << indent
                << "q_xx (F12b C or) = "   << scprintf("%12.10f", - q_xx_f12bor_C * 2.0)
                << "  q_yy (F12b C or) = " << scprintf("%12.10f", - q_yy_f12bor_C * 2.0)
                << "  q_zz (F12b C or) = " << scprintf("%12.10f", - q_zz_f12bor_C * 2.0)
                << std::endl << std::endl;
    }
#endif

  }

  template <typename T>
  typename SingleReference_R12Intermediates<T>::TArray2
  SingleReference_R12Intermediates<T>::rdm1() {

#define ENABLE_SRR12_RDM1 1

#if ENABLE_SRR12_RDM1
//    if (0) {
//      {
//      typedef TiledArray::Array<T,2> Array;
//      Array FiA = _2("<i|F|A'>");
//      Array FAi = _2("<A'|F|i>");
//      std::cout << "<i|F|A'>:" << std::endl << FiA << std::endl;
//      std::cout << "<A'|F|i>:" << std::endl << FAi << std::endl;
//      Array FiA_2(FiA.get_world(), FiA.trange());
//      FiA_2("i,A'") = FAi("A',i");
//      std::cout << "<i|F|A'>=Perm(<A'|F|i>):" << std::endl << FiA_2 << std::endl;
//      }
//      {
//      typedef TiledArray::Array<T,4> Array;
//      Array g_ij_ab = _4("<i j|g|a b>");
//      Array g_ab_ij = _4("<a b|g|i j>");
//      std::cout << "<i j|g|a b>:" << std::endl << g_ij_ab << std::endl;
//      std::cout << "<a b|g|i j>:" << std::endl << g_ab_ij << std::endl;
//      Array g_ij_ab_2(g_ij_ab.get_world(), g_ij_ab.trange());
//      g_ij_ab_2("i,j,a,b") = g_ab_ij("a,b,i,j");
//      std::cout << "<i j|g|a b>=Perm(<a b|g|i j>):" << std::endl << g_ij_ab_2 << std::endl;
//      Array should_be_zero = g_ij_ab("i,j,a,b") - g_ab_ij("a,b,i,j");
//      std::cout << "<i j|g|a b> - Perm(<a b|g|i j>):" << std::endl << should_be_zero << std::endl;
//      const double max_nonzero = norminf(should_be_zero("i,j,a,b"));
//      std::cout << "|| <i j|g|a b> - Perm(<a b|g|i j>) ||_\infty = " << max_nonzero << std::endl;
//      }
//      {
//      typedef TiledArray::Array<T,2> Array;
//      Array mu_z_ij = _2("<i|mu_z|j>");
//      Array gamma_ij = _2("<i|gamma|j>");
//      const double mu_z_e = dot(mu_z_ij("i,j"), gamma_ij("i,j"));
//      double mu_z_n = 0.0;
//      Ref<Molecule> mol = r12world_->basis()->molecule();
//      for(int a=0; a<mol->natom(); ++a) {
//        mu_z_n += mol->Z(a) * mol->r(a, 2);
//      }
//      std::cout << "mu_z = " << -mu_z_e+mu_z_n << std::endl;
//      }
//    }

    // can only ask for T1 with i in bra!
    // since we computed T1 CABS, they are expressed in terms of all virtuals = A'
    // if you turn off vir-CABS coupling, use a' (i.e. CABS only)
//    TArray2 T1iA = _2("<i|T1|A'>");
//    //t1_cabs_.print("T1(RefSCMatrix)");
//    //std::cout << "T1(cabs)\n" << T1iA << std::endl;
//    TArray2 T1ia = _2("<i|T1|a>");
//    //std::cout << "T1(cabs) => i by a block\n" << T1ia << std::endl;
//
//    // recompute E2(CABS) = T1_cabs . H1
//    const double E2_cabs = 2.0 * dot(T1iA("i,A'"), _2("<i|F|A'>"));
//    std::cout << "E2_cabs (recomputed) = " << E2_cabs << std::endl;

#if 0
    // recompute T1_cabs and re-recompute E2_cabs
    {
      typedef TiledArray::Array<T,2> Array;
      Array Fii = _2("<i|F|j>");
      Array FAA = _2("<A'|F|B'>");
      // this computes Z_i^A' = T_i^B' F_B'^A' - F_i^j T_j^A'
      detail::_CABS_singles_h0t1<double> cabs_singles_rhs_eval(FAA, Fii);

      TA::ConjugateGradientSolver<Array, detail::_CABS_singles_h0t1<double> > cg_solver;
      Array FiA = _2("<i|F|A'>");
      Array minus_FiA = -1.0 * FiA("i,A'");
      Array T1_recomp = FiA;

      // make preconditioner: Delta_iA = <i|F|i> - <A'|F|A'>
      typedef detail::diag_precond2<double> pceval_type; //!< evaluator of preconditioner
      typedef TA::Array<T, 2, LazyTensor<T, 2, pceval_type > > TArray2d;
      TArray2d Delta_iA(FiA.get_world(), FiA.trange());
      pceval_type Delta_iA_gen(TA::array_to_eigen(Fii),
                               TA::array_to_eigen(FAA));

      // construct local tiles
      for(auto t=Delta_iA.trange().tiles().begin();
          t!=Delta_iA.trange().tiles().end();
          ++t)
        if (Delta_iA.is_local(*t)) {
          std::array<std::size_t, 2> index; std::copy(t->begin(), t->end(), index.begin());
          madness::Future < typename TArray2d::value_type >
            tile((LazyTensor<T, 2, pceval_type >(&Delta_iA, index, &Delta_iA_gen)
                ));

          // Insert the tile into the array
          Delta_iA.set(*t, tile);
        }
      Array preconditioner = Delta_iA("i,A'");

#if 0
      std::cout << "FiA:\n" << FiA << std::endl;
      std::cout << "Fii:\n" << Fii << std::endl;
      std::cout << "FAA:\n" << FAA << std::endl;
      std::cout << "preconditioner:\n" << preconditioner << std::endl;
#endif

      // solves CABS singles equations T_i^B' F_B'^A' - F_i^j T_j^A' = -F_i^A' using CG
      auto resnorm = cg_solver(cabs_singles_rhs_eval,
                               minus_FiA,
                               T1_recomp,
                               preconditioner,
                               1e-10);
      std::cout << "Converged CG to " << resnorm << std::endl;
      const double E2_cabs = 2.0 * dot(T1_recomp("i,A'"), _2("<i|F|A'>")); // 2 accounts for spin
      std::cout << "E2_cabs (re-recomputed) = " << E2_cabs << std::endl;

      // compute the (2)_S dipole moment
      {
        const double mu_z_e = 2*dot(_2("<i|mu_z|A'>"), T1_recomp("i,A'"))
                            - dot(_2("<i|mu_z|j>"), T1_recomp("i,A'") * T1_recomp("j,A'") )
                            + dot(_2("<A'|mu_z|B'>"), T1_recomp("i,A'") * T1_recomp("i,B'") );
        std::cout << "Mu_z (2)_S = " << 2*mu_z_e << std::endl; // 2 accounts for spin degeneracy
        std::cout << "Mu_z ref = " << (dot(_2("<i|mu_z|j>"),_2("<i|gamma|j>")))
                  << std::endl; // 2 is included in gamma!
      }

      // compute the (2)_S quadrupole moment
      {
        const double q_zz_e = 2*dot(_2("<i|q_zz|A'>"), T1_recomp("i,A'"))
                            - dot(_2("<i|q_zz|j>"), T1_recomp("i,A'") * T1_recomp("j,A'") )
                            + dot(_2("<A'|q_zz|B'>"), T1_recomp("i,A'") * T1_recomp("i,B'") );
        std::cout << "Q_zz (2)_S = " << 2*q_zz_e << std::endl; // 2 accounts for spin degeneracy
        std::cout << "Q_zz ref = " << (dot(_2("<i|q_zz|j>"),_2("<i|gamma|j>")))
                  << std::endl; // 2 is included in gamma!
      }

      // compute orbital rotation multipliers in the (2)_S Lagrangian
      if (1) {
        TArray2 Tia = T1_recomp("j,A'") * _2("<A'|I|a>");
        TArray2 Xia_1 = 2* (_2("<i|F|j>") * Tia("j,a") - T1_recomp("i,B'") * _2("<B'|F|a>"))
                        + 2 * _2("<i|F|C'>") * T1_recomp("j,C'") * Tia("j,a")
                        - 2 * (2 * _4("<j i|g|B' a>") - _4("<j i|g|a B'>") + 2 * _4("<j a|g|B' i>") - _4("<j a|g|i B'>")) * T1_recomp("j,B'")
                        - 2 * (2 * _4("<B' i|g|C' a>") - _4("<B' i|g|a C'>")) * T1_recomp("j,B'") * T1_recomp("j,C'")
                        + 2 * (2 * _4("<j i|g|k a>") - _4("<j i|g|a k>")) * T1_recomp("j,B'") * T1_recomp("k,B'");
        std::cout << "Xia_1, should not be 0:\n" << Xia_1 << std::endl;

        TArray2 Faa = _2("<a|F|b>");
        TArray4 G_ij_ab = _4("<i j|g|a b>");
        TArray4 G_ia_jb = _4("<i a|g|j b>");
        detail::_OrbResponse<double> response_lhs_eval(Faa, Fii, G_ij_ab, G_ia_jb);

        TA::ConjugateGradientSolver<Array, detail::_OrbResponse<double> > cg_solver;
        Array kappa = Xia_1;

        // make preconditioner: Delta_ia = <i|F|i> - <a|F|a>
        typedef detail::diag_precond2<double> pceval_type; //!< evaluator of preconditioner
        typedef TA::Array<T, 2, LazyTensor<T, 2, pceval_type > > TArray2d;
        TArray2d Delta_ia(Xia_1.get_world(), Xia_1.trange());
        pceval_type Delta_ia_gen(TA::array_to_eigen(Fii),
                                 TA::array_to_eigen(Faa));

        // construct local tiles
        for(auto t=Delta_ia.trange().tiles().begin();
            t!=Delta_ia.trange().tiles().end();
            ++t)
          if (Delta_ia.is_local(*t)) {
            std::array<std::size_t, 2> index; std::copy(t->begin(), t->end(), index.begin());
            madness::Future < typename TArray2d::value_type >
              tile((LazyTensor<T, 2, pceval_type >(&Delta_ia, index, &Delta_ia_gen)
                  ));

            // Insert the tile into the array
            Delta_ia.set(*t, tile);
          }
        Array preconditioner = Delta_ia("i,a");

        // solves orbital response using CG
        auto resnorm = cg_solver(response_lhs_eval,
                                 Xia_1,
                                 kappa,
                                 preconditioner,
                                 1e-10);
        std::cout << "Converged CG to " << resnorm << std::endl;

        // verify the solution:
        {
          TArray2 res = Xia_1;
          response_lhs_eval(kappa, res);
          std::cout << "should be zero = " << TA::expressions::norm2(res("i,a") - Xia_1("i,a"));
        }

        const double mu_z_e = dot(kappa("i,a"), _2("<i|mu_z|a>"));
        std::cout << "mu_z_e (orb response) = " << 2*mu_z_e << std::endl;

      }
    }

    // now recompute T1_cabs and E2_cabs using only CABS for perturbation!
    {
      std::cout << "computing E2_cabs due to CABS only as the first-order space\n";
      typedef TiledArray::Array<T,2> Array;
      Array Fii = _2("<i|F|j>");
      Array FAA = _2("<a'|F|b'>");
      // this computes Z_i^a' = T_i^b' F_b'^a' - F_i^j T_j^a'
      detail::_CABS_singles_h0t1<double> cabs_singles_rhs_eval(FAA, Fii);

      TA::ConjugateGradientSolver<Array, detail::_CABS_singles_h0t1<double> > cg_solver;
      Array FiA = _2("<i|F|a'>");
      Array minus_FiA = -1.0 * FiA("i,a'");
      Array T1_recomp = FiA;

      // make preconditioner: Delta_iA = <i|F|i> - <a'|F|a'>
      typedef detail::diag_precond2<double> pceval_type; //!< evaluator of preconditioner
      typedef TA::Array<T, 2, LazyTensor<T, 2, pceval_type > > TArray2d;
      TArray2d Delta_iA(FiA.get_world(), FiA.trange());
      pceval_type Delta_iA_gen(TA::array_to_eigen(Fii),
                               TA::array_to_eigen(FAA));

      // construct local tiles
      for(auto t=Delta_iA.trange().tiles().begin();
          t!=Delta_iA.trange().tiles().end();
          ++t)
        if (Delta_iA.is_local(*t)) {
          std::array<std::size_t, 2> index; std::copy(t->begin(), t->end(), index.begin());
          madness::Future < typename TArray2d::value_type >
            tile((LazyTensor<T, 2, pceval_type >(&Delta_iA, index, &Delta_iA_gen)
                ));

          // Insert the tile into the array
          Delta_iA.set(*t, tile);
        }
      Array preconditioner = Delta_iA("i,a'");

#if 0
      std::cout << "FiA:\n" << FiA << std::endl;
      std::cout << "Fii:\n" << Fii << std::endl;
      std::cout << "FAA:\n" << FAA << std::endl;
      std::cout << "preconditioner:\n" << preconditioner << std::endl;
#endif

      // solves CABS singles equations T_i^B' F_B'^A' - F_i^j T_j^A' = -F_i^A' using CG
      auto resnorm = cg_solver(cabs_singles_rhs_eval,
                               minus_FiA,
                               T1_recomp,
                               preconditioner,
                               1e-10);
      std::cout << "Converged CG to " << resnorm << std::endl;
      const double E2_cabs = 2.0 * dot(T1_recomp("i,a'"), _2("<i|F|a'>")); // 2 accounts for spin
      std::cout << "E2_cabs (re-recomputed) = " << E2_cabs << std::endl;

      // compute the (2)_S dipole moment
      {
        const double mu_z_e = 2*dot(_2("<i|mu_z|a'>"), T1_recomp("i,a'"))
                            - dot(_2("<i|mu_z|j>"), T1_recomp("i,a'") * T1_recomp("j,a'") )
                            + dot(_2("<a'|mu_z|b'>"), T1_recomp("i,a'") * T1_recomp("i,b'") );
        std::cout << "Mu_z (2)_S = " << 2*mu_z_e << std::endl; // 2 accounts for spin degeneracy
      }

      // compute orbital rotation multipliers in the (2)_S Lagrangian
      if (1) {
        // only include the first-order terms in the Lagrangian derivative
        TArray2 Xia_1 = - 2 * T1_recomp("i,b'") * _2("<b'|F|a>")
                      - 2 * (2 * _4("<j i|g|b' a>") - _4("<j i|g|a b'>") + 2 * _4("<j a|g|b' i>") - _4("<j a|g|i b'>")) * T1_recomp("j,b'")
                      - 2 * (2 * _4("<b' i|g|c' a>") - _4("<b' i|g|a c'>")) * T1_recomp("j,b'") * T1_recomp("j,c'")
                      + 2 * (2 * _4("<j i|g|k a>") - _4("<j i|g|a k>")) * T1_recomp("j,b'") * T1_recomp("k,b'");
        std::cout << "Xia_1, should not be 0:\n" << Xia_1 << std::endl;

        TArray2 Faa = _2("<a|F|b>");
        TArray4 G_ij_ab = _4("<i j|g|a b>");
        TArray4 G_ia_jb = _4("<i a|g|j b>");
        detail::_OrbResponse<double> response_lhs_eval(Faa, Fii, G_ij_ab, G_ia_jb);

        TA::ConjugateGradientSolver<Array, detail::_OrbResponse<double> > cg_solver;
        Array kappa = Xia_1;

        // make preconditioner: Delta_ia = <i|F|i> - <a|F|a>
        typedef detail::diag_precond2<double> pceval_type; //!< evaluator of preconditioner
        typedef TA::Array<T, 2, LazyTensor<T, 2, pceval_type > > TArray2d;
        TArray2d Delta_ia(Xia_1.get_world(), Xia_1.trange());
        pceval_type Delta_ia_gen(TA::array_to_eigen(Fii),
                                 TA::array_to_eigen(Faa));

        // construct local tiles
        for(auto t=Delta_ia.trange().tiles().begin();
            t!=Delta_ia.trange().tiles().end();
            ++t)
          if (Delta_ia.is_local(*t)) {
            std::array<std::size_t, 2> index; std::copy(t->begin(), t->end(), index.begin());
            madness::Future < typename TArray2d::value_type >
              tile((LazyTensor<T, 2, pceval_type >(&Delta_ia, index, &Delta_ia_gen)
                  ));

            // Insert the tile into the array
            Delta_ia.set(*t, tile);
          }
        Array preconditioner = Delta_ia("i,a");

        // solves orbital response using CG
        auto resnorm = cg_solver(response_lhs_eval,
                                 Xia_1,
                                 kappa,
                                 preconditioner,
                                 1e-10);
        std::cout << "Converged CG to " << resnorm << std::endl;
        const double mu_z_e = dot(kappa("i,a"), _2("<i|mu_z|a>"));
        std::cout << "mu_z_e (orb reponse) = " << 2*mu_z_e << std::endl;

      }
    }

    // compute T1_obs and E2_obs
    {
      std::cout << "computing E2_obs (virtuals is the first-order space)\n";
      typedef TiledArray::Array<T,2> Array;
      Array Fii = _2("<i|F|j>");
      Array FAA = _2("<a|F|b>");
      // this computes Z_i^a = T_i^b F_b^a - F_i^j T_j^a
      detail::_CABS_singles_h0t1<double> cabs_singles_rhs_eval(FAA, Fii);

      TA::ConjugateGradientSolver<Array, detail::_CABS_singles_h0t1<double> > cg_solver;
      Array FiA = _2("<i|F|a>");
      Array minus_FiA = -1.0 * FiA("i,a");
      Array T1_obs = FiA;

      // make preconditioner: Delta_iA = <i|F|i> - <a|F|a>
      typedef detail::diag_precond2<double> pceval_type; //!< evaluator of preconditioner
      typedef TA::Array<T, 2, LazyTensor<T, 2, pceval_type > > TArray2d;
      TArray2d Delta_iA(FiA.get_world(), FiA.trange());
      pceval_type Delta_iA_gen(TA::array_to_eigen(Fii),
                               TA::array_to_eigen(FAA));

      // construct local tiles
      for(auto t=Delta_iA.trange().tiles().begin();
          t!=Delta_iA.trange().tiles().end();
          ++t)
        if (Delta_iA.is_local(*t)) {
          std::array<std::size_t, 2> index; std::copy(t->begin(), t->end(), index.begin());
          madness::Future < typename TArray2d::value_type >
            tile((LazyTensor<T, 2, pceval_type >(&Delta_iA, index, &Delta_iA_gen)
                ));

          // Insert the tile into the array
          Delta_iA.set(*t, tile);
        }
      Array preconditioner = Delta_iA("i,a");

#if 0
      std::cout << "FiA:\n" << FiA << std::endl;
      std::cout << "Fii:\n" << Fii << std::endl;
      std::cout << "FAA:\n" << FAA << std::endl;
      std::cout << "preconditioner:\n" << preconditioner << std::endl;
#endif

      // solves CABS singles equations T_i^B' F_B'^A' - F_i^j T_j^A' = -F_i^A' using CG
      auto resnorm = cg_solver(cabs_singles_rhs_eval,
                               minus_FiA,
                               T1_obs,
                               preconditioner,
                               1e-10);
      std::cout << "Converged CG to " << resnorm << std::endl;
      // same as the lagrangian since solved for T1 exactly
      const double E2_obs = 2.0 * dot(T1_obs("i,a"), _2("<i|F|a>")); // 2 accounts for spin
      std::cout << "E2_obs = " << E2_obs << std::endl;
      const double E0_obs = dot(_2("<i|F|j>"), _2("<i|gamma|j>")); // 2 included in gamma
      std::cout << "E0_obs = " << scprintf("%15.10lf",E0_obs) << std::endl;

      // compute the (2)_S dipole moment
      {
        const double mu_z_e = 2*dot(_2("<i|mu_z|a>"), T1_obs("i,a"))
                            - dot(_2("<i|mu_z|j>"), T1_obs("i,a") * T1_obs("j,a") )
                            + dot(_2("<a|mu_z|b>"), T1_obs("i,a") * T1_obs("i,b") );
        std::cout << scprintf("Mu_z (2)_S OBS = %25.15lf", 2*mu_z_e) << std::endl; // 2 accounts for spin degeneracy
        std::cout << scprintf("Mu_z ref = %25.15lf", (dot(_2("<i|mu_z|j>"),_2("<i|gamma|j>")))) << std::endl;
      }

      // compute orbital rotation multipliers in the (2)_S OBS Lagrangian
      if (1) {
        // only include the first-order terms in the Lagrangian derivative
        TArray2 Xia_1 = _2("<i|F|j>") * T1_obs("j,a") - T1_obs("i,b") * _2("<b|F|a>")
                        + 2 * _2("<i|F|c>") * T1_obs("j,c") * T1_obs("j,a")
                        + 2 * T1_obs("i,c") * T1_obs("j,c") * _2("<j|F|a>")
                        ;
        std::cout << "Xia_1, should not be 0:\n" << Xia_1 << std::endl;
        TArray2 Fia = _2("<i|F|a>");
        std::cout << "Fia, should be close to Xia_1:\n" << Fia << std::endl;

        TArray2 Faa = _2("<a|F|b>");
        TArray4 G_ij_ab = _4("<i j|g|a b>");
        TArray4 G_ia_jb = _4("<i a|g|j b>");
        detail::_OrbResponse<double> response_lhs_eval(Faa, Fii, G_ij_ab, G_ia_jb);

        TA::ConjugateGradientSolver<Array, detail::_OrbResponse<double> > cg_solver;
        Array kappa = Xia_1;

        // make preconditioner: Delta_ia = <i|F|i> - <a|F|a>
        typedef detail::diag_precond2<double> pceval_type; //!< evaluator of preconditioner
        typedef TA::Array<T, 2, LazyTensor<T, 2, pceval_type > > TArray2d;
        TArray2d Delta_ia(Xia_1.get_world(), Xia_1.trange());
        pceval_type Delta_ia_gen(TA::array_to_eigen(Fii),
                                 TA::array_to_eigen(Faa));

        // construct local tiles
        for(auto t=Delta_ia.trange().tiles().begin();
            t!=Delta_ia.trange().tiles().end();
            ++t)
          if (Delta_ia.is_local(*t)) {
            std::array<std::size_t, 2> index; std::copy(t->begin(), t->end(), index.begin());
            madness::Future < typename TArray2d::value_type >
              tile((LazyTensor<T, 2, pceval_type >(&Delta_ia, index, &Delta_ia_gen)
                  ));

            // Insert the tile into the array
            Delta_ia.set(*t, tile);
          }
        Array preconditioner = Delta_ia("i,a");

        // solves orbital response using CG
        auto resnorm = cg_solver(response_lhs_eval,
                                 Xia_1,
                                 kappa,
                                 preconditioner,
                                 1e-10);
        std::cout << "Converged CG to " << resnorm << std::endl;
        std::cout << scprintf("E_2 (orb reponse) = %25.15lf", 2*dot(kappa("i,a"), _2("<i|F|a>"))) << std::endl;
        const double mu_z_e = dot(kappa("i,a"), _2("<i|mu_z|a>"));
        std::cout << scprintf("mu_z_e (orb reponse) = %25.15lf", 2*mu_z_e) << std::endl;

      }
    }

    if (1) {
      for(int i=0; i<10; ++i) {
        TArray4 should_be_zero = _4("<b j|g|a i>") - _4("<b i|g|a j>");
        std::cout << "should be 0: " << TA::expressions::norminf(should_be_zero("b,j,a,i")) << std::endl;
      }
    }
#endif

    /// this is just an example of how to compute the density
    TArray2 r2_i_j;
    //r2_i_j("i,j") = _4("<i j|r|p q>") * _4("<k_F(p) j|r|p q>");

    return r2_i_j;
#else // ENABLE_SRR12_RDM1
    MPQC_ASSERT(false); // not converted yet to the TiledArray expressions branch
    return TArray2();
#endif
  }

  template <typename T>
  typename SingleReference_R12Intermediates<T>::TArray4
  SingleReference_R12Intermediates<T>::V_spinfree(bool symmetrize_p1_p2) {

    TArray4 V_ij_mn;
    V_ij_mn("i1,i2,m1,m2") = _4("<i1 i2|gr|m1 m2>") - _4("<i1 i2|r|p1 p2>") * _4("<m1 m2|g|p1 p2>")
                    - _4("<i1 i2|r|m3_gamma(m) a'>") * _4("<m1 m2|g|m3 a'>");

    if (symmetrize_p1_p2)
      V_ij_mn("i1,i2,m1,m2") = 0.5 * (V_ij_mn("i1,i2,m1,m2") + V_ij_mn("i2,i1,m2,m1"));

    return V_ij_mn;
  }

  template <typename T>
  typename SingleReference_R12Intermediates<T>::TArray4
  SingleReference_R12Intermediates<T>::X_spinfree(bool symmetrize_p1_p2) {

    TArray4 X_ij_kl;
    X_ij_kl("i1,i2,j1,j2") = _4("<i1 i2|r2|j1 j2>") - _4("<i1 i2|r|p1 p2>") * _4("<j1 j2|r|p1 p2>")
                    - _4("<i1 i2|r|m3_gamma(m) a'>") * _4("<j1 j2|r|m3 a'>");

    if (symmetrize_p1_p2)
      X_ij_kl("i1,i2,j1,j2") = 0.5 * (X_ij_kl("i1,i2,j1,j2") + X_ij_kl("i2,i1,j2,j1"));

    return X_ij_kl;
  }

  template <typename T>
  typename SingleReference_R12Intermediates<T>::TArray4
  SingleReference_R12Intermediates<T>::B_spinfree(bool symmetrize_p1_p2) {

    TArray4 B_ij_kl;
    B_ij_kl("i1,i2,j1,j2") =

    // everything seems scaled up by factor of 2 relative to Eq.(12) in J. Chem. Phys. 135, 214105 (2011),
    // due to including particle 1 and particle 2 contributions?

    // diag                      Q
        _4("<i1 i2|rTr|j1 j2>") + 2.0 * _4("<i1 i2|r2|j1 j2_hJ(p')>")

        //           rKr
            - 2.0 * _4("<i1 i2|r|r' s'>") * _4("<j1 j2|r|r' s'_K(p')>")

            //           rFr
            - 2.0 * _4("<i1 i2|r|r s>") * _4("<j1 j2|r|r s_F(p)>")

            //           rFr_2, extra 2 due to bra-ket symmetrization
            - 4.0 * _4("<i1 i2|r|r s>") * _4("<j1 j2|r|r s_F(a')>")

            //           rFGr
            - _4("<i1 i2|r|n_gamma(m) b'>") * _4("<j1 j2|r|n b'_F(a')>")

            //           rFGr_2
            - _4("<i1 i2|r|n_gamma(m) a'>") * _4("<j1 j2|r|n_F(p') a'>");

    B_ij_kl("i1,i2,j1,j2") = 0.5
        * (B_ij_kl("i1,i2,j1,j2") + B_ij_kl("i2,i1,j2,j1"));
    B_ij_kl("i1,i2,j1,j2") = 0.5
        * (B_ij_kl("i1,i2,j1,j2") + B_ij_kl("j1,j2,i1,i2"));

    if (symmetrize_p1_p2)
      B_ij_kl("i1,i2,j1,j2") = 0.5
          * (B_ij_kl("i1,i2,j1,j2") + B_ij_kl("i2,i1,j2,j1"));

    return B_ij_kl;
  }


}; // end of namespace sc

#endif // end of header guard


// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
