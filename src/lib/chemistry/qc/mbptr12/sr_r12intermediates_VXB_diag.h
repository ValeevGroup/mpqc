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

#include <TiledArray/eigen.h>
#include <TiledArray/algebra/conjgrad.h>
#include <TiledArray/expressions.h>

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

  inline TA::Array<double, 2 > XaiAddToXam(const TA::Array<double, 2 >& Xam,
                                           const TA::Array<double, 2 >& Xai) {
    MPQC_ASSERT(Xam.size() == 1);

    typedef TA::Array<double, 2> TArray2;
    typename TArray2::value_type tile_X =
    static_cast<const typename TArray2::value_type&>(Xam.find(0)).clone();
    typename TArray2::value_type tile_Xai = Xai.find(0);
    const std::size_t zero_cols = tile_X.range().size()[1] - tile_Xai.range().size()[1];

    std::array<std::size_t, 2> i = {{0, 0}};
    std::size_t ix = 0;
    for(i[0] = tile_X.range().start()[0];
      i[0] < tile_X.range().finish()[0]; ++i[0]) {
      for(i[1] = zero_cols + tile_X.range().start()[1];
        i[1] < tile_X.range().finish()[1]; ++i[1], ++ix) {
        tile_X[i] +=  tile_Xai[ix];
      }
    }

    TArray2 X(Xam.get_world(), Xam.trange());
    X.set(0, tile_X);
    return X;
  }

  template <typename T>
  typename SingleReference_R12Intermediates<T>::TArray2
  SingleReference_R12Intermediates<T>::Bpk_qk(const char* p, const char* q) {

    const double C_0 = 1.0 / 2.0;
    const double C_1 = 1.0 / 4.0;
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

    TArray2 B_pq =
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
         -  ( B_C1 * (_4(pk_Ab.c_str()) * _4(qk_AFb.c_str())
                    + _4(kp_Ab.c_str()) * _4(kq_AFb.c_str()))
            + B_C2 * (_4(kp_Ab.c_str()) * _4(qk_AFb.c_str())
                    + _4(pk_Ab.c_str()) * _4(kq_AFb.c_str()))

            + B_C1 * (_4(qk_Ab.c_str()) * _4(pk_AFb.c_str())
                    + _4(kq_Ab.c_str()) * _4(kp_AFb.c_str()))
            + B_C2 * (_4(kq_Ab.c_str()) * _4(pk_AFb.c_str())
                    + _4(qk_Ab.c_str()) * _4(kp_AFb.c_str()))
            )
          ;
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

    TArray4 Bpr_qs =
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
            )
          ;
    return Bpr_qs;
  }

  template <typename T>
  typename SingleReference_R12Intermediates<T>::TArray4
  SingleReference_R12Intermediates<T>::Vpq_rs(const char* p, const char* q,
                                              const char* r, const char* s)
  {

    const double C_0 = 0.5;
    const double C_1 = 0.25;

    const double V_C1 = (0.5 * C_0 + 1.5 * C_1);
    const double V_C2 = (0.5 * C_0 - 1.5 * C_1);

    std::string gr_pqrs = std::string("<") + p + " " + q + "|gr|" + r + " " + s + ">";
    std::string gr_qprs = std::string("<") + q + " " + p + "|gr|" + r + " " + s + ">";

    std::string rpq_pq = std::string("<") + p + " " + q + "|r|p q>";
    std::string rqp_pq = std::string("<") + q + " " + p + "|r|p q>";
    std::string grs_pq = std::string("<p q|g|") + r + " " + s + ">";

    std::string rpq_apn = std::string("<") + p + " " + q + "|r|a' n>";
    std::string rqp_apn = std::string("<") + q + " " + p + "|r|a' n>";
    std::string grs_apn = std::string("<a' n|g|") + r + " " + s + ">";

    std::string rpq_nap = std::string("<") + p + " " + q + "|r|n a'>";
    std::string rqp_nap = std::string("<") + q + " " + p + "|r|n a'>";
    std::string grs_nap = std::string("<n a'|g|") + r + " " + s + ">";

    TArray4 Vpq_rs =  V_C1 * _4(gr_pqrs.c_str()) + V_C2 * _4(gr_qprs.c_str())
                   - (V_C1 * _4(rpq_pq.c_str()) + V_C2 * _4(rqp_pq.c_str()))
                     * _4(grs_pq.c_str())
                   - (V_C1 * _4(rpq_apn.c_str()) + V_C2 * _4(rqp_apn.c_str()))
                     * _4(grs_apn.c_str())
                   - (V_C1 * _4(rpq_nap.c_str()) + V_C2 * _4(rqp_nap.c_str()))
                     * _4(grs_nap.c_str())
                   ;
    return Vpq_rs;
  }

  template <typename T>
  typename SingleReference_R12Intermediates<T>::TArray2
  SingleReference_R12Intermediates<T>::Vrk_sk(const char* r, const char* s)
  {

    const double C_0 = 0.5;
    const double C_1 = 0.25;

    const double V_C1 = (0.5 * C_0 + 1.5 * C_1);
    const double V_C2 = (0.5 * C_0 - 1.5 * C_1);

    std::string gr_rksk = std::string("<") + r + " k|gr|" + s + " l>";
    std::string gr_krsk = std::string("<k ") + r + "|gr|" + s + " l>";

    std::string r_rkpq = std::string("<") + r + " k|r|p q>";
    std::string r_krpq = std::string("<k ") + r + "|r|p q>";
    std::string g_skpq = std::string("<p q|g|") + s + " k>";

    std::string r_rkapn = std::string("<") + r + " k|r|a' n>";
    std::string r_krapn = std::string("<k ") + r + "|r|a' n>";
    std::string g_skapn = std::string("<a' n|g|") + s + " k>";

    std::string r_rknap = std::string("<") + r + " k|r|n a'>";
    std::string r_krnap = std::string("<k ") + r + "|r|n a'>";
    std::string g_sknap = std::string("<n a'|g|") + s + " k>";

    TArray2 Vrk_sk = (V_C1 * _4(gr_rksk.c_str()) + V_C2 * _4(gr_krsk.c_str()))
                     * _2("<l|I|k")
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

    TArray2 V_ij_ij_cabs = dotket(ij_xy("<i j|g|m a'>"), ij_xy("<i j|r|m a'>"));
    // don't need this in closed shell!
    //TArray2 V_ij_ij_cabs1 = dotket(ij_xy("<i j|g|a' m>"), ij_xy("<i j|r|a' m>"));

    TArray2 V_ij_ij = take(ij_xy("<i j|gr|p q>"), ij) - dotket(ij_xy("<i j|g|p q>"), ij_xy("<i j|r|p q>"))
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
          R1 = K_bj("b,j") * A_bjai("b,j,a,i");
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
    // D^m_n =  t^m_A' * t^A'_n
    TArray2 D_e2_mn = TmA("m,A'") * TmA("n,A'");
    // D^A'_B' = t^A'_m * t^m_B'
    TArray2 D_e2_AB = TmA("m,A'") * TmA("m,B'");

    TArray4 g_aAmn = _4("<a A'|g|m n>");
    TArray4 g_ammn = _4("<a m|g|m1 n>");
    TArray2 Xam_E2 = 2.0 * (
                     - _2("<A'|F|a>") * TmA("m,A'")
                     + Tma("n,a") * _2("<m|F|n>")
                     + (TmA("n,A'") * Tma("n,a")) * _2("<m|F|A'>")
                     //
                     - ( 2.0 * _4("<a n|g|m A'>") - _4("<a n|g|A' m>")
                       + 2.0 * g_aAmn("a,A',m,n") - g_aAmn("a,A',n,m")
                       ) * TmA("n,A'")
                     //
                     - (2.0 * _4("<a A'|g|m B'>") - _4("<a A'|g|B' m>"))
                       * D_e2_AB("A',B'")
                     //
                     + (2.0 * g_ammn("a,n1,m,n2") - g_ammn("a,n1,n2,m"))
                       * D_e2_mn("n1,n2")
                     );
    return Xam_E2;
  }

  // Xam contribution from MP2
  template <typename T>
  typename SingleReference_R12Intermediates<T>::TArray2
  SingleReference_R12Intermediates<T>::Xam_mp2(const TArray4& T2_ijab,
                                               const TArray2& Dij,
                                               const TArray2& Dab) {

    // Xai contribution from MP2 part
    TArray2 Xai_mp2 = 2.0 * (
                      // derivatives of 1/4 g^kl_cd T^cd_kl
                      // g^al_cd T^cd_il
                      - _4("<a l|g|c d>")
                        * (2.0 * T2_ijab("i,l,c,d") - T2_ijab("l,i,c,d"))
                      );
    TArray2 Xam_mp2 = 2.0 * (
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

//    const std::size_t zero_cols = Xam_mp2.range().size()[1] - Xai_mp2.range().size()[1];
//    for(auto t = X_mp2.begin(); t != X_mp2.end(); ++t) {
//      typename TArray2::range_type::index i = t.index();
//      typename TArray2::value_type tile = Xam_mp2.find(i);
//
//      if(i[1] >= zero_cols) {
//        std::array<std::size_t, 2> source_index = {{i[0], i[1] - zero_cols}};
//        const typename TArray2::value_type source_t = Xai_mp2.find(source_index);
//        *t = typename TArray2::value_type(tile.range(), tile.begin(), source_t.begin(), std::plus<double>());
//      } else {
//        *t = tile.clone();
//      }
//    }
    TArray2 X_mp2 = XaiAddToXam(Xam_mp2, Xai_mp2);

    return X_mp2;
  }

  // Xai contribution from MP2 F12 coulping part
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

    // Xai contribution from MP2 part & MP2 F12 coulping part
    TArray2 Xai_mp2f12 =
                        // derivatives of 1/4 g^kl_cd \tilde{T}^cd_kl
                        // g^al_cd \tilde{T}^cd_il
                        -  _4("<a l|g|c d>")
                           * (R_C1 * A_ijab("i,l,c,d") + R_C2 * A_ijab("l,i,c,d"))

                        // F^a'_b R^al_a'c \tilde{T}^bc_il
                      - _4("<a l|r|b_F(a') c>")
                        * (  R_C1 * T2_ijab("i,l,b,c") + R_C2 * T2_ijab("l,i,b,c")
                           + RR_C1 * A_ijab("i,l,b,c") + RR_C2 * A_ijab("l,i,b,c"))
                      - _4("<a l|r|c b_F(a')>")
                        * (  R_C1 * T2_ijab("i,l,c,b") + R_C2 * T2_ijab("l,i,c,b")
                           + RR_C1 * A_ijab("i,l,c,b") + RR_C2 * A_ijab("l,i,c,b"));
    TArray2 Xam_mp2f12 =
                        // derivatives of 1/4 g^kl_cd \tilde{T}^cd_kl
                        // g_ic^kl \tilde{T}^Ac_kl
                          (R_C1 * A_ijab("k,l,a,c") + R_C2 * A_ijab("k,l,c,a"))
                          * _4("<m c|g|k l>")

                        // derivatives of 1/2 F^c_d D^d_c - 1/2 F^l_k D^k_l
                        // g^ac_md D^d_c
                        - (2.0 * _4("<a c|g|m d>") - _4("<a c|g|d m>"))
                          * Dab("d,c")
                        // g^al_mk D^k_l
                        + (2.0 * _4("<a l|g|m k>") - _4("<a l|g|k m>"))
                          * Dij("k,l")

                        // derivatives of 1/2 F^a'_b R^kl_a'c \tilde{T}^bc_kl
                        //   1/2 F^a'_m R^kl_a'b \tilde{T}^ab_kl
                        // & 1/2 F^a'_b R^kl_ma' \tilde{T}^ab_kl
                        + (  R_C1 * T2_ijab("k,l,a,b") + R_C2 * T2_ijab("k,l,b,a")
                           + RR_C1 * A_ijab("k,l,a,b") + RR_C2 * A_ijab("k,l,b,a"))
                          * (_4("<k l|r|m_F(a') b>") + _4("<k l|r|m b_F(a')>"))

                       //   g^aa'_mb 1/2 R^kl_a'c \tilde{T}^bc_kl
                       // & g^ab_ma' 1/2 R^kl_a'c \tilde{T}^bc_kl
                        - (  2.0 * _4("<a a'|g|m b>") - _4("<a a'|g|b m>")
                           + 2.0 * _4("<a b|g|m a'>") - _4("<a b|g|a' m>"))
                          * RT_apb("a',b");

    //std::cout << std::endl << "Xam_mp2f12:\n" << Xam_mp2f12 << std::endl;
    TArray2 X_mp2f12 = XaiAddToXam(Xam_mp2f12, Xai_mp2f12);

    return X_mp2f12;
  }

  // Xam contribution from F12 V part
  template <typename T>
  typename SingleReference_R12Intermediates<T>::TArray2
  SingleReference_R12Intermediates<T>::Xam_V(const double C_0, const double C_1) {

    const double R_C1 = (0.5 * C_0 + 1.5 * C_1);
    const double R_C2 = (0.5 * C_0 - 1.5 * C_1);

    TArray4 gr_ak_ij = _4("<a k|gr|i j>");
    TArray4 rpq_kl = _4("<p q|r|k l>");
    TArray4 rapn_kl = _4("<a' n|r|k l>");
    TArray4 rnap_kl = _4("<n a'|r|k l>");

    TArray4 raap_kl = _4("<a a'|r|k l>");
    TArray4 rmap_kl = _4("<m a'|r|k l>");
    TArray2 Ikl = _2("<k|I|l>");

    TArray4 rpq_ak = _4("<p q|r|a k>");
    TArray2 Xai_V = // 1/2 R^ik_AC g^AC_ak
                  - ( (R_C1 * gr_ak_ij("a,k,i,l") + R_C2 * gr_ak_ij("a,k,l,i"))
                      * Ikl("k,l")
                    - _4("<a k|g|p q>")
                       * ( R_C1 * rpq_kl("p,q,i,k") + R_C2 * rpq_kl("p,q,k,i") )
                    - _4("<a k|g|a' n>")
                       * ( R_C1 * rapn_kl("a',n,i,k") + R_C2 * rapn_kl("a',n,k,i"))
                    - _4("<a k|g|n a'>")
                       * ( R_C1 * rnap_kl("n,a',i,k") + R_C2 * rnap_kl("n,a',k,i"))
                    )
                    // 1/2 R^ak_AC g^AC_ik
                  - ( (R_C1 * gr_ak_ij("a,l,i,k") + R_C2 * gr_ak_ij("a,l,k,i"))
                      * Ikl("k,l")
                  - _4("<i k|g|p q>")
                    * ( R_C1 * rpq_ak("p,q,a,k") + R_C2 * rpq_ak("q,p,a,k"))
                  - _4("<i k|g|a' n>")
                     * ( R_C1 * _4("<a' n|r|a k>") + R_C2 * _4("<a' n|r|k a>"))
                  - _4("<i k|g|n a'>")
                     * ( R_C1 * _4("<n a'|r|a k>") + R_C2 * _4("<n a'|r|k a>"))
                    );
    TArray2 Xam_V =  // 1/2 R^kl_ac' g^mc'_kl
                    ( R_C1 * raap_kl("a,a',k,l") + R_C2 * raap_kl("a,a',l,k"))
                    * _4("<k l|g|m a'>")
                    // 1/2 R^kl_mc' g^ac'_kl
                  + ( R_C1 * rmap_kl("m,a',k,l") + R_C2 * rmap_kl("m,a',l,k"))
                    * _4("<k l|g|a a'>");
    TArray2 X_V = XaiAddToXam(Xam_V, Xai_V);

    return X_V;
  }

  // Xam contribution from F12 X part
  template <typename T>
  typename SingleReference_R12Intermediates<T>::TArray2
  SingleReference_R12Intermediates<T>::Xam_X(const double C_0, const double C_1) {

    const double RR_C1 = 0.5 * C_0 * C_0 + 1.5 * C_1 * C_1;
    const double RR_C2 = 0.5 * C_0 * C_0 - 1.5 * C_1 * C_1;

    TArray4 r2_akjl = _4("<a k|r2|j l>");
    TArray4 r2_kajl = _4("<k a|r2|j l>");
    TArray2 F_ij = _2("<i|F|j>");

    TArray4 r_pqjk = _4("<p q|r|j k>");
    TArray4 r_apnjk = _4("<a' n|r|j k>");
    TArray4 r_napjk = _4("<n a'|r|j k>");

    // - 1/2 (F^i_i + F^k_k) R^ak_alpha beta R_ik^alpha beta
    TArray2 X_ai =  - (RR_C1 * r2_akjl("a,k,j,l") + RR_C2 * r2_kajl("k,a,j,l"))
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
    TArray4 r_aapkl = _4("<a a'|r|k l>");
    TArray4 r_klmap = _4("<k l|r|m a'>");
    TArray2 X_am = (RR_C1 * r_aapkl("a,b',k,l") + RR_C2 * r_aapkl("a,b',l,k"))
                    * ( r_klmap("j,l,m,b'") * F_ij("j,k")
                      + r_klmap("k,j,m,b'") * F_ij("j,l"));
    TArray2 X = XaiAddToXam(X_am, X_ai);

    return X;
  }

  // Xam contribution from F12 B part
  template <typename T>
  typename SingleReference_R12Intermediates<T>::TArray2
  SingleReference_R12Intermediates<T>::Xam_B(const double C_0, const double C_1) {

    const double B_C1 = 0.5 * C_0 * C_0 + 1.5 * C_1 * C_1;
    const double B_C2 = 0.5 * C_0 * C_0 - 1.5 * C_1 * C_1;

    // 1st parts of B_ai:
    // - R^ak_A'B' F^A'_C' R^ik_C'B' (A': all virtual orbital index)

    TArray4 r_akPQ = _4("<a k|r|p' q'>");
    TArray4 r_kaPQ = _4("<k a|r|p' q'>");
    TArray4 rik_PKQ = _4("<i k|r|p'_K(r') q'>");

    TArray4 r_akPn = _4("<a k|r|p' n>");
    TArray4 r_kaPn = _4("<k a|r|p' n>");
    TArray4 rik_PFn = _4("<i k|r|p'_F(r') n>");

    TArray4 r_akmA = _4("<a k|r|m a'>");
    TArray4 r_kamA = _4("<k a|r|m a'>");
    TArray4 rik_mFA = _4("<i k|r|m_F(n) a'>");

    TArray4 r_akpq = _4("<a k|r|p b>");
    TArray4 r_kapq = _4("<k a|r|p b>");
    TArray4 rik_pFb = _4("<i k|r|p_F(r) b>");

    TArray4 rik_nFA = _4("<i k|r|n_F(p') a'>");
    //
    TArray4 r_ikmA = _4("<i k|r|m a'>");
    TArray4 rak_nFA = _4("<a k|r|n_F(p') a'>");
    TArray4 rka_nFA = _4("<k a|r|n_F(p') a'>");

    TArray4 r_akAb  = _4("<a k|r|a' b>");
    TArray4 r_kaAb  = _4("<k a|r|a' b>");
    TArray4 rik_AFb = _4("<i k|r|a'_F(q) b>");
    //
    TArray4 r_ikAb  = _4("<i k|r|a' b>");
    TArray4 rak_AFb = _4("<a k|r|a'_F(q) b>");
    TArray4 rka_AFb = _4("<k a|r|a'_F(q) b>");

    TArray2 B_ai = //          diag
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

    // 2nd parts of B_ai:
    //    1/2 R^ab'_kl F^c'_m R^c'b'_kl
    // + (1/2 R^mb'_kl F^C'_a R^C'b'_kl + 1/2 R^mb'_kl F^c'_b' R^ab'_kl)
    TArray4 r_abpkl = _4("<a b'|r|k l>");
    TArray4 raFpbp_kl = _4("<a_F(c') b'|r|k l>");
    TArray4 raFbp_kl = _4("<a_F(c) b'|r|k l>");
    TArray4 rabpFp_kl = _4("<a b'_F(c')|r|k l>") ;

    TArray2 B_am = (B_C1 * r_abpkl("a,b',k,l") + B_C2 * r_abpkl("a,b',l,k"))
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

    TArray2 X_B = XaiAddToXam(B_am, B_ai);

    return X_B;
  }

  template <typename T>
  typename SingleReference_R12Intermediates<T>::TArray2
  SingleReference_R12Intermediates<T>::rdm1() {

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

    bool compute_dipole = false;
    bool compute_quadrupole = false;
    bool compute_EFG = true;

    // Obtain property integrals which is needed in multiple places
    // dipole integrals
    TArray2 mu_z_mn = (compute_dipole ? _2("<m|mu_z|n>"): TArray2());
    TArray2 mu_z_am = (compute_dipole ? _2("<a|mu_z|m>") : TArray2());

    // quadrupole integrals
    TArray2 q_xx_mn = (compute_quadrupole ? _2("<m|q_xx|n>") : TArray2());
    TArray2 q_yy_mn = (compute_quadrupole ? _2("<m|q_yy|n>") : TArray2());
    TArray2 q_zz_mn = (compute_quadrupole ? _2("<m|q_zz|n>") : TArray2());
    TArray2 Qxx_mn = (compute_quadrupole ?
                      q_xx_mn("m,n") - (q_zz_mn("m,n") + q_yy_mn("m,n")) * 0.5 : TArray2());
    TArray2 Qyy_mn = (compute_quadrupole ?
                      q_yy_mn("m,n") - (q_zz_mn("m,n") + q_xx_mn("m,n")) * 0.5 : TArray2());
    TArray2 Qzz_mn = (compute_quadrupole ?
                      q_zz_mn("m,n") - (q_xx_mn("m,n") + q_yy_mn("m,n")) * 0.5 : TArray2());
    TArray2 Qxz_mn = (compute_quadrupole ? _2("<m|q_xz|n>") * 1.5 : TArray2());
    TArray2 Qxy_mn = (compute_quadrupole ? _2("<m|q_xy|n>") * 1.5 : TArray2());
    TArray2 Qyz_mn = (compute_quadrupole ? _2("<m|q_yz|n>") * 1.5 : TArray2());

    TArray2 q_xx_am = (compute_quadrupole ? _2("<a|q_xx|m>") : TArray2());
    TArray2 q_yy_am = (compute_quadrupole ? _2("<a|q_yy|m>") : TArray2());
    TArray2 q_zz_am = (compute_quadrupole ? _2("<a|q_zz|m>") : TArray2());
    TArray2 Qxx_am = (compute_quadrupole ?
                      q_xx_am("a,m") - (q_zz_am("a,m") + q_yy_am("a,m")) * 0.5 : TArray2());
    TArray2 Qyy_am = (compute_quadrupole ?
                      q_yy_am("a,m") - (q_zz_am("a,m") + q_xx_am("a,m")) * 0.5 : TArray2());
    TArray2 Qzz_am = (compute_quadrupole ?
                      q_zz_am("a,m") - (q_xx_am("a,m") + q_yy_am("a,m")) * 0.5 : TArray2());
    TArray2 Qxy_am = (compute_quadrupole ? _2("<a|q_xy|m>") * 1.5 : TArray2());
    TArray2 Qxz_am = (compute_quadrupole ? _2("<a|q_xz|m>") * 1.5 : TArray2());
    TArray2 Qyz_am = (compute_quadrupole ? _2("<a|q_yz|m>") * 1.5 : TArray2());

    // electric field gradient integrals
    TArray2 v_xx_mn = (compute_EFG ? _2("<m|ddphi_xx|n>") : TArray2());
    TArray2 v_yy_mn = (compute_EFG ? _2("<m|ddphi_yy|n>") : TArray2());
    TArray2 v_zz_mn = (compute_EFG ? _2("<m|ddphi_zz|n>") : TArray2());
    TArray2 v_xy_mn = (compute_EFG ? _2("<m|ddphi_xy|n>") : TArray2());
    TArray2 v_xz_mn = (compute_EFG ? _2("<m|ddphi_xz|n>") : TArray2());
    TArray2 v_yz_mn = (compute_EFG ? _2("<m|ddphi_yz|n>") : TArray2());

    TArray2 v_xx_am = (compute_EFG ? _2("<a|ddphi_xx|m>") : TArray2());
    TArray2 v_yy_am = (compute_EFG ? _2("<a|ddphi_yy|m>") : TArray2());
    TArray2 v_zz_am = (compute_EFG ? _2("<a|ddphi_zz|m>") : TArray2());
    TArray2 v_xy_am = (compute_EFG ? _2("<a|ddphi_xy|m>") : TArray2());
    TArray2 v_xz_am = (compute_EFG ? _2("<a|ddphi_xz|m>") : TArray2());
    TArray2 v_yz_am = (compute_EFG ? _2("<a|ddphi_yz|m>") : TArray2());

    // compute HF contribution
    Ref<Molecule> mol = r12world_->basis()->molecule();
    const int natom = mol->natom();

    TArray2 Imn = _2("<m|I|n>");

    // HF eletric dipole
    if (compute_dipole) {
      // Nuclear contribution to dipole
      double mu_z_n = 0.0;
      for(int a = 0; a < natom; ++a) {
        const double x = mol->r(a, 0);
        const double y = mol->r(a, 1);
        const double z = mol->r(a, 2);
        const double Z = mol->Z(a);
        mu_z_n += Z * z;
      }

      // SCF contribution to electronic electric dipole
      const double mu_z_scf = dot(mu_z_mn("m,n"), Imn("m,n"));

      std::cout << std::endl
                << "mu_z (HF=SCF+N) = " << scprintf("%12.10f", mu_z_n - mu_z_scf * 2.0) //electron charge = -1, hence the minus
                << std::endl << std::endl;
    }

    // HF eletric quadrupole
    if (compute_quadrupole) {
      // Nuclear electric quadrupole (traceless)
      double q_xx_n = 0.0, q_yy_n = 0.0, q_zz_n = 0.0;
      double q_xy_n = 0.0, q_xz_n = 0.0, q_yz_n = 0.0;
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
      const double q_xx_scf = dot(Qxx_mn("m,n"), Imn("m,n"));
      const double q_yy_scf = dot(Qyy_mn("m,n"), Imn("m,n"));
      const double q_zz_scf = dot(Qzz_mn("m,n"), Imn("m,n"));
      const double q_xy_scf = dot(Qxy_mn("m,n"), Imn("m,n"));
      const double q_xz_scf = dot(Qxz_mn("m,n"), Imn("m,n"));
      const double q_yz_scf = dot(Qyz_mn("m,n"), Imn("m,n"));

      // HF electronic electric quadrupoles
      std::cout << std::endl
                << "traceless quadrupole moment (HF=SCF+N)" << std::endl
                << "q_xx (HF) = " << scprintf("%12.10f", q_xx_n - q_xx_scf * 2.0)
                << "  q_yy (HF) = " << scprintf("%12.10f", q_yy_n - q_yy_scf * 2.0)
                << "  q_zz (HF) = " << scprintf("%12.10f", q_zz_n - q_zz_scf * 2.0)
                << std::endl
                << "q_xy (HF) = " << scprintf("%12.10f", q_xy_n - q_xy_scf * 2.0)
                << "  q_xz (HF) = " << scprintf("%12.10f", q_xz_n - q_xz_scf * 2.0)
                << "  q_yz (HF) = " << scprintf("%12.10f", q_yz_n- q_yz_scf * 2.0)
                << std::endl;
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

      std::cout << std::endl
                << "electric gradient (HF=SCF+N)" << std::endl
                << "v_xx (HF) = " << scprintf("%12.10f", v_xx_n + v_xx_scf * 2.0)
                << "  v_yy (HF) = " << scprintf("%12.10f", v_yy_n + v_yy_scf * 2.0)
                << "  v_zz (HF) = " << scprintf("%12.10f", v_zz_n + v_zz_scf * 2.0)
                << std::endl
                << "v_xy (HF) = " << scprintf("%12.10f", v_xy_n + v_xy_scf * 2.0)
                << "  v_xz (HF) = " << scprintf("%12.10f", v_xz_n + v_xz_scf * 2.0)
                << "  v_yz (HF) = " << scprintf("%12.10f", v_yz_n + v_yz_scf * 2.0)
                << std::endl;
    }

    // compute integrals needed for orbital relaxation
    // i.e. solve Abnam Dbn = Xam (close-shell formula)
    TArray4 g_mnab = _4("<m n|g|a b>");
    TArray4 A_bnam = - _4("<b n|g|a m>") - g_mnab("m,n,b,a") + 4.0 * g_mnab("n,m,b,a")
                     + _2("<b|F|a>") * Imn("m,n") - _2("<a|I|b>") * _2("<m|F|n>");

    // Make preconditioner: Delta_am = 1 / (<a|F|a> - <m|F|m>) for
    // solving k_bn A_bnam = X_am
    TArray2 mFmn = - _2("<m|F|n>");
    TArray2 mFab = - _2("<a|F|b>");
    typedef detail::diag_precond2<double> pceval_type; //!< evaluator of preconditioner
    pceval_type Delta_am_gen(TA::array_to_eigen(mFab), TA::array_to_eigen(mFmn));

    typedef TA::Array<T, 2, LazyTensor<T, 2, pceval_type > > TArray2d;
    TArray2 Iam = _2("<a|I|m>");
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
    TArray2 preconditioner = Delta_am("a,m");

    detail::Orbital_relaxation_Abjai<double> Orbital_relaxation_Abnam(A_bnam);
    TA::ConjugateGradientSolver<TiledArray::Array<T,2>,
                                detail::Orbital_relaxation_Abjai<double> > cg_solver2;

    // compute CABS singles contribution
#if 1
    {
    TArray2 TmA = _2("<m|T1|A'>");
    TArray2 Tma = _2("<m|T1|a>");

    // density from CABS Singles contribution
    // D^m_n =  t^m_A' * t^A'_n
    TArray2 D_e2_mn = TmA("m,A'") * TmA("n,A'");
    // D^A'_B' = t^A'_m * t^m_B'
    TArray2 D_e2_AB = TmA("m,A'") * TmA("m,B'");
    // D^A'_m = t^A'_m

    // CABS singles orbital response
    TArray2 Xam_E2 = Xam_CabsSingles(TmA, Tma);
    TArray2 Dbn_E2(Xam_E2.get_world(), Xam_E2.trange());
    // solve k_bn A_bnam = X_am
    auto resnorm_E2 = cg_solver2(Orbital_relaxation_Abnam,
                                 Xam_E2,
                                 Dbn_E2,
                                 preconditioner,
                                 1e-10);

    // compute CABS singles contribution to dipole
    if (compute_dipole) {
      TArray2 mu_z_AB = _2("<A'|mu_z|B'>");
      TArray2 mu_z_mA = _2("<m|mu_z|A'>");

      const double mu_z_e2 = - dot(mu_z_mn("m,n"), D_e2_mn("m,n"))
                             + dot(mu_z_AB("A',B'"), D_e2_AB("A',B'"))
                             + dot(mu_z_mA("m,A'"), TmA("m,A'")) * 2.0;
      std::cout << std::endl << "mu_z (E2) = "
                << scprintf("%12.10f", - mu_z_e2 * 2.0) << std::endl;

      // CABS singles orbital response contribution to dipole
      const double mu_z_E2 = dot(mu_z_am("a,m"), Dbn_E2("a,m"));
      std::cout << std::endl << std::endl
                << "mu_z (E2 orbital response) = " << scprintf("%12.10f", - mu_z_E2 * 2.0)
                << std::endl << std::endl;
    }
    world_.gop.fence();

    // CABS singles contribution to quadrupoles
    if (compute_quadrupole) {
      TArray2 q_xx_AB = _2("<A'|q_xx|B'>");
      TArray2 q_yy_AB = _2("<A'|q_yy|B'>");
      TArray2 q_zz_AB = _2("<A'|q_zz|B'>");
      TArray2 Qxx_AB = q_xx_AB("A',B'") - (q_zz_AB("A',B'") + q_yy_AB("A',B'")) * 0.5;
      TArray2 Qyy_AB = q_yy_AB("A',B'") - (q_xx_AB("A',B'") + q_zz_AB("A',B'")) * 0.5;
      TArray2 Qzz_AB = q_zz_AB("A',B'") - (q_xx_AB("A',B'") + q_yy_AB("A',B'")) * 0.5;
      TArray2 Qxy_AB = _2("<A'|q_xy|B'>") * 1.5;
      TArray2 Qxz_AB = _2("<A'|q_xz|B'>") * 1.5;
      TArray2 Qyz_AB = _2("<A'|q_yz|B'>") * 1.5;

      TArray2 q_xx_mA = _2("<m|q_xx|A'>");
      TArray2 q_yy_mA = _2("<m|q_yy|A'>");
      TArray2 q_zz_mA = _2("<m|q_zz|A'>");
      TArray2 Qxx_mA = q_xx_mA("m,A'") - (q_zz_mA("m,A'")  + q_yy_mA("m,A'") ) * 0.5;
      TArray2 Qyy_mA = q_yy_mA("m,A'") - (q_xx_mA("m,A'")  + q_zz_mA("m,A'") ) * 0.5;
      TArray2 Qzz_mA = q_zz_mA("m,A'") - (q_xx_mA("m,A'")  + q_yy_mA("m,A'") ) * 0.5;
      TArray2 Qxy_mA = _2("<m|q_xy|A'>") * 1.5;
      TArray2 Qxz_mA = _2("<m|q_xz|A'>") * 1.5;
      TArray2 Qyz_mA = _2("<m|q_yz|A'>") * 1.5;

      const double q_xx_e2 = - dot(Qxx_mn("m,n"), D_e2_mn("m,n"))
                             + dot(Qxx_AB("A',B'"), D_e2_AB("A',B'"))
                             + dot(Qxx_mA("m,A'"), TmA("m,A'")) * 2.0;
      const double q_yy_e2 = - dot(Qyy_mn("m,n"), D_e2_mn("m,n"))
                             + dot(Qyy_AB("A',B'"), D_e2_AB("A',B'"))
                             + dot(Qyy_mA("m,A'"), TmA("m,A'")) * 2.0;
      const double q_zz_e2 = - dot(Qzz_mn("m,n"), D_e2_mn("m,n"))
                             + dot(Qzz_AB("A',B'"), D_e2_AB("A',B'"))
                             + dot(Qzz_mA("m,A'"), TmA("m,A'")) * 2.0;
      const double q_xy_e2 = - dot(Qxy_mn("m,n"), D_e2_mn("m,n"))
                             + dot(Qxy_AB("A',B'"), D_e2_AB("A',B'"))
                             + dot(Qxy_mA("m,A'"), TmA("m,A'")) * 2.0;
      const double q_xz_e2 = - dot(Qxz_mn("m,n"), D_e2_mn("m,n"))
                             + dot(Qxz_AB("A',B'"), D_e2_AB("A',B'"))
                             + dot(Qxz_mA("m,A'"), TmA("m,A'")) * 2.0;
      const double q_yz_e2 = - dot(Qyz_mn("m,n"), D_e2_mn("m,n"))
                             + dot(Qyz_AB("A',B'"), D_e2_AB("A',B'"))
                             + dot(Qyz_mA("m,A'"), TmA("m,A'")) * 2.0;

      const double q_xx_e2or = dot(Qxx_am("a,m"), Dbn_E2("a,m"));
      const double q_yy_e2or = dot(Qyy_am("a,m"), Dbn_E2("a,m"));
      const double q_zz_e2or = dot(Qzz_am("a,m"), Dbn_E2("a,m"));
      const double q_xy_e2or = dot(Qxy_am("a,m"), Dbn_E2("a,m"));
      const double q_xz_e2or = dot(Qxz_am("a,m"), Dbn_E2("a,m"));
      const double q_yz_e2or = dot(Qyz_am("a,m"), Dbn_E2("a,m"));

      std::cout << std::endl
                << "traceless quadrupole moment (E2)" << std::endl
                << "q_xx (E2) = " << scprintf("%12.10f", - q_xx_e2 * 2.0)
                << "  q_yy (E2) = " << scprintf("%12.10f", - q_yy_e2 * 2.0)
                << "  q_zz (E2) = " << scprintf("%12.10f", - q_zz_e2 * 2.0)
                << std::endl
                << "q_xy (E2) = " << scprintf("%12.10f", - q_xy_e2 * 2.0)
                << "  q_xz (E2) = " << scprintf("%12.10f", - q_xz_e2 * 2.0)
                << "  q_yz (E2) = " << scprintf("%12.10f", - q_yz_e2 * 2.0)
                << std::endl;
      std::cout << std::endl
                << "traceless quadrupole moment (E2 orbital response)" << std::endl
                << "q_xx (E2 or) = " << scprintf("%12.10f", - q_xx_e2or * 2.0)
                << "  q_yy (E2 or) = " << scprintf("%12.10f", - q_yy_e2or * 2.0)
                << "  q_zz (E2 or) = " << scprintf("%12.10f", - q_zz_e2or * 2.0)
                << std::endl
                << "q_xy (E2 or) = " << scprintf("%12.10f", - q_xy_e2or * 2.0)
                << "  q_xz (E2 or) = " << scprintf("%12.10f", - q_xz_e2or * 2.0)
                << "  q_yz (E2 or) = " << scprintf("%12.10f", - q_yz_e2or * 2.0)
                << std::endl << std::endl;
    }
    world_.gop.fence();

    // CABS singles contribution to
    // electric field gradient near the 1st nucleus
    if (compute_EFG) {
      TArray2 v_xx_mA = _2("<m|ddphi_xx|A'>");
      TArray2 v_yy_mA = _2("<m|ddphi_yy|A'>");
      TArray2 v_zz_mA = _2("<m|ddphi_zz|A'>");
      TArray2 v_xy_mA = _2("<m|ddphi_xy|A'>");
      TArray2 v_xz_mA = _2("<m|ddphi_xz|A'>");
      TArray2 v_yz_mA = _2("<m|ddphi_yz|A'>");

      TArray2 v_xx_AB = _2("<A'|ddphi_xx|B'>");
      TArray2 v_yy_AB = _2("<A'|ddphi_yy|B'>");
      TArray2 v_zz_AB = _2("<A'|ddphi_zz|B'>");
      TArray2 v_xy_AB = _2("<A'|ddphi_xy|B'>");
      TArray2 v_xz_AB = _2("<A'|ddphi_xz|B'>");
      TArray2 v_yz_AB = _2("<A'|ddphi_yz|B'>");

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

      std::cout << std::endl
                << "electric gradient (E2)" << std::endl
                << "v_xx (E2) = " << scprintf("%12.10f", v_xx_e2 * 2.0)
                << "  v_yy (E2) = " << scprintf("%12.10f", v_yy_e2 * 2.0)
                << "  v_zz (E2) = " << scprintf("%12.10f", v_zz_e2 * 2.0)
                << std::endl
                << "v_xy (E2) = " << scprintf("%12.10f", v_xy_e2 * 2.0)
                << "  v_xz (E2) = " << scprintf("%12.10f", v_xz_e2 * 2.0)
                << "  v_yz (E2) = " << scprintf("%12.10f", v_yz_e2 * 2.0)
                << std::endl;
      std::cout << std::endl
                << "electric gradient (E2 orbital response)" << std::endl
                << "v_xx (E2 or) = " << scprintf("%12.10f", v_xx_e2or * 2.0)
                << "  v_yy (E2 or) = " << scprintf("%12.10f", v_yy_e2or * 2.0)
                << "  v_zz (E2 or) = " << scprintf("%12.10f", v_zz_e2or * 2.0)
                << std::endl
                << "v_xy (E2 or) = " << scprintf("%12.10f", v_xy_e2or * 2.0)
                << "  v_xz (E2 or) = " << scprintf("%12.10f", v_xz_e2or * 2.0)
                << "  v_yz (E2 or) = " << scprintf("%12.10f", v_yz_e2or * 2.0)
                << std::endl << std::endl;
    }
    world_.gop.fence();

    }
    world_.gop.fence();
#endif

    // compute integrals needed for MP2 and following calculations of properties
    // dipole integrals
    TArray2 mu_z_ij = (compute_dipole ? _2("<i|mu_z|j>") : TArray2());
    TArray2 mu_z_ab = (compute_dipole ? _2("<a|mu_z|b>") : TArray2());

    // quadrupole integrasl
    TArray2 q_xx_ij = (compute_quadrupole ? _2("<i|q_xx|j>") : TArray2());
    TArray2 q_yy_ij = (compute_quadrupole ? _2("<i|q_yy|j>") : TArray2());
    TArray2 q_zz_ij = (compute_quadrupole ? _2("<i|q_zz|j>") : TArray2());
    TArray2 Qxx_ij = (compute_quadrupole ?
                      q_xx_ij("i,j") - (q_zz_ij("i,j") + q_yy_ij("i,j")) * 0.5 : TArray2());
    TArray2 Qyy_ij = (compute_quadrupole ?
                      q_yy_ij("i,j") - (q_zz_ij("i,j") + q_xx_ij("i,j")) * 0.5 : TArray2());
    TArray2 Qzz_ij = (compute_quadrupole ?
                      q_zz_ij("i,j") - (q_xx_ij("i,j") + q_yy_ij("i,j")) * 0.5 : TArray2());
    TArray2 Qxy_ij = (compute_quadrupole ? _2("<i|q_xy|j>") * 1.5 : TArray2());
    TArray2 Qxz_ij = (compute_quadrupole ? _2("<i|q_xz|j>") * 1.5 : TArray2());
    TArray2 Qyz_ij = (compute_quadrupole ? _2("<i|q_yz|j>") * 1.5 : TArray2());

    TArray2 q_xx_ab = (compute_quadrupole ? _2("<a|q_xx|b>") : TArray2());
    TArray2 q_yy_ab = (compute_quadrupole ? _2("<a|q_yy|b>") : TArray2());
    TArray2 q_zz_ab = (compute_quadrupole ? _2("<a|q_zz|b>") : TArray2());
    TArray2 Qxx_ab = (compute_quadrupole ?
                      q_xx_ab("a,b") - (q_zz_ab("a,b") + q_yy_ab("a,b")) * 0.5 : TArray2());
    TArray2 Qyy_ab = (compute_quadrupole ?
                      q_yy_ab("a,b") - (q_zz_ab("a,b") + q_xx_ab("a,b")) * 0.5 : TArray2());
    TArray2 Qzz_ab = (compute_quadrupole ?
                      q_zz_ab("a,b") - (q_xx_ab("a,b") + q_yy_ab("a,b")) * 0.5 : TArray2());
    TArray2 Qxy_ab = (compute_quadrupole ? _2("<a|q_xy|b>") * 1.5 : TArray2());
    TArray2 Qxz_ab = (compute_quadrupole ? _2("<a|q_xz|b>") * 1.5 : TArray2());
    TArray2 Qyz_ab = (compute_quadrupole ? _2("<a|q_yz|b>") * 1.5 : TArray2());

    // electric field gradient integrals
    TArray2 v_xx_ij = (compute_EFG ? _2("<i|ddphi_xx|j>") : TArray2());
    TArray2 v_yy_ij = (compute_EFG ? _2("<i|ddphi_yy|j>") : TArray2());
    TArray2 v_zz_ij = (compute_EFG ? _2("<i|ddphi_zz|j>") : TArray2());
    TArray2 v_xy_ij = (compute_EFG ? _2("<i|ddphi_xy|j>") : TArray2());
    TArray2 v_xz_ij = (compute_EFG ? _2("<i|ddphi_xz|j>") : TArray2());
    TArray2 v_yz_ij = (compute_EFG ? _2("<i|ddphi_yz|j>") : TArray2());

    TArray2 v_xx_ab = (compute_EFG ? _2("<a|ddphi_xx|b>") : TArray2());
    TArray2 v_yy_ab = (compute_EFG ? _2("<a|ddphi_yy|b>") : TArray2());
    TArray2 v_zz_ab = (compute_EFG ? _2("<a|ddphi_zz|b>") : TArray2());
    TArray2 v_xy_ab = (compute_EFG ? _2("<a|ddphi_xy|b>") : TArray2());
    TArray2 v_xz_ab = (compute_EFG ? _2("<a|ddphi_xz|b>") : TArray2());
    TArray2 v_yz_ab = (compute_EFG ? _2("<a|ddphi_yz|b>") : TArray2());

    // compute MP2 and its orbital response contribution

    // compute Delta_ijab = - 1 / (- <i|F|i> - <j|F|j> + <a|F|a> + <b|F|b>)
    // which is needed for MP2 amplitudes
    TArray2 mFij = - _2("<i|F|j>");
    TArray4 g_ijab = _4("<i j|g|a b>");

    typedef detail::diag_precond4<double> pc4eval_type;
    typedef TA::Array<T, 4, LazyTensor<T, 4, pc4eval_type > > TArray4d;
    TArray4d Delta_ijab(g_ijab.get_world(), g_ijab.trange());
    pc4eval_type Delta_ijab_gen(TA::array_to_eigen(mFij), TA::array_to_eigen(mFij),
                                TA::array_to_eigen(mFab),TA::array_to_eigen(mFab));
    // construct local tiles
    for(auto t = Delta_ijab.trange().tiles().begin();
        t != Delta_ijab.trange().tiles().end(); ++t)
      if (Delta_ijab.is_local(*t)) {
        std::array<std::size_t, 4> index;
        std::copy(t->begin(), t->end(), index.begin());
        madness::Future < typename TArray4d::value_type >
          tile((LazyTensor<T, 4, pc4eval_type >(&Delta_ijab, index, &Delta_ijab_gen)
              ));

        // Insert the tile into the array
        Delta_ijab.set(*t, tile);
      }

    // MP2 amplitues:
    TArray4 T2_ijab = g_ijab("i,j,a,b") * Delta_ijab("i,j,a,b");

#if 1
    // MP2 density
    TArray2 D_mp2_ij = (2.0 * T2_ijab("i,k,a,b") - T2_ijab("k,i,a,b"))
                          * T2_ijab("j,k,a,b");
    TArray2 D_mp2_ab = (2.0 * T2_ijab("i,j,a,c") - T2_ijab("i,j,c,a"))
                          * T2_ijab("i,j,b,c");
    // MP2 orbital response
    TArray2 X_mp2 = Xam_mp2(T2_ijab, D_mp2_ij, D_mp2_ab);
    TArray2 Dbn_mp2(X_mp2.get_world(), X_mp2.trange());
    auto resnorm_mp2 = cg_solver2(Orbital_relaxation_Abnam,
                                  X_mp2,
                                  Dbn_mp2,
                                  preconditioner,
                                  1e-10);
    // MP2 dipole
    if (compute_dipole) {
      // MP2 density contribution to dipole
      const double mu_z_mp2 = - dot(mu_z_ij("i,j"), D_mp2_ij("i,j"))
                            + dot(mu_z_ab("a,b"), D_mp2_ab("a,b"));
      std::cout << std::endl << "mu_z (MP2) = "
                << scprintf("%12.10f", - mu_z_mp2 * 2.0) << std::endl;

      // MP2 orbital response contribution to dipole
      const double mu_z_mp2or = dot(mu_z_am("a,m"), Dbn_mp2("a,m"));
      std::cout << std::endl
                << "mu_z (MP2 orbital response) = "<< scprintf("%12.10f", - mu_z_mp2or * 2.0)
                << std::endl << std::endl;
    }

    // MP2 contribution to quadrupoles
    if (compute_quadrupole) {
      const double q_xx_mp2 = - dot(Qxx_ij("i,j"), D_mp2_ij("i,j"))
                              + dot(Qxx_ab("a,b"), D_mp2_ab("a,b"));
      const double q_yy_mp2 = - dot(Qyy_ij("i,j"), D_mp2_ij("i,j"))
                              + dot(Qyy_ab("a,b"), D_mp2_ab("a,b"));
      const double q_zz_mp2 = - dot(Qzz_ij("i,j"), D_mp2_ij("i,j"))
                              + dot(Qzz_ab("a,b"), D_mp2_ab("a,b"));
      const double q_xy_mp2 = - dot(Qxy_ij("i,j"), D_mp2_ij("i,j"))
                              + dot(Qxy_ab("a,b"), D_mp2_ab("a,b"));
      const double q_xz_mp2 = - dot(Qxz_ij("i,j"), D_mp2_ij("i,j"))
                              + dot(Qxz_ab("a,b"), D_mp2_ab("a,b"));
      const double q_yz_mp2 = - dot(Qyz_ij("i,j"), D_mp2_ij("i,j"))
                              + dot(Qyz_ab("a,b"), D_mp2_ab("a,b"));

      const double q_xx_mp2or = dot(Qxx_am("a,m"), Dbn_mp2("a,m"));
      const double q_yy_mp2or = dot(Qyy_am("a,m"), Dbn_mp2("a,m"));
      const double q_zz_mp2or = dot(Qzz_am("a,m"), Dbn_mp2("a,m"));
      const double q_xy_mp2or = dot(Qxy_am("a,m"), Dbn_mp2("a,m"));
      const double q_xz_mp2or = dot(Qxz_am("a,m"), Dbn_mp2("a,m"));
      const double q_yz_mp2or = dot(Qyz_am("a,m"), Dbn_mp2("a,m"));

      std::cout << std::endl
                << "traceless quadrupole moment (MP2)" << std::endl
                << "q_xx (MP2) = " << scprintf("%12.10f", - q_xx_mp2 * 2.0)
                << "  q_yy (MP2) = " << scprintf("%12.10f", - q_yy_mp2 * 2.0)
                << "  q_zz (MP2) = " << scprintf("%12.10f", - q_zz_mp2 * 2.0)
                << std::endl
                << "q_xy (MP2) = " << scprintf("%12.10f", - q_xy_mp2 * 2.0)
                << "  q_xz (MP2) = " << scprintf("%12.10f", - q_xz_mp2 * 2.0)
                << "  q_yz (MP2) = " << scprintf("%12.10f", - q_yz_mp2 * 2.0)
                << std::endl;
      std::cout << std::endl
                << "traceless quadrupole moment (MP2 orbital response)" << std::endl
                << "q_xx (MP2 or) = " << scprintf("%12.10f", - q_xx_mp2or * 2.0)
                << "  q_yy (MP2 or) = " << scprintf("%12.10f", - q_yy_mp2or * 2.0)
                << "  q_zz (MP2 or) = " << scprintf("%12.10f", - q_zz_mp2or * 2.0)
                << std::endl
                << "q_xy (MP2 or) = " << scprintf("%12.10f", - q_xy_mp2or * 2.0)
                << "  q_xz (MP2 or) = " << scprintf("%12.10f", - q_xz_mp2or * 2.0)
                << "  q_yz (MP2 or) = " << scprintf("%12.10f", - q_yz_mp2or * 2.0)
                << std::endl << std::endl;
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

      std::cout << std::endl
                << "electric gradient (MP2)" << std::endl
                << "v_xx (MP2) = " << scprintf("%12.10f", - v_xx_mp2 * 2.0)
                << "  v_yy (MP2) = " << scprintf("%12.10f", - v_yy_mp2 * 2.0)
                << "  v_zz (MP2) = " << scprintf("%12.10f", - v_zz_mp2 * 2.0)
                << std::endl
                << "v_xy (MP2) = " << scprintf("%12.10f", - v_xy_mp2 * 2.0)
                << "  v_xz (MP2) = " << scprintf("%12.10f", - v_xz_mp2 * 2.0)
                << "  v_yz (MP2) = " << scprintf("%12.10f", - v_yz_mp2 * 2.0)
                << std::endl;
      std::cout << std::endl
                << "electric gradient (MP2 orbital response)" << std::endl
                << "v_xx (MP2 or) = " << scprintf("%12.10f", - v_xx_mp2or * 2.0)
                << "  v_yy (MP2 or) = " << scprintf("%12.10f", - v_yy_mp2or * 2.0)
                << "  v_zz (MP2 or) = " << scprintf("%12.10f", - v_zz_mp2or * 2.0)
                << std::endl
                << "v_xy (MP2 or) = " << scprintf("%12.10f", - v_xy_mp2or * 2.0)
                << "  v_xz (MP2 or) = " << scprintf("%12.10f", - v_xz_mp2or * 2.0)
                << "  v_yz (MP2 or) = " << scprintf("%12.10f", - v_yz_mp2or * 2.0)
                << std::endl << std::endl;
    }
#endif

    // compute integrals needed for F12 related calculations of properties
    TArray2 mu_z_apb = (compute_dipole ? _2("<a'|mu_z|b>") : TArray2());

    TArray2 q_xx_apb = (compute_quadrupole ? _2("<a'|q_xx|b>") : TArray2());
    TArray2 q_yy_apb = (compute_quadrupole ? _2("<a'|q_yy|b>") : TArray2());
    TArray2 q_zz_apb = (compute_quadrupole ? _2("<a'|q_zz|b>") : TArray2());

    TArray2 Qxx_apb = (compute_quadrupole ?
                       q_xx_apb("a',b") - (q_zz_apb("a',b") + q_yy_apb("a',b")) * 0.5 : TArray2());
    TArray2 Qyy_apb = (compute_quadrupole ?
                       q_yy_apb("a',b") - (q_zz_apb("a',b") + q_xx_apb("a',b")) * 0.5 : TArray2());
    TArray2 Qzz_apb = (compute_quadrupole ?
                       q_zz_apb("a',b") - (q_xx_apb("a',b") + q_yy_apb("a',b")) * 0.5 : TArray2());
    TArray2 Qxy_apb = (compute_quadrupole ? _2("<a'|q_xy|b>") * 1.5 : TArray2());
    TArray2 Qxz_apb = (compute_quadrupole ? _2("<a'|q_xz|b>") * 1.5 : TArray2());
    TArray2 Qyz_apb = (compute_quadrupole ? _2("<a'|q_yz|b>") * 1.5 : TArray2());

    TArray2 v_xx_apb = (compute_EFG ? _2("<a'|ddphi_xx|b>") : TArray2());
    TArray2 v_yy_apb = (compute_EFG ? _2("<a'|ddphi_yy|b>") : TArray2());
    TArray2 v_zz_apb = (compute_EFG ? _2("<a'|ddphi_zz|b>") : TArray2());
    TArray2 v_xy_apb = (compute_EFG ? _2("<a'|ddphi_xy|b>") : TArray2());
    TArray2 v_xz_apb = (compute_EFG ? _2("<a'|ddphi_xz|b>") : TArray2());
    TArray2 v_yz_apb = (compute_EFG ? _2("<a'|ddphi_yz|b>") : TArray2());

    // singlet and triplet coefficients for F12 and coupling terms
    const double C_0 = 1.0 / 2.0;
    const double C_1 = 1.0 / 4.0;
    // compute coefficients needed in the F12 and coupling calculations
    const double R_C1 = (0.5 * C_0 + 1.5 * C_1);
    const double R_C2 = (0.5 * C_0 - 1.5 * C_1);
    const double RR_C1 = 0.5 * C_0 * C_0 + 1.5 * C_1 * C_1;
    const double RR_C2 = 0.5 * C_0 * C_0 - 1.5 * C_1 * C_1;

    // MP2 F12 coupling part
#if 1
    {
    TArray4 C_ijab = _4("<i j|r|a_F(a') b>") + _4("<i j|r|a b_F(a')>");
    TArray4 A_ijab = C_ijab("i,j,a,b") * Delta_ijab("i,j,a,b");

    // Coupling and F12 part of MP2F12 density
    TArray2 D_mp2f12_ij =  (R_C1 * A_ijab("i,k,a,b") + R_C2 * A_ijab("k,i,a,b"))
                           * T2_ijab("j,k,a,b")
                         // A_jk^ab T_ab^ik
                         + (R_C1 * A_ijab("j,k,a,b") + R_C2 * A_ijab("k,j,a,b"))
                           * T2_ijab("i,k,a,b")
                         // A^ik_ab A^ab_jk
                         + (RR_C1 * A_ijab("i,k,a,b") + RR_C2 * A_ijab("k,i,a,b"))
                           * A_ijab("j,k,a,b");
    TArray2 D_mp2f12_ab =  (R_C1 * T2_ijab("i,j,a,c") + R_C2 * T2_ijab("i,j,c,a"))
                           * A_ijab("i,j,b,c")
                         // A^ac_ij T^ij_bc
                         + (R_C1 * T2_ijab("i,j,b,c") + R_C2 * T2_ijab("i,j,c,b"))
                           * A_ijab("i,j,a,c")
                         //  A_ac^ij A_ij^bc
                         + (RR_C1 * A_ijab("i,j,a,c") + RR_C2 * A_ijab("i,j,c,a"))
                           * A_ijab("i,j,b,c");

    // 1/2 R^kl_a'c \tilde{T}^bc_kl
    TArray2 RT_apb = _4("<a' c|r|k l>")
                     * ( R_C1 * T2_ijab("k,l,b,c") + R_C2 * T2_ijab("k,l,c,b")
                       + RR_C1 * A_ijab("k,l,b,c") + RR_C2 * A_ijab("k,l,c,b")
                       );

    // MP2 F12 coupling contribution to orbital response
    TArray2 Xmp2f12_contri = Xam_Cmp2f12(C_0, C_1,T2_ijab, A_ijab,
                                         D_mp2f12_ij, D_mp2f12_ab, RT_apb);
    TArray2 Dbn_mp2f12(Xmp2f12_contri.get_world(), Xmp2f12_contri.trange());
    auto resnorm_mp2f12 = cg_solver2(Orbital_relaxation_Abnam,
                                     Xmp2f12_contri,
                                     Dbn_mp2f12,
                                     preconditioner,
                                     1e-10);

    // MP2-F12 coupling contribution to dipole
    if (compute_dipole) {
      // MP2 F12 coupling density contribution to dipole
      const double mu_z_mp2f12 = - dot(mu_z_ij("i,j"), D_mp2f12_ij("i,j"))
                                 + dot(mu_z_ab("a,b"), D_mp2f12_ab("a,b"))
                                 + dot(mu_z_apb("a',b"), RT_apb("a',b")) * 2.0;
      std::cout << std::endl << "mu_z (MP2F12 coupling) = "
                << scprintf("%12.10f", - mu_z_mp2f12 * 2.0) << std::endl;

      // MP2-F12 coupling orbital response contribution to dipole
      const double mu_z_Xam_mp2f12 = dot(mu_z_am("a,m"), Dbn_mp2f12("a,m"));
      std::cout << std::endl
                << "mu_z (MP2-F12 coupling orbital response) = " << scprintf("%12.10f", - mu_z_Xam_mp2f12 * 4.0)
                << std::endl << std::endl;;
    }

    // MP2-F12 coupling contribution to quadrupoles
    if (compute_quadrupole) {
      const double q_xx_mp2f12C = - dot(Qxx_ij("i,j"), D_mp2f12_ij("i,j"))
                                  + dot(Qxx_ab("a,b"), D_mp2f12_ab("a,b"))
                                  + dot(Qxx_apb("a',b"), RT_apb("a',b")) * 2.0;
      const double q_yy_mp2f12C = - dot(Qyy_ij("i,j"), D_mp2f12_ij("i,j"))
                                  + dot(Qyy_ab("a,b"), D_mp2f12_ab("a,b"))
                                  + dot(Qyy_apb("a',b"), RT_apb("a',b")) * 2.0;
      const double q_zz_mp2f12C = - dot(Qzz_ij("i,j"), D_mp2f12_ij("i,j"))
                                  + dot(Qzz_ab("a,b"), D_mp2f12_ab("a,b"))
                                  + dot(Qzz_apb("a',b"), RT_apb("a',b")) * 2.0;
      const double q_xy_mp2f12C = - dot(Qxy_ij("i,j"), D_mp2f12_ij("i,j"))
                                  + dot(Qxy_ab("a,b"), D_mp2f12_ab("a,b"))
                                  + dot(Qxy_apb("a',b"), RT_apb("a',b")) * 2.0;
      const double q_xz_mp2f12C = - dot(Qxz_ij("i,j"), D_mp2f12_ij("i,j"))
                                  + dot(Qxz_ab("a,b"), D_mp2f12_ab("a,b"))
                                  + dot(Qxz_apb("a',b"), RT_apb("a',b")) * 2.0;
      const double q_yz_mp2f12C = - dot(Qyz_ij("i,j"), D_mp2f12_ij("i,j"))
                                  + dot(Qyz_ab("a,b"), D_mp2f12_ab("a,b"))
                                  + dot(Qyz_apb("a',b"), RT_apb("a',b")) * 2.0;

      const double q_xx_mp2f12Cor = dot(Qxx_am("a,m"), Dbn_mp2f12("a,m"));
      const double q_yy_mp2f12Cor = dot(Qyy_am("a,m"), Dbn_mp2f12("a,m"));
      const double q_zz_mp2f12Cor = dot(Qzz_am("a,m"), Dbn_mp2f12("a,m"));
      const double q_xy_mp2f12Cor = dot(Qxy_am("a,m"), Dbn_mp2f12("a,m"));
      const double q_xz_mp2f12Cor = dot(Qxz_am("a,m"), Dbn_mp2f12("a,m"));
      const double q_yz_mp2f12Cor = dot(Qyz_am("a,m"), Dbn_mp2f12("a,m"));

      std::cout << std::endl
                << "traceless quadrupole moment (MP2F12 coupling)" << std::endl
                << "q_xx (MP2F12 C) = " << scprintf("%12.10f", - q_xx_mp2f12C * 2.0)
                << "  q_yy (MP2F12 C) = " << scprintf("%12.10f", - q_yy_mp2f12C * 2.0)
                << "  q_zz (MP2F12 C) = " << scprintf("%12.10f", - q_zz_mp2f12C * 2.0)
                << std::endl
                << "q_xy (MP2F12 C) = " << scprintf("%12.10f", - q_xy_mp2f12C * 2.0)
                << "  q_xz (MP2F12 C) = " << scprintf("%12.10f", - q_xz_mp2f12C * 2.0)
                << "  q_yz (MP2F12 C) = " << scprintf("%12.10f", - q_yz_mp2f12C * 2.0)
                << std::endl;
      std::cout << std::endl
                << "traceless quadrupole moment (MP2F12 coupling orbital response)" << std::endl
                << "q_xx (MP2F12 C or) = " << scprintf("%12.10f", - q_xx_mp2f12Cor * 2.0)
                << "  q_yy (MP2F12 C or) = " << scprintf("%12.10f", - q_yy_mp2f12Cor * 2.0)
                << "  q_zz (MP2F12 C or) = " << scprintf("%12.10f", - q_zz_mp2f12Cor * 2.0)
                << std::endl
                << "q_xy (MP2F12 C or) = " << scprintf("%12.10f", - q_xy_mp2f12Cor * 2.0)
                << "  q_xz (MP2F12 C or) = " << scprintf("%12.10f", - q_xz_mp2f12Cor * 2.0)
                << "  q_yz (MP2F12 C or) = " << scprintf("%12.10f", - q_yz_mp2f12Cor * 2.0)
                << std::endl << std::endl;
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

      std::cout << std::endl
                << "electric gradient (MP2F12 coupling)" << std::endl
                << "v_xx (MP2F12 C) = " << scprintf("%12.10f", v_xx_mp2f12C * 2.0)
                << "  v_yy (MP2F12 C) = " << scprintf("%12.10f", v_yy_mp2f12C * 2.0)
                << "  v_zz (MP2F12 C) = " << scprintf("%12.10f", v_zz_mp2f12C * 2.0)
                << std::endl
                << "v_xy (MP2F12 C) = " << scprintf("%12.10f", v_xy_mp2f12C * 2.0)
                << "  v_xz (MP2F12 C) = " << scprintf("%12.10f", v_xz_mp2f12C * 2.0)
                << "  v_yz (MP2F12 C) = " << scprintf("%12.10f", v_yz_mp2f12C * 2.0)
                << std::endl;
      std::cout << std::endl
                << "electric gradient (MP2F12 coupling orbital response)" << std::endl
                << "v_xx (MP2F12 C or) = " << scprintf("%12.10f", v_xx_mp2f12Cor * 2.0)
                << "  v_yy (MP2F12 C or) = " << scprintf("%12.10f", v_yy_mp2f12Cor * 2.0)
                << "  v_zz (MP2F12 C or) = " << scprintf("%12.10f", v_zz_mp2f12Cor * 2.0)
                << std::endl
                << "v_xy (MP2F12 C or) = " << scprintf("%12.10f", v_xy_mp2f12Cor * 2.0)
                << "  v_xz (MP2F12 C or) = " << scprintf("%12.10f", v_xz_mp2f12Cor * 2.0)
                << "  v_yz (MP2F12 C or) = " << scprintf("%12.10f", v_yz_mp2f12Cor * 2.0)
                << std::endl << std::endl;
    }

    }
    world_.gop.fence();
#endif

    // F12 contribution
#if 1
    // F12 density from X and B terms
    TArray4 r2_ijkl = _4("<i j|r2|k l>");
    TArray4 r_ijpq = _4("<i j|r|p q>");
    TArray4 r_ijapn = _4("<i j|r|a' n>");
    TArray4 r_ijnap = _4("<i j|r|n a'>");

    // Dij = 1/2 R^ik_A'B' R^A'B'_kl (A': all virtual)
    TArray2 D_f12_ij =  (RR_C1 * r2_ijkl("i,k,j,l") + RR_C2 * r2_ijkl("k,i,j,l"))
                        * _2("<k|I|l>")
                      - (RR_C1 * r_ijpq("i,k,p,q") + RR_C2 * r_ijpq("k,i,p,q"))
                        * r_ijpq("j,k,p,q")
                      - (RR_C1 * r_ijapn("i,k,a',n") + RR_C2 * r_ijapn("k,i,a',n"))
                        * r_ijapn("j,k,a',n")
                      - (RR_C1 * r_ijnap("i,k,n,a'") + RR_C2 * r_ijnap("k,i,n,a'"))
                        * r_ijnap("j,k,n,a'");

    // DA'B' = 1/2 R^A'C'_kl R^kl_B'C' (A': all virtual)
    TArray4 r_acpkl = _4("<a c'|r|k l>");
    TArray2 D_f12_ab = (RR_C1 * r_acpkl("a,c',k,l") + RR_C2 * r_acpkl("a,c',l,k"))
                       * r_acpkl("b,c',k,l");

    TArray4 r_apcpkl = _4("<a' c'|r|k l>");
    TArray2 D_f12_apbp = (RR_C1 * r_acpkl("c,a',l,k") + RR_C2 * r_acpkl("c,a',k,l"))
                         * r_acpkl("c,b',l,k")
                       + (RR_C1 * r_apcpkl("a',c',k,l") + RR_C2 * r_apcpkl("a',c',l,k"))
                         * r_apcpkl("b',c',k,l");

    TArray2 D_f12_apb = (RR_C1 * r_apcpkl("a',c',k,l") + RR_C2 * r_apcpkl("a',c',l,k"))
                        * r_acpkl("b,c',k,l");

    // F12 orbital response contribution to dipole
    // X and B density contribution to Xam
    TArray4 g_abmc = _4("<a b|g|m c>");
    TArray2 gdf12_am =   (2.0 * _4("<a k|g|m l>") - _4("<a k|g|l m>"))
                         * D_f12_ij("k,l")
                       //
                       - (2.0 * g_abmc("a,b,m,c") - g_abmc("b,a,m,c"))
                         * D_f12_ab("b,c")
                       //
                       - (2.0 * _4("<a b'|g|m c'>") - _4("<a b'|g|c' m>"))
                         * D_f12_apbp("b',c'")
                       //
                       - (2.0 * _4("<a b'|g|m c>") - _4("<a b'|g|c m>"))
                         * D_f12_apb("b',c")
                       - (2.0 * _4("<a b|g|m c'>") - _4("<a b|g|c' m>"))
                         * D_f12_apb("c',b");

    // V contribution to F12 Xam
    TArray2 Xam_Vcontri = Xam_V(C_0,C_1);

    // X contribution to F12 Xam
    TArray2 Xam_Xcontri = Xam_X(C_0,C_1);

    // B contribution to F12 Xam
    TArray2 Xam_Bcontri = Xam_B(C_0,C_1);

    // CC F12 coupling contribution to Xai
#if 0
    // CC CT2
    TArray4 T2_ijab = _4("<i j|T2|a b>");
    // 1/2 R^kl_a'c T2^bc_kl
    TArray2 RT2_apb = _4("<a' c|r|k l>")
                      * (V_C1 * T2_ijab("k,l,b,c") + V_C2 * T2_ijab("k,l,c,b"));
    // test for RT2_apb
//    TArray2 mu_z_apb = _2("<a'|mu_z|b>");
//    const double mu_z_RT2 = dot(mu_z_apb("a',b"), RT2_apb("a',b"));
//    std::cout << std::endl << "***  "
//              << "mu_z (RT2) = " << - mu_z_RT2 * 4.0
//              << "  ***"<< std::endl;

    TArray2 CT2_ai =
                   // 1/2 F^a'_i R^kl_a'b T2^ab_kl
                 - (V_C1 * T2_ijab("k,l,a,b") + V_C2 * T2_ijab("k,l,b,a"))
                   * _4("<k l|r|i_F(a') b>")
                   //   g^aa'_ib 1/2 R^kl_a'c T2^bc_kl
                   // & g^ab_ia' 1/2 R^kl_a'c T2^bc_kl
                 + ( 2.0 * _4("<a a'|g|i b>") - _4("<a a'|g|b i>")
                   + 2.0 * _4("<a b|g|i a'>") - _4("<a b|g|a' i>")
                    ) * RT2_apb("a',b")

                    // 1/2 F^a'_b R^kl_ia' T2^ab_kl
                  - (V_C1 * T2_ijab("k,l,a,b") + V_C2 * T2_ijab("k,l,b,a"))
                    * _4("<k l|r|i b_F(a')>")

                    // F^a'_b R^al_a'c T2^bc_il
                  + _4("<a l|r|b_F(a') c>")
                    * (V_C1 * T2_ijab("i,l,b,c") + V_C2 * T2_ijab("l,i,b,c"))
                  + _4("<a l|r|c b_F(a')>")
                    * (V_C1 * T2_ijab("i,l,c,b") + V_C2 * T2_ijab("l,i,c,b"))
                 ;

    // VT1 & VT2 coupling contribution to F12 Xai

    const char* a = "a";
    const char* c = "c";
    const char* d = "d";

    const char* i = "i";
    const char* j = "j";
    const char* k = "k";
    const char* l = "l";
    // can not use p q n a' for Vpq_rs
    TArray4 V_alic = Vpq_rs(a,l,i,c);
    TArray4 V_ilac = Vpq_rs(i,l,a,c);
    // can not use k l n a' for Vrk_sk
    TArray2 V_ai = Vrk_sk(a,c);
    TArray2 V_ji = Vrk_sk(j,i);
    TArray2 Tia = _2("<i|T1|a>");
    TArray2 VT1_ai =
                    //   1/2 R^al_A'B' g^A'B'_ic t1^c_l
                    // + 1/2 R^il_A'B' g^A'B'_ac t1^c_l
                   ( V_alic("a,l,i,c") + V_ilac("i,l,a,c")
//                   // test for V_alic("a,l,i,c")
//                     V_C1 * _4("<a l|gr|i c>") + V_C2 * _4("<l a|gr|i c>")
//                   - (V_C1 * _4("<a l|r|p q>") + V_C2 * _4("<l a|r|p q>"))
//                     * _4("<p q|g|i c>")
//                   - (V_C1 * _4("<a l|r|a' n>") + V_C2 * _4("<l a|r|a' n>"))
//                     * _4("<a' n|g|i c>")
//                   - (V_C1 * _4("<a l|r|n a'>") + V_C2 * _4("<l a|r|n a'>"))
//                     * _4("<n a'|g|i c>")
//                   // test for V_ilac("i,l,a,c")
//                    V_C1 * _4("<a c|gr|i l>") + V_C2 * _4("<a c|gr|l i>")
//                   - _4("<a c|g|p q>")
//                     * (V_C1 * _4("<p q|r|i l>") + V_C2 * _4("<p q|r|l i>"))
//                   - _4("<a c|g|a' n>")
//                     * (V_C1 * _4("<a' n|r|i l>") + V_C2 * _4("<a' n|r|l i>"))
//                   - _4("<a c|g|n a'>")
//                     * (V_C1 * _4("<n a'|r|i l>") + V_C2 * _4("<n a'|r|l i>"))
                   ) * Tia("l,c")

                   // 1/2 R^ak_A'B' g^A'B'_ck t1^c_i
                   + V_ai("a,c") * Tia("i,c")
//                   // test
//                  ( (V_C1 * _4("<a k|gr|c l>") + V_C2 * _4("<k a|gr|c l>"))
//                    * _2("<k|I|l>")
//                  - (V_C1 * _4("<a k|r|p q>") + V_C2 * _4("<k a|r|p q>"))
//                    * _4("<p q|g|c k>")
//                  - (V_C1 * _4("<a k|r|a' n>") + V_C2 * _4("<k a|r|a' n>"))
//                    * _4("<a' n|g|c k>")
//                  - (V_C1 * _4("<a k|r|n a'>") + V_C2 * _4("<k a|r|n a'>"))
//                    * _4("<n a'|g|c k>")
//                  )

                  // 1/2 R^jk_A'B' g^A'B'_ik t1^a_j
                  - Tia("j,a") * V_ji("j,i")
//                  // test
//                  ( (V_C1 * _4("<j k|gr|i l>") + V_C2 * _4("<k j|gr|i l>"))
//                    * _2("<k|I|l>")
//                  - (V_C1 * _4("<j k|r|p q>") + V_C2 * _4("<k j|r|p q>"))
//                    * _4("<p q|g|i k>")
//                  - (V_C1 * _4("<j k|r|a' n>") + V_C2 * _4("<k j|r|a' n>"))
//                    * _4("<a' n|g|i k>")
//                  - (V_C1 * _4("<j k|r|n a'>") + V_C2 * _4("<k j|r|n a'>"))
//                    * _4("<n a'|g|i k>")
//                  )

                   // R^kl_ib' g^ab'_kc t1^c_l
                 - ( _4("<a b'|g|k c>")
                     * (V_C1 * _4("<k l|r|i b'>") + V_C2 * _4("<k l|r|b' i>"))
                   + _4("<a b'|g|c k>")
                     * (V_C1 * _4("<l k|r|i b'>") + V_C2 * _4("<l k|r|b' i>"))
                   // R^kl_ab' g^ib'_kc t1^c_l
                   + (V_C1 * _4("<a b'|r|k l>") + V_C2 * _4("<b' a|r|k l>"))
                     * _4("<k c|g|i b'>")
                   + (V_C1 * _4("<a b'|r|l k>") + V_C2 * _4("<b' a|r|l k>"))
                     * _4("<c k|g|i b'>")
                   ) * Tia("l,c")
                   ;
    //std::cout << "VT1_ai \n" << VT1_ai << std::endl;

    // can not use p q n a' for Vpq_rs
    TArray4 V_alcd = Vpq_rs(a,l,c,d);
    TArray4 V_klid = Vpq_rs(k,l,i,d);
    TArray4 T2 = _4("<i j|T2|a b>");
    TArray2 VT2_ai =
                     // 1/4 R^al_A'B' g^A'B'_il T^cd_il
                     V_alcd("a,l,c,d") * T2("i,l,c,d")
                     // 1/4 R^kl_A'B' g^A'B'_id T^ad_kl
                   - T2("k,l,a,d") * V_klid("k,l,i,d")
                     // 1/4 R^kl_ab' g^ib'_cd T^cd_kl
                   - ( V_C1 * _4("<a b'|r|k l>") + V_C2 * _4("<b' a|r|k l>"))
                     * T2("k,l,c,d") * _4("<c d|g|i b'>")
                     // 1/4 R^kl_ib' g^ab'_cd T^cd_kl
                   - _4("<a b'|g|c d>") * T2("k,l,c,d")
                     * ( V_C1 * _4("<k l|r|i b'>") + V_C2 * _4("<k l|r|b' i>"))
                   ;

#endif

    // F12 contribution to orbital response
    TArray2 Xam_f12 = 2.0 * ( Xam_Vcontri("a,m")
                            - Xam_Xcontri("a,m")
                            + Xam_Bcontri("a,m")
                            + gdf12_am("a,m")  // F12 density terms
                            // CT2_ai("a,i") + VT1_ai("a,i")+ VT2_ai("a,i") // CC F12 coupling terms
                            );

    TArray2 Dbn_f12(Xam_f12.get_world(), Xam_f12.trange());
    auto resnorm_f12 = cg_solver2(Orbital_relaxation_Abnam,
                                  Xam_f12,
                                  Dbn_f12,
                                  preconditioner,
                                  1e-10);

    if (compute_dipole) {
      TArray2 mu_z_apbp = _2("<a'|mu_z|b'>");
      // F12 density contribution to dipole
      const double mu_z_f12 = - dot(mu_z_ij("i,j"), D_f12_ij("i,j"))
                              + dot(mu_z_ab("a,b"), D_f12_ab("a,b"))
                              + dot(mu_z_apbp("a',b'"), D_f12_apbp("a',b'"))
                              + dot(mu_z_apb("a',b"), D_f12_apb("a',b")) * 2.0;
      std::cout << std::endl
                << "** mu_z (F12) = " << - mu_z_f12 * 2.0
                << std::endl;
      const double mu_z_Xam_f12 = dot(mu_z_am("a,m"), Dbn_f12("a,m"));
      std::cout << std::endl
                << "mu_z (F12 orbital response) = " << scprintf("%12.10f", - mu_z_Xam_f12 * 2.0)
                << std::endl << std::endl;
    }
    world_.gop.fence();

    // F12 contribution to quadrupoles
    if (compute_quadrupole) {
      TArray2 q_xx_apbp = _2("<a'|q_xx|b'>");
      TArray2 q_yy_apbp = _2("<a'|q_yy|b'>");
      TArray2 q_zz_apbp = _2("<a'|q_zz|b'>");

      TArray2 Qxx_apbp = q_xx_apbp("a',b'") - (q_zz_apbp("a',b'") + q_yy_apbp("a',b'")) * 0.5;
      TArray2 Qyy_apbp = q_yy_apbp("a',b'") - (q_zz_apbp("a',b'") + q_xx_apbp("a',b'")) * 0.5;
      TArray2 Qzz_apbp = q_zz_apbp("a',b'") - (q_xx_apbp("a',b'") + q_yy_apbp("a',b'")) * 0.5;
      TArray2 Qxy_apbp = _2("<a'|q_xy|b'>") * 1.5;
      TArray2 Qxz_apbp = _2("<a'|q_xz|b'>") * 1.5;
      TArray2 Qyz_apbp = _2("<a'|q_yz|b'>") * 1.5;

      const double q_xx_f12 = - dot(Qxx_ij("i,j"), D_f12_ij("i,j"))
                              + dot(Qxx_ab("a,b"), D_f12_ab("a,b"))
                              + dot(Qxx_apbp("a',b'"), D_f12_apbp("a',b'"))
                              + dot(Qxx_apb("a',b"), D_f12_apb("a',b")) * 2.0;
      const double q_yy_f12 = - dot(Qyy_ij("i,j"), D_f12_ij("i,j"))
                              + dot(Qyy_ab("a,b"), D_f12_ab("a,b"))
                              + dot(Qyy_apbp("a',b'"), D_f12_apbp("a',b'"))
                              + dot(Qyy_apb("a',b"), D_f12_apb("a',b")) * 2.0;
      const double q_zz_f12 = - dot(Qzz_ij("i,j"), D_f12_ij("i,j"))
                              + dot(Qzz_ab("a,b"), D_f12_ab("a,b"))
                              + dot(Qzz_apbp("a',b'"), D_f12_apbp("a',b'"))
                              + dot(Qzz_apb("a',b"), D_f12_apb("a',b")) * 2.0;
      const double q_xy_f12 = - dot(Qxy_ij("i,j"), D_f12_ij("i,j"))
                              + dot(Qxy_ab("a,b"), D_f12_ab("a,b"))
                              + dot(Qxy_apbp("a',b'"), D_f12_apbp("a',b'"))
                              + dot(Qxy_apb("a',b"), D_f12_apb("a',b")) * 2.0;
      const double q_xz_f12 = - dot(Qxz_ij("i,j"), D_f12_ij("i,j"))
                              + dot(Qxz_ab("a,b"), D_f12_ab("a,b"))
                              + dot(Qxz_apbp("a',b'"), D_f12_apbp("a',b'"))
                              + dot(Qxz_apb("a',b"), D_f12_apb("a',b")) * 2.0;
      const double q_yz_f12 = - dot(Qyz_ij("i,j"), D_f12_ij("i,j"))
                              + dot(Qyz_ab("a,b"), D_f12_ab("a,b"))
                              + dot(Qyz_apbp("a',b'"), D_f12_apbp("a',b'"))
                              + dot(Qyz_apb("a',b"), D_f12_apb("a',b")) * 2.0;

      const double q_xx_f12or = dot(Qxx_am("a,m"), Dbn_f12("a,m"));
      const double q_yy_f12or = dot(Qyy_am("a,m"), Dbn_f12("a,m"));
      const double q_zz_f12or = dot(Qzz_am("a,m"), Dbn_f12("a,m"));
      const double q_xy_f12or = dot(Qxy_am("a,m"), Dbn_f12("a,m"));
      const double q_xz_f12or = dot(Qxz_am("a,m"), Dbn_f12("a,m"));
      const double q_yz_f12or = dot(Qyz_am("a,m"), Dbn_f12("a,m"));

      std::cout << "traceless quadrupole moment (F12)" << std::endl
                << "q_xx (F12) = " << scprintf("%12.10f", - q_xx_f12 * 2.0)
                << "  q_yy (F12) = " << scprintf("%12.10f", - q_yy_f12 * 2.0)
                << "  q_zz (F12) = " << scprintf("%12.10f", - q_zz_f12 * 2.0)
                << std::endl
                << "q_xy (F12) = " << scprintf("%12.10f", - q_xy_f12 * 2.0)
                << "  q_xz (F12) = " << scprintf("%12.10f", - q_xz_f12 * 2.0)
                << "  q_yz (F12) = " << scprintf("%12.10f", - q_yz_f12 * 2.0)
                << std::endl;
      std::cout << "traceless quadrupole moment (F12 orbital response)" << std::endl
                << "q_xx (F12 or) = " << scprintf("%12.10f", - q_xx_f12or * 2.0)
                << "  q_yy (F12 or) = " << scprintf("%12.10f", - q_yy_f12or * 2.0)
                << "  q_zz (F12 or) = " << scprintf("%12.10f", - q_zz_f12or * 2.0)
                << std::endl
                << "q_xy (F12 or) = " << scprintf("%12.10f", - q_xy_f12or * 2.0)
                << "  q_xz (F12 or) = " << scprintf("%12.10f", - q_xz_f12or * 2.0)
                << "  q_yz (F12 or) = " << scprintf("%12.10f", - q_yz_f12or * 2.0)
                << std::endl << std::endl;
    }
    world_.gop.fence();

    if (compute_EFG) {
      TArray2 v_xx_apbp = _2("<a'|ddphi_xx|b'>");
      TArray2 v_yy_apbp = _2("<a'|ddphi_yy|b'>");
      TArray2 v_zz_apbp = _2("<a'|ddphi_zz|b'>");
      TArray2 v_xy_apbp = _2("<a'|ddphi_xy|b'>");
      TArray2 v_xz_apbp = _2("<a'|ddphi_xz|b'>");
      TArray2 v_yz_apbp = _2("<a'|ddphi_yz|b'>");

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

      std::cout << "electric gradient (F12)" << std::endl
                << "v_xx (F12) = " << scprintf("%12.10f", v_xx_f12 * 2.0)
                << "  v_yy (F12) = " << scprintf("%12.10f", v_yy_f12 * 2.0)
                << "  v_zz (F12) = " << scprintf("%12.10f", v_zz_f12 * 2.0)
                << std::endl
                << "v_xy (F12) = " << scprintf("%12.10f", v_xy_f12 * 2.0)
                << "  v_xz (F12) = " << scprintf("%12.10f", v_xz_f12 * 2.0)
                << "  v_yz (F12) = " << scprintf("%12.10f", v_yz_f12 * 2.0)
                << std::endl;
      std::cout << "electric gradient (F12 orbital response)" << std::endl
                << "v_xx (F12 or) = " << scprintf("%12.10f", v_xx_f12or * 2.0)
                << "  v_yy (F12 or) = " << scprintf("%12.10f", v_yy_f12or * 2.0)
                << "  v_zz (F12 or) = " << scprintf("%12.10f", v_zz_f12or * 2.0)
                << std::endl
                << "v_xy (F12 or) = " << scprintf("%12.10f", v_xy_f12or * 2.0)
                << "  v_xz (F12 or) = " << scprintf("%12.10f", v_xz_f12or * 2.0)
                << "  v_yz (F12 or) = " << scprintf("%12.10f", v_yz_f12or * 2.0)
                << std::endl << std::endl;
    }
    world_.gop.fence();
#endif

    /// this is just an example of how to compute the density
    TArray2 r2_i_j = _4("<i j|r|p q>") * _4("<k_F(p) j|r|p q>");

    return r2_i_j;
  }


  template <typename T>
  typename SingleReference_R12Intermediates<T>::TArray4
  SingleReference_R12Intermediates<T>::V_spinfree(bool symmetrize_p1_p2) {

    TArray4 V_ij_mn = _4("<i1 i2|gr|m1 m2>") - _4("<i1 i2|r|p1 p2>") * _4("<m1 m2|g|p1 p2>")
                    - _4("<i1 i2|r|m3_gamma(m) a'>") * _4("<m1 m2|g|m3 a'>");

    if (symmetrize_p1_p2)
      return 0.5 * (V_ij_mn("i1,i2,m1,m2") + V_ij_mn("i2,i1,m2,m1"));
    else
      return V_ij_mn;
  }

  template <typename T>
  typename SingleReference_R12Intermediates<T>::TArray4
  SingleReference_R12Intermediates<T>::X_spinfree(bool symmetrize_p1_p2) {

    TArray4 X_ij_kl = _4("<i1 i2|r2|j1 j2>") - _4("<i1 i2|r|p1 p2>") * _4("<j1 j2|r|p1 p2>")
                    - _4("<i1 i2|r|m3_gamma(m) a'>") * _4("<j1 j2|r|m3 a'>");

    if (symmetrize_p1_p2)
      return 0.5 * (X_ij_kl("i1,i2,j1,j2") + X_ij_kl("i2,i1,j2,j1"));
    else
      return X_ij_kl;
  }

  template <typename T>
  typename SingleReference_R12Intermediates<T>::TArray4
  SingleReference_R12Intermediates<T>::B_spinfree(bool symmetrize_p1_p2) {

    TArray4 B_ij_kl =

        // everything seems scaled up by factor of 2 relative to Eq.(12) in J. Chem. Phys. 135, 214105 (2011),
        // due to including particle 1 and particle 2 contributions?

        // diag                      Q
        _4("<i1 i2|rTr|j1 j2>") + 2.0 * _4("<i1 i2|r2|j1 j2_hJ(p')>")

        //           rKr
        -2.0 * _4("<i1 i2|r|r' s'>") * _4("<j1 j2|r|r' s'_K(p')>")

        //           rFr
        -2.0 * _4("<i1 i2|r|r s>") * _4("<j1 j2|r|r s_F(p)>")

        //           rFr_2, extra 2 due to bra-ket symmetrization
        -4.0 * _4("<i1 i2|r|r s>") * _4("<j1 j2|r|r s_F(a')>")

        //           rFGr
        - _4("<i1 i2|r|n_gamma(m) b'>") * _4("<j1 j2|r|n b'_F(a')>")

        //           rFGr_2
        - _4("<i1 i2|r|n_gamma(m) a'>") * _4("<j1 j2|r|n_F(p') a'>");

    B_ij_kl = 0.5 * (B_ij_kl("i1,i2,j1,j2") + B_ij_kl("i2,i1,j2,j1"));
    B_ij_kl = 0.5 * (B_ij_kl("i1,i2,j1,j2") + B_ij_kl("j1,j2,i1,i2"));

    if (symmetrize_p1_p2)
      return 0.5 * (B_ij_kl("i1,i2,j1,j2") + B_ij_kl("i2,i1,j2,j1"));
    else
      return B_ij_kl;
  }


}; // end of namespace sc

#endif // end of header guard


// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
