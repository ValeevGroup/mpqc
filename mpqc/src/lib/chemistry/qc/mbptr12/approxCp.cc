//
// approxC_gbc.cc
//
// Copyright (C) 2006 Edward Valeev
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

#include <stdexcept>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#include <scconfig.h>
#include <util/misc/formio.h>
#include <util/misc/regtime.h>
#include <util/class/class.h>
#include <util/state/state.h>
#include <util/state/state_text.h>
#include <util/state/state_bin.h>
#include <math/scmat/local.h>
#include <math/scmat/matrix.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/mbptr12/blas.h>
#include <chemistry/qc/mbptr12/distarray4.h>
#include <chemistry/qc/mbptr12/vxb_eval_info.h>
#include <chemistry/qc/mbptr12/pairiter.h>
#include <chemistry/qc/mbptr12/r12int_eval.h>
#include <chemistry/qc/mbptr12/creator.h>
#include <chemistry/qc/mbptr12/container.h>
#include <chemistry/qc/mbptr12/compute_tbint_tensor.h>
#include <chemistry/qc/mbptr12/contract_tbint_tensor.h>
#include <chemistry/qc/mbptr12/twoparticlecontraction.h>
#include <chemistry/qc/mbptr12/utils.h>
#include <chemistry/qc/mbptr12/utils.impl.h>
#include <chemistry/qc/mbptr12/print.h>
#include <chemistry/qc/mbptr12/debug.h>

using namespace std;
using namespace sc;

#define DEBUG_PRINT_ALL_B_CONTRIBUTIONS 1
#define INCLUDE_Q 1
#define INCLUDE_P 1
# define INCLUDE_P_PHP 1
# define INCLUDE_Q_JK 1
# define INCLUDE_P_pFp 1
# define INCLUDE_P_mFP 1
# define INCLUDE_P_pFA 1
# define INCLUDE_P_mFm 1

void R12IntEval::compute_BCp_() {
  if (evaluated_)
    return;

  const bool vbs_eq_obs = r12info()->obs_eq_vbs();
  const bool abs_eq_obs = r12info()->obs_eq_ribs();
  const unsigned int maxnabs = r12info()->maxnabs();

  const unsigned int nf12 = corrfactor()->nfunctions();

  // some combinations are not implemented yet or are not sane
  if (ansatz()->projector() == LinearR12::Projector_3)
    throw FeatureNotImplemented(
                                "B(C') cannot be evaluated yet when using ansatz 3",
                                __FILE__, __LINE__);
  if (this->dk() > 0)
    throw FeatureNotImplemented(
                                "B(C') cannot be evaluated yet when df > 0",
                                __FILE__, __LINE__);

  Timer tim_B_app_C("B(app. C') intermediate");
  ExEnv::out0() << endl << indent
      << "Entered B(app. C', gbc) intermediate evaluator" << endl;
  ExEnv::out0() << incindent;

  for (int s = 0; s < nspincases2(); s++) {
    const SpinCase2 spincase2 = static_cast<SpinCase2> (s);
    const SpinCase1 spin1 = case1(spincase2);
    const SpinCase1 spin2 = case2(spincase2);

    Ref<SingleRefInfo> refinfo = r12info()->refinfo();
    Ref<OrbitalSpace> occ1 = refinfo->occ(spin1);
    Ref<OrbitalSpace> occ2 = refinfo->occ(spin2);
    Ref<OrbitalSpace> orbs1 = refinfo->orbs(spin1);
    Ref<OrbitalSpace> orbs2 = refinfo->orbs(spin2);
    Ref<OrbitalSpace> xspace1 = xspace(spin1);
    Ref<OrbitalSpace> xspace2 = xspace(spin2);
    Ref<OrbitalSpace> vir1 = vir(spin1);
    Ref<OrbitalSpace> vir2 = vir(spin2);
    Ref<OrbitalSpace> cabs1 = r12info()->ribs_space(spin1);
    Ref<OrbitalSpace> cabs2 = r12info()->ribs_space(spin2);

    bool empty_vir_space = vir1->rank() == 0 || vir2->rank() == 0;

#if INCLUDE_Q
    // if can only use 1 RI index, F can be resolved by the OBS
    Ref<OrbitalSpace> hJ_x1, hJ_x2;
    hJ_x1 = hj_x_P(spin1);
    hJ_x2 = hj_x_P(spin2);
    std::string Qlabel = prepend_spincase(spincase2, "Q(C') intermediate");
    Timer tim_Q(Qlabel);
    ExEnv::out0() << endl << indent << "Entered " << Qlabel << " evaluator"
        << endl;
    ExEnv::out0() << incindent;

    // compute Q = F12^2 (note F2_only = true in compute_X_ calls)
    RefSCMatrix Q;
    compute_X_(Q, spincase2, xspace1, xspace2, xspace1, hJ_x2, true);
    if (xspace1 != xspace2) {
      compute_X_(Q, spincase2, xspace1, xspace2, hJ_x1, xspace2, true);
    } else {
      Q.scale(2.0);
      if (spincase2 == AlphaBeta) {
        symmetrize<false> (Q, Q, xspace1, xspace1);
      }
    }

    ExEnv::out0() << decindent;
    ExEnv::out0() << indent << "Exited " << Qlabel << " evaluator" << endl;
    tim_Q.exit();

    if (debug_ >= DefaultPrintThresholds::mostO4
        || DEBUG_PRINT_ALL_B_CONTRIBUTIONS) {
      globally_sum_intermeds_();
      std::string label = prepend_spincase(spincase2, "Q(C') contribution");
      Q.print(label.c_str());
    }
    B_[s].accumulate(Q);
    Q = 0;
#endif // INCLUDE_Q
#if INCLUDE_P
    // compute P
    // WARNING implemented only using CABS/CABS+ approach when OBS!=ABS

    const LinearR12::ABSMethod absmethod = r12info()->abs_method();
    if (!abs_eq_obs && absmethod != LinearR12::ABS_CABS && absmethod
        != LinearR12::ABS_CABSPlus) {
      throw FeatureNotImplemented(
                                  "R12IntEval::compute_BCp_() -- approximation C must be used with absmethod=cabs/cabs+ if OBS!=ABS",
                                  __FILE__, __LINE__);
    }

    std::string Plabel = prepend_spincase(spincase2, "P(C') intermediate");
    Timer tim_P(Plabel);
    ExEnv::out0() << endl << indent << "Entered " << Plabel << " evaluator"
        << endl;
    ExEnv::out0() << incindent;

    Ref<OrbitalSpace> ribs1, ribs2;
    if (abs_eq_obs) {
      ribs1 = orbs1;
      ribs2 = orbs2;
    } else {
      ribs1 = r12info()->ribs_space();
      ribs2 = r12info()->ribs_space();
    }
    RefSCMatrix P = B_[s].clone();  P.assign(0.0);

#if INCLUDE_P_PHP
    {
      Ref<OrbitalSpace> hribs1, hribs2;
      hribs1 = h_P_P(spin1);
      hribs2 = h_P_P(spin2);
      // R_klPm h_PQ R_Qmij
      compute_FxF_(P, spincase2, xspace1, xspace2, xspace1, xspace2, occ1,
                   occ2, ribs1, ribs2, hribs1, hribs2);
      if (debug_ >= DefaultPrintThresholds::allO4
          || DEBUG_PRINT_ALL_B_CONTRIBUTIONS)
        P.print("P(incl R_klPm h_PQ R_Qmij)");
    }
#endif // INCLUDE_P_PHP

#if INCLUDE_Q_JK
    {
      Ref<OrbitalSpace> j_x_p1 = J_i_p(spin1);
      Ref<OrbitalSpace> j_x_p2 = J_i_p(spin2);
      // R_klPm R_Pmip (J)_pj
      {
        using namespace LinearR12;
        Ref<TwoParticleContraction> dircontract =
          new Direct_Contraction(occ1->rank(),ribs2->rank(),2.0);
        std::vector<std::string> tforms_bra;
        {
          R12TwoBodyIntKeyCreator tformkey_creator(r12info()->moints_runtime4(),xspace1,occ1,xspace2,ribs2,
                                                   r12info()->corrfactor(),true);
          fill_container(tformkey_creator,tforms_bra);
        }
        std::vector<std::string> tforms_ket;
        {
          R12TwoBodyIntKeyCreator tformkey_creator(r12info()->moints_runtime4(),xspace1,occ1,j_x_p2,ribs2,
                                                   r12info()->corrfactor(),true);
          fill_container(tformkey_creator,tforms_ket);
        }
        // contract
        contract_tbint_tensor<
          ManyBodyTensors::I_to_T,
          ManyBodyTensors::I_to_T,
          ManyBodyTensors::I_to_T,
          true,true,false>
          (
            P, corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_f12(),
            xspace1, xspace2,
            occ1, ribs2,
            xspace1, j_x_p2,
            occ1, ribs2,
            dircontract,
            spincase2!=AlphaBeta, tforms_bra, tforms_ket
          );
      }
      if (xspace1 != xspace2)
        abort();
      if (debug_ >= DefaultPrintThresholds::allO4
          || DEBUG_PRINT_ALL_B_CONTRIBUTIONS)
        P.print("P(incl R_klPm R_Pmip J_pj)");
    }
#endif // INCLUDE_Q_JK

#if INCLUDE_P_pFp
    {
      Ref<OrbitalSpace> forbs1, forbs2;
      forbs1 = hj_p_p(spin1);
      forbs2 = hj_p_p(spin2);
      // R_klpa hJ_pq R_qaij
      compute_FxF_(P, spincase2, xspace1, xspace2, xspace1, xspace2, vir1,
                   vir2, orbs1, orbs2, forbs1, forbs2);
      if (debug_ >= DefaultPrintThresholds::allO4
          || DEBUG_PRINT_ALL_B_CONTRIBUTIONS)
        P.print("P(incl R_klpa hJ_pq R_qaij)");
    }
#endif // INCLUDE_P_pFp
#if INCLUDE_P_mFP
    // R_klmA F_mp R_pAij + R_klpA F_pm R_mAij
    {
      Ref<OrbitalSpace> focc1, focc2;
      focc1 = hj_m_p(spin1);
      focc2 = hj_m_p(spin2);
      RefSCMatrix Ptmp;
      compute_FxF_(Ptmp, spincase2, xspace1, xspace2, xspace1, xspace2, cabs1,
                   cabs2, occ1, occ2, focc1, focc2);
      // bra-ket symmetrization will take care of the conjugate term
      Ptmp.scale(2.0);
      P.accumulate(Ptmp);
      if (debug_ >= DefaultPrintThresholds::allO4
          || DEBUG_PRINT_ALL_B_CONTRIBUTIONS)
        P.print("P(incl 2 R_klmA hJ_mp R_pAij)");
    }
#endif // INCLUDE_P_mFP
#if INCLUDE_P_pFA
    {
      Ref<OrbitalSpace> forbs1 = hj_p_A(spin1);
      Ref<OrbitalSpace> forbs2 = hj_p_A(spin2);
      RefSCMatrix Ptmp;
      compute_FxF_(Ptmp, spincase2, xspace1, xspace2, xspace1, xspace2, vir1,
                   vir2, orbs1, orbs2, forbs1, forbs2);
      Ptmp.scale(2.0);
      P.accumulate(Ptmp);
      if (debug_ >= DefaultPrintThresholds::allO4
          || DEBUG_PRINT_ALL_B_CONTRIBUTIONS)
        P.print("P(incl R_klpb hJ_pA R_Abij)");
    }
#endif // INCLUDE_P_pFA
    P.scale(-1.0);

#if INCLUDE_P_mFm
    {
      Ref<OrbitalSpace> focc1 = hj_m_m(spin1);
      Ref<OrbitalSpace> focc2 = hj_m_m(spin2);
      // R_klmA F_mn R_nAij
      compute_FxF_(P, spincase2, xspace1, xspace2, xspace1, xspace2, cabs1,
                   cabs2, occ1, occ2, focc1, focc2);
      if (debug_ >= DefaultPrintThresholds::allO4
          || DEBUG_PRINT_ALL_B_CONTRIBUTIONS)
        P.print("P(incl R_klmA F_mn R_nAij)");
    }
#endif // INCLUDE_P_mFm
    ExEnv::out0() << decindent;
    ExEnv::out0() << indent << "Exited " << Plabel << " evaluator" << endl;
    tim_P.exit();

    if (debug_ >= DefaultPrintThresholds::mostO4
        || DEBUG_PRINT_ALL_B_CONTRIBUTIONS) {
      std::string label = prepend_spincase(spincase2, "P(C') contribution");
      P.print(label.c_str());
    }

    B_[s].accumulate(P);
    P = 0;
#endif // INCLUDE_P
    // Bra-Ket symmetrize the B(C') contribution
    B_[s].scale(0.5);
    RefSCMatrix B_t = B_[s].t();
    B_[s].accumulate(B_t);
  } // end of spinpair loop

  globally_sum_intermeds_();

  ExEnv::out0() << decindent;
  ExEnv::out0() << indent << "Exited B(app. C') intermediate evaluator" << endl;

  tim_B_app_C.exit();
}

////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
