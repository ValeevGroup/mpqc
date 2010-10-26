//
// approxC_GenRefansatz2.cc
//
// Copyright (C) 2008 Martin Torheyden
//
// Author: Martin Torheyden <mtorhey@vt.edu>
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
#include <util/class/class.h>
#include <util/state/state.h>
#include <util/state/state_text.h>
#include <util/state/state_bin.h>
#include <math/scmat/local.h>
#include <math/scmat/matrix.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/basis/integral.h>
#include <math/scmat/blas.h>
#include <chemistry/qc/mbptr12/r12wfnworld.h>
#include <chemistry/qc/mbptr12/pairiter.h>
#include <chemistry/qc/mbptr12/r12int_eval.h>
#include <chemistry/qc/mbptr12/creator.h>
#include <chemistry/qc/mbptr12/container.h>
#include <chemistry/qc/mbptr12/compute_tbint_tensor.h>
#include <chemistry/qc/mbptr12/contract_tbint_tensor.h>
#include <chemistry/qc/mbptr12/contract_tbint_tensors_to_obtensor.h>
#include <chemistry/qc/mbptr12/twoparticlecontraction.h>
#include <chemistry/qc/mbptr12/utils.h>
#include <chemistry/qc/mbptr12/utils.impl.h>
#include <chemistry/qc/mbptr12/print.h>

using namespace std;
using namespace sc;

#define INCLUDE_Q 1
#define INCLUDE_P 1
#define INCLUDE_P_PKP 1
#define INCLUDE_P_pFp 1
#define INCLUDE_P_pFA 1
#define INCLUDE_P_pgammaFgammap 1
#define INCLUDE_P_gammaF_p_A 1
#define INCLUDE_P_Fgamma_P_p 1

void R12IntEval::compute_BC_GenRefansatz2_() {
  if (evaluated_)
    return;

  Ref<R12IntEval> thisref(this);

  const bool vbs_eq_obs = r12world()->basis()->equiv(r12world()->basis_vir());
  //const bool abs_eq_obs = r12world()->basis()->equiv(r12world()->basis_ri());
  const unsigned int maxnabs = r12world()->r12tech()->maxnabs();

  const unsigned int nf12 = corrfactor()->nfunctions();
  Timer timer("B(app. C) general reference Ansatz2 intermediate");

  ExEnv::out0() << endl << indent
          << "Entered B(app. C) general reference Ansatz2 intermediate evaluator" << endl;
  ExEnv::out0() << incindent;

  for(int s=0; s<nspincases2(); s++) {
    const SpinCase2 spincase2 = static_cast<SpinCase2>(s);
    const SpinCase1 spin1 = case1(spincase2);
    const SpinCase1 spin2 = case2(spincase2);

    Ref<OrbitalSpace> occ1 = occ(spin1);
    Ref<OrbitalSpace> occ2 = occ(spin2);
    Ref<OrbitalSpace> orbs1 = orbs(spin1);
    Ref<OrbitalSpace> orbs2 = orbs(spin2);
    Ref<OrbitalSpace> GG1space = GGspace(spin1);
    Ref<OrbitalSpace> GG2space = GGspace(spin2);
    Ref<OrbitalSpace> xspace1 = xspace(spin1);
    Ref<OrbitalSpace> xspace2 = xspace(spin2);
    Ref<OrbitalSpace> vir1 = vir(spin1);
    Ref<OrbitalSpace> vir2 = vir(spin2);
    bool empty_vir_space = vir1->rank()==0 || vir2->rank()==0;

#if INCLUDE_Q
    // if can only use 1 RI index, h+J can be resolved by the OBS
    Ref<OrbitalSpace> hj_x1, hj_x2;
    if (maxnabs > 1) {
        hj_x1 = hj_x_P(spin1);
        hj_x2 = hj_x_P(spin2);
    }
    else {
        hj_x1 = hj_x_p(spin1);
        hj_x2 = hj_x_p(spin2);
    }
    std::string Qlabel = prepend_spincase(spincase2,"Q(C) Ansatz 2 intermediate");
    Timer Qtimer(Qlabel.c_str());
    ExEnv::out0() << endl << indent
              << "Entered " << Qlabel << " evaluator" << endl;
    ExEnv::out0() << incindent;

    // compute Q = F12^2 (note F2_only = true in compute_X_ calls)
    RefSCMatrix Q;
    compute_X_(Q,spincase2,GG1space,GG2space,
               GG1space,hj_x2,true);
    if (GG1space != GG2space) {
        compute_X_(Q,spincase2,GG1space,GG2space,
               hj_x1,GG2space,true);
    }
    else {
        Q.scale(2.0);
        if (spincase2 == AlphaBeta) {
          symmetrize<false>(Q,Q,GG1space,GG1space);
        }
    }

    //Q.scale(-1.0);

    ExEnv::out0() << decindent;
    ExEnv::out0() << indent << "Exited " << Qlabel << " evaluator" << endl;

    Qtimer.exit();

    if (debug_ >= DefaultPrintThresholds::mostO4) {
      std::string label = prepend_spincase(spincase2, "Q(C) contribution");
      Q.print(label.c_str());
    }
    B_[s].accumulate(Q);
    Q = 0;
#endif  /* INCLUDE_Q */
#if INCLUDE_P
    const std::string Plabel = prepend_spincase(spincase2,
                                                "P(C) Ansatz 2 intermediate");
    Timer Ptimer;
    ExEnv::out0() << endl << indent << "Entered " << Plabel << " evaluator"
        << endl;
    ExEnv::out0() << incindent;

    Ref<OrbitalSpace> ribs1 = r12world()->ribs_space();
    Ref<OrbitalSpace> ribs2 = r12world()->ribs_space();
    Ref<OrbitalSpace> cabs1 = r12world()->cabs_space(spin1);
    Ref<OrbitalSpace> cabs2 = r12world()->cabs_space(spin2);
    RefSCMatrix P = B_[s].clone();
    P.assign(0.0);

#if INCLUDE_P_PKP
    {
      RefSCMatrix Ptmp = P.clone();
      Ptmp.assign(0.0);
      Ref<OrbitalSpace> kribs1 = K_P_P(spin1);
      Ref<OrbitalSpace> kribs2 = K_P_P(spin2);
      compute_FxF_(Ptmp, spincase2, GG1space, GG2space, GG1space, GG2space,
                   ribs1, ribs2, ribs1, ribs2, kribs1, kribs2);
      if (debug_ >= DefaultPrintThresholds::allO4)
        Ptmp.print("R_pqp'q' k_p's' R_s'q'rs term");
      Ptmp.scale(-1.0);
      P.accumulate(Ptmp);
    }
#endif  /* INCLUDE_P_PKP */

#if INCLUDE_P_pFp
    {
      RefSCMatrix Ptmp;
      Ptmp = P.clone();
      Ptmp.assign(0.0);
      Ref<OrbitalSpace> f_p_p1 = F_p_p(spin1);
      Ref<OrbitalSpace> f_p_p2 = F_p_p(spin2);
      compute_FxF_(Ptmp, spincase2, GG1space, GG2space, GG1space, GG2space,
                   orbs1, orbs2, orbs1, orbs2, f_p_p1, f_p_p2);
      Ptmp.scale(-1.0);
      P.accumulate(Ptmp);
      if (debug_ >= DefaultPrintThresholds::allO4)
        Ptmp.print(
                   "contribution = -bar{r}^{q_3 p_3}_{v w} f^{q_2}_{q_3} bar{r}^{r s}_{q_2 p_2}.");
    }
#endif  /* INCLUDE_P_pFp */
#if INCLUDE_P_pFA
    {
      RefSCMatrix Ptmp;
      Ptmp = P.clone();
      Ptmp.assign(0.0);
      Ref<OrbitalSpace> f_p_A1 = F_p_A(spin1);
      Ref<OrbitalSpace> f_p_A2 = F_p_A(spin2);
      compute_FxF_(Ptmp, spincase2, GG1space, GG2space, GG1space, GG2space,
                   orbs1, orbs2, orbs1, orbs2, f_p_A1, f_p_A2);
      Ptmp.scale(-2.0);
      P.accumulate(Ptmp);
      if (debug_ >= DefaultPrintThresholds::allO4)
        Ptmp.print(
                   "contribution = -bar{r}^{alpha p_3}_{v w} f^{q_2}_{alpha} bar{r}^{r s}_{q_2 p_3} - bar{r}^{q_3 p_3}_{v w} f^{alpha}_{q_3} bar{r}^{r s}_{alpha p_3}");
    }
#endif  /* INCLUDE_P_pFA */
#if INCLUDE_P_pgammaFgammap
    {
      RefSCMatrix Ptmp;
      Ptmp = P.clone();
      Ptmp.assign(0.0);
      Ref<OrbitalSpace> gammafgamma_p_p1;
      Ref<OrbitalSpace> gammafgamma_p_p2;
      gammafgamma_p_p1 = gammaFgamma_p_p(spin1);
      gammafgamma_p_p2 = gammaFgamma_p_p(spin2);
      compute_FxF_(Ptmp, spincase2, GG1space, GG2space, GG1space, GG2space,
                   cabs1, cabs2, orbs1, orbs2, gammafgamma_p_p1,
                   gammafgamma_p_p2);
      P.accumulate(Ptmp);
      if (debug_ >= DefaultPrintThresholds::allO4)
        Ptmp.print(
                   "contribution = bar{r}^{p_3 alpha}_{v w} gamma^{q_3}_{p_3} f^{q_2}_{q_3} gamma^{p_2}_{q_2} bar{r}^{r s}_{p_2 alpha}");
    }
#endif  /* INCLUDE_P_pgammaFgammap */
#if INCLUDE_P_gammaF_p_A
    {
      // this contribution has to be computed "manually"
      RefSCMatrix Ptmp = P.clone();
      Ptmp.assign(0.0);
      RefSCMatrix Ptmp_p_A = P.clone();
      Ptmp_p_A.assign(0.0);
      Ref<OrbitalSpace> gamma_p_p1;
      Ref<OrbitalSpace> gamma_p_p2;
      gamma_p_p1 = gamma_p_p(spin1);
      gamma_p_p2 = gamma_p_p(spin2);
      Ref<OrbitalSpace> f_A_A1 = F_A_A(spin1);
      Ref<OrbitalSpace> f_A_A2 = F_A_A(spin2);

      std::vector<std::string> tforms_bra_p_A;
      {
        R12TwoBodyIntKeyCreator tform_creator(moints_runtime4(), GG1space,
                                              gamma_p_p1, GG2space, cabs2,
                                              corrfactor(), true);
        fill_container(tform_creator, tforms_bra_p_A);
      }
      std::vector<std::string> tforms_ket_p_A;
      {
        R12TwoBodyIntKeyCreator tform_creator(moints_runtime4(), GG1space,
                                              orbs1, GG2space, f_A_A2,
                                              corrfactor(), true);
        fill_container(tform_creator, tforms_ket_p_A);
      }
      contract_tbint_tensor<true, true> (Ptmp_p_A,
                                         corrfactor()->tbint_type_f12(),
                                         corrfactor()->tbint_type_f12(), -1.0,
                                         GG1space, GG2space, gamma_p_p1, cabs2,
                                         GG1space, GG2space, orbs1, f_A_A2,
                                         spincase2 != AlphaBeta,
                                         tforms_bra_p_A, tforms_ket_p_A);
      Ptmp.accumulate(Ptmp_p_A);
      P.accumulate(Ptmp_p_A);
      if (spincase2 == AlphaBeta) {
        RefSCMatrix Ptmp_A_p = P.clone();
        Ptmp_A_p.assign(0.0);
        std::vector<std::string> tforms_bra_A_p;
        {
          R12TwoBodyIntKeyCreator tform_creator(moints_runtime4(), GG1space,
                                                cabs1, GG2space, gamma_p_p2,
                                                corrfactor(), true);
          fill_container(tform_creator, tforms_bra_A_p);
        }
        std::vector<std::string> tforms_ket_A_p;
        {
          R12TwoBodyIntKeyCreator tform_creator(moints_runtime4(), GG1space,
                                                f_A_A1, GG2space, orbs2,
                                                corrfactor(), true);
          fill_container(tform_creator, tforms_ket_A_p);
        }
        contract_tbint_tensor<true, true> (Ptmp_A_p,
                                           corrfactor()->tbint_type_f12(),
                                           corrfactor()->tbint_type_f12(),
                                           -1.0, GG1space, GG2space, cabs1,
                                           gamma_p_p2, GG1space, GG2space,
                                           f_A_A1, orbs2, spincase2
                                               != AlphaBeta, tforms_bra_A_p,
                                           tforms_ket_A_p);
        Ptmp.accumulate(Ptmp_A_p);
        P.accumulate(Ptmp_A_p);
      }
      if (debug_ >= DefaultPrintThresholds::allO4)
        Ptmp.print(
                   "contribution = bar{r}^{p_3 alpha_3}_{v w} gamma^{p_2}_{p_3} f^{alpha_2}_{alpha_3} bar{r}^{r s}_{p_2 alpha_2}");
    }
#endif  /* INCLUDE_P_gammaF_p_A */
#if INCLUDE_P_Fgamma_P_p
    {
      RefSCMatrix Ptmp = P.clone();
      Ptmp.assign(0.0);
      Ref<OrbitalSpace> fgamma_p_P1 = Fgamma_p_P(spin1);
      Ref<OrbitalSpace> fgamma_p_P2 = Fgamma_p_P(spin2);
      compute_FxF_(Ptmp, spincase2, GG1space, GG2space, GG1space, GG2space,
                   cabs1, cabs2, orbs1, orbs2, fgamma_p_P1, fgamma_p_P2);
      Ptmp.scale(-2.0);
      P.accumulate(Ptmp);
      if (debug_ >= DefaultPrintThresholds::allO4)
        Ptmp.print(
                   "contribution = -bar{r}^{kappa alpha}_{v w} f^{q_2}_{kappa} gamma^{p_2}_{q_2} bar{r}^{r s}_{p_2 alpha} - bar{r}^{p_3 alpha}_{v w} gamma^{q_3}_{p_3} f^{kappa}_{q_3} bar{r}^{r s}_{kappa alpha}");
    }
#endif  /* INCLUDE_P_Fgamma_P_p */

    if (debug_ >= DefaultPrintThresholds::mostO4) {
      std::string label = prepend_spincase(spincase2, "P(C) contribution");
      P.print(label.c_str());
    }

    ExEnv::out0() << decindent;
    ExEnv::out0() << indent << "Exited " << Plabel << " evaluator" << endl;
    Ptimer.exit();

    // Accumulate P into B
    B_[s].accumulate(P);
    P = 0;

    Ptimer.exit();
#endif /* INCLUDE_P */

    // Bra-Ket symmetrize the B(C) contribution
    B_[s].scale(0.5);
    RefSCMatrix B_t = B_[s].t();
    B_[s].accumulate(B_t);
  }

  ExEnv::out0() << decindent;
  ExEnv::out0() << endl << indent
          << "Exited B(app. C) general reference Ansatz2 intermediate evaluator" << endl;
  timer.exit();
}
