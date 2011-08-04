//
// approxC_ansatz1.cc
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
#include <util/misc/timer.h>
#include <util/class/class.h>
#include <util/state/state.h>
#include <util/state/state_text.h>
#include <util/state/state_bin.h>
#include <math/scmat/local.h>
#include <math/scmat/matrix.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/mbptr12/r12wfnworld.h>
#include <math/mmisc/pairiter.h>
#include <chemistry/qc/mbptr12/r12int_eval.h>
#include <chemistry/qc/mbptr12/creator.h>
#include <chemistry/qc/mbptr12/container.h>
#include <chemistry/qc/mbptr12/compute_tbint_tensor.h>
#include <chemistry/qc/mbptr12/twoparticlecontraction.h>
#include <chemistry/qc/lcao/utils.h>
#include <chemistry/qc/lcao/utils.impl.h>
#include <util/misc/print.h>

using namespace std;
using namespace sc;

#define INCLUDE_Q 1
#define INCLUDE_P 1
#define INCLUDE_P_PKP 1
#define INCLUDE_P_PFP 1
#define INCLUDE_P_pFp 1
#define INCLUDE_P_mFP 1
#define INCLUDE_P_pFA 1
#define INCLUDE_P_mFm 1

void R12IntEval::compute_BC_ansatz1_() {
  if (evaluated_)
  return;

  const bool vbs_eq_obs = r12world()->basis()->equiv(r12world()->basis_vir());
  //const bool abs_eq_obs = r12world()->basis()->equiv(r12world()->basis_ri());
  const unsigned int maxnabs = r12world()->r12tech()->maxnabs();

  const unsigned int nf12 = corrfactor()->nfunctions();

  Timer timer("B(app. C) Ansatz1 intermediate");
  ExEnv::out0() << endl << indent
        << "Entered B(app. C) intermediate evaluator" << endl;
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
    std::string Qlabel = prepend_spincase(spincase2,"Q(C) Ansatz1 intermediate");
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

    ExEnv::out0() << decindent;
    ExEnv::out0() << indent << "Exited " << Qlabel << " evaluator" << endl;

    Qtimer.exit();


        if (debug_ >= DefaultPrintThresholds::mostO4) {
#if 0
        {
#endif
          std::string label = prepend_spincase(spincase2,"Q(C) contribution");
          ExEnv::out0() << indent << __FILE__ << ": "<<__LINE__<<"\n";
          Q.print(label.c_str());
        }
        B_[s].accumulate(Q); Q = 0;
#endif  /* INCLUDE_Q */
#if INCLUDE_P
        // compute P
        const std::string Plabel = prepend_spincase(spincase2,"P(C) Ansatz1 intermediate");
        Timer Ptimer(Plabel.c_str());
        ExEnv::out0() << endl << indent
              << "Entered " << Plabel << " evaluator" << endl;
        ExEnv::out0() << incindent;

        Ref<OrbitalSpace> ribs1 = r12world()->ribs_space();
        Ref<OrbitalSpace> ribs2 = r12world()->ribs_space();
        Ref<OrbitalSpace> cabs1 = r12world()->cabs_space(spin1);
        Ref<OrbitalSpace> cabs2 = r12world()->cabs_space(spin2);
        RefSCMatrix P;

#if INCLUDE_P_PKP
        {
          Ref<OrbitalSpace> kribs1 = K_P_P(spin1);
          Ref<OrbitalSpace> kribs2 = K_P_P(spin2);
          compute_FxF_(P,spincase2,
                       GG1space,GG2space,
                       GG1space,GG2space,
                       ribs1,ribs2,
                       ribs1,ribs2,
                       kribs1,kribs2);
          if (debug_ >= DefaultPrintThresholds::allO4)
            P.print("P(incl) R_pqp'q' k_p's' R_s'q'rs");
        }
#endif  /* INCLUDE_P_PKP */
#if INCLUDE_P_pFA
        {
          RefSCMatrix Ptmp;
          Ref<OrbitalSpace> f_p_cabs1 = F_p_A(spin1);
          Ref<OrbitalSpace> f_p_cabs2 = F_p_A(spin2);
          compute_FxF_(Ptmp,spincase2,
                       GG1space,GG2space,
                       GG1space,GG2space,
                       cabs1,cabs2,
                       orbs1,orbs2,
                       f_p_cabs1,f_p_cabs2);
          Ptmp.scale(2.0);
          P.accumulate(Ptmp);
          if (debug_ >= DefaultPrintThresholds::allO4)
            P.print("P(incl) 2 R_pqa'b' f_a'u R_ub'rs");
        }
#endif  /* INCLUDE_P_pFA */
#if INCLUDE_P_PFP
        {
          Ref<OrbitalSpace> fribs1 = F_P_P(spin1);
          Ref<OrbitalSpace> fribs2 = F_P_P(spin2);
          compute_FxF_(P,spincase2,
                       GG1space,GG2space,
                       GG1space,GG2space,
                       orbs1,orbs2,
                       ribs1,ribs2,
                       fribs1,fribs2);
          if (debug_ >= DefaultPrintThresholds::allO4)
            P.print("P(incl) R_pqp'u f_p'q' R_q'urs");
        }
#endif  /* INCLUDE_P_PFP */
#if INCLUDE_P_pFp
        {
          Ref<OrbitalSpace> forbs1 = F_p_p(spin1);
          Ref<OrbitalSpace> forbs2 = F_p_p(spin2);
          compute_FxF_(P,spincase2,
                       GG1space,GG2space,
                       GG1space,GG2space,
                       cabs1,cabs2,
                       orbs1,orbs2,
                       forbs1,forbs2);
          if (debug_ >= DefaultPrintThresholds::allO4)
            P.print("P(incl) R_pqub' f_uv R_vb'rs");
        }
#endif  /* INCLUDE_P_pFp */

        P.scale(-1.0);

        ExEnv::out0() << decindent;
        ExEnv::out0() << indent << "Exited " << Plabel << " evaluator" << endl;
        Ptimer.exit();


        if (debug_ >= DefaultPrintThresholds::mostO4) {
#if 0
        {
#endif
          std::string label = prepend_spincase(spincase2,"P(C) contribution");
          ExEnv::out0() << indent << __FILE__ << ": "<<__LINE__<<"\n";
          P.print(label.c_str());
        }
#endif  /* INCLUDE_P */

        // Accumulate P into B
        B_[s].accumulate(P); P=0;

        // Bra-Ket symmetrize the B(C) contribution
        B_[s].scale(0.5);
        RefSCMatrix B_t = B_[s].t();
        B_[s].accumulate(B_t);

  }

  timer.exit();
}
