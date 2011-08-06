//
// approxB.cc
//
// Copyright (C) 2005 Edward Valeev
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
#include <math/distarray4/distarray4.h>
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
#define INCLUDE_P_AKB 1
#define INCLUDE_P_AKb 1
#define INCLUDE_P_aKB 1

void
R12IntEval::compute_BB_()
{
  if (evaluated_)
    return;

  const bool vbs_eq_obs = r12world()->obs_eq_vbs();
  const bool abs_eq_obs = r12world()->obs_eq_ribs();

  Timer tim_B_app_B("B(app. B) intermediate");
  ExEnv::out0() << endl << indent
  << "Entered B(app. B) intermediate evaluator" << endl;
  ExEnv::out0() << incindent;

  for(int s=0; s<nspincases2(); s++) {
    const SpinCase2 spincase2 = static_cast<SpinCase2>(s);
    const SpinCase1 spin1 = case1(spincase2);
    const SpinCase1 spin2 = case2(spincase2);

    Ref<OrbitalSpace> occ1 = occ(spin1);
    Ref<OrbitalSpace> occ2 = occ(spin2);
    Ref<OrbitalSpace> orbs1 = orbs(spin1);
    Ref<OrbitalSpace> orbs2 = orbs(spin2);
    Ref<OrbitalSpace> vir1 = vir(spin1);
    Ref<OrbitalSpace> vir2 = vir(spin2);
    const Ref<OrbitalSpace>& xspace1 = xspace(spin1);
    const Ref<OrbitalSpace>& xspace2 = xspace(spin2);

#if INCLUDE_Q
    const bool include_Kp = abs_eq_obs;

    std::string Qlabel = prepend_spincase(spincase2,"Q intermediate");
    Timer tim_Q(Qlabel);
    ExEnv::out0() << endl << indent
                  << "Entered " << Qlabel << " evaluator" << endl;
    ExEnv::out0() << incindent;

    // compute Q
    RefSCMatrix Q;
    if (include_Kp) {
	if (vbs_eq_obs) {
	    Ref<OrbitalSpace> Kxp2 = K_x_p(spin2);
	    compute_X_(Q,spincase2,xspace1,xspace2,
		       xspace1,Kxp2);
	} // VBS == OBS
	else { // VBX != OBS
	    Ref<OrbitalSpace> Kxm2 = K_x_m(spin2);
	    compute_X_(Q,spincase2,xspace1,xspace2,
		       xspace1,Kxm2);
	    Ref<OrbitalSpace> Kxa2 = K_x_a(spin2);
	    compute_X_(Q,spincase2,xspace1,xspace2,
		       xspace1,Kxa2);
	}
    }
    if (!abs_eq_obs) {
	Ref<OrbitalSpace> KxA2 = K_x_A(spin2);
	compute_X_(Q,spincase2,xspace1,xspace2,
                               xspace1,KxA2);
    }
    if (xspace1 != xspace2) {
	if (include_Kp) {
	    if (vbs_eq_obs) {
		Ref<OrbitalSpace> Kxp1 = K_x_p(spin1);
		compute_X_(Q,spincase2,xspace1,xspace2,
			   Kxp1,xspace2);
	    } // VBS == OBS
	    else { // VBS != OBS
		Ref<OrbitalSpace> Kxm1 = K_x_m(spin1);
		compute_X_(Q,spincase2,xspace1,xspace2,
			   Kxm1,xspace2);
		Ref<OrbitalSpace> Kxa1 = K_x_a(spin1);
		compute_X_(Q,spincase2,xspace1,xspace2,
			   Kxa1,xspace2);
	    }
	}
	if (!abs_eq_obs) {
	    Ref<OrbitalSpace> KxA1 = K_x_A(spin1);
	    compute_X_(Q,spincase2,xspace1,xspace2,
		       KxA1,xspace2);
	}
    }
    else {
	Q.scale(2.0);
	if (spincase2 == AlphaBeta) {
	    symmetrize<false>(Q,Q,xspace1,xspace1);
	}
    }

    ExEnv::out0() << decindent;
    ExEnv::out0() << indent << "Exited " << Qlabel << " evaluator" << endl;
    tim_Q.exit();

    if (debug_ >= DefaultPrintThresholds::mostO4) {
      std::string label = prepend_spincase(spincase2,"B(Q) contribution");
      Q.print(label.c_str());
    }
    BB_[s].accumulate(Q); Q = 0;
#endif // INCLUDE_Q

#if INCLUDE_P
    // compute P
    // WARNING implemented only using CABS/CABS+ approach
    if (!abs_eq_obs && !omit_P()) {

      std::string Plabel = prepend_spincase(spincase2,"P intermediate");
      Timer tim_P(Plabel);
      ExEnv::out0() << endl << indent
                    << "Entered " << Plabel << " evaluator" << endl;
      ExEnv::out0() << incindent;

      Ref<OrbitalSpace> cabs1 = r12world()->cabs_space(spin1);
      Ref<OrbitalSpace> cabs2 = r12world()->cabs_space(spin2);

      RefSCMatrix P;
      if (r12world()->r12tech()->maxnabs() < 2) {

	  if (vbs_eq_obs) {
	      Ref<OrbitalSpace> kvir1_obs = K_a_p(spin1);
	      Ref<OrbitalSpace> kvir2_obs = K_a_p(spin2);

	      // R_klpB K_pa R_aBij
	      compute_FxF_(P,spincase2,
			   xspace1,xspace2,
			   xspace1,xspace2,
			   cabs1,cabs2,
			   vir1,vir2,
			   kvir1_obs,kvir2_obs);
	  } // VBS == OBS
	  else { // VBS != OBS
	      Ref<OrbitalSpace> Kma1 = K_m_a(spin1);
	      Ref<OrbitalSpace> Kma2 = K_m_a(spin2);

	      // R_klmB K_ma R_aBij
	      compute_FxF_(P,spincase2,
			   xspace1,xspace2,
			   xspace1,xspace2,
			   cabs1,cabs2,
			   occ1,occ2,
			   Kma1,Kma2);

	      Ref<OrbitalSpace> Kaa1 = K_a_a(spin1);
	      Ref<OrbitalSpace> Kaa2 = K_a_a(spin2);

	      // R_klbB K_ba R_aBij
	      compute_FxF_(P,spincase2,
			   xspace1,xspace2,
			   xspace1,xspace2,
			   cabs1,cabs2,
			   vir1,vir2,
			   Kaa1,Kaa2);
	  }

      }
      else {

#if INCLUDE_P_AKB || INCLUDE_P_AKb
        Ref<OrbitalSpace> kcabs1 = K_A_P(spin1);
        Ref<OrbitalSpace> kcabs2 = K_A_P(spin2);
#endif

#if INCLUDE_P_AKB
        // R_klPB K_PA R_ABij
        compute_FxF_(P,spincase2,
                     xspace1,xspace2,
                     xspace1,xspace2,
                     cabs1,cabs2,
                     cabs1,cabs2,
                     kcabs1,kcabs2);

#endif // INCLUDE_P_AKB

#if INCLUDE_P_AKb
        // R_klPb K_PA R_Abij
        compute_FxF_(P,spincase2,
                     xspace1,xspace2,
                     xspace1,xspace2,
                     vir1,vir2,
                     cabs1,cabs2,
                     kcabs1,kcabs2);

#endif // INCLUDE_P_AKb

#if INCLUDE_P_aKB
        Ref<OrbitalSpace> kvir1_ribs = K_a_P(spin1);
        Ref<OrbitalSpace> kvir2_ribs = K_a_P(spin2);
        // R_klPB K_Pa R_aBij
        compute_FxF_(P,spincase2,
                     xspace1,xspace2,
                     xspace1,xspace2,
                     cabs1,cabs2,
                     vir1,vir2,
                     kvir1_ribs,kvir2_ribs);

#endif // INCLUDE_P_aKB
      }

      P.scale(-1.0);

      ExEnv::out0() << decindent;
      ExEnv::out0() << indent << "Exited " << Plabel << " evaluator" << endl;
      tim_P.exit();

      BB_[s].accumulate(P); P = 0;
    }
#endif // INCLUDE_P

    // Bra-Ket symmetrize the B(B) contribution
    BB_[s].scale(0.5);
    RefSCMatrix BB_t = BB_[s].t();
    BB_[s].accumulate(BB_t);
  }

  globally_sum_intermeds_();

  ExEnv::out0() << decindent;
  ExEnv::out0() << indent << "Exited B(app. B) intermediate evaluator" << endl;

  tim_B_app_B.exit();
}

////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
