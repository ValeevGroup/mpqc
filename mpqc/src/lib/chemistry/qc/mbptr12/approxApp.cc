//
// approxApp.cc
//
// Copyright (C) 2006 Edward Valeev
//
// Author: Edward Valeev <edward.valeev@chemistry.gatech.edu>
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
#include <util/misc/timer.h>
#include <util/class/class.h>
#include <util/state/state.h>
#include <util/state/state_text.h>
#include <util/state/state_bin.h>
#include <math/scmat/local.h>
#include <math/scmat/matrix.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/mbptr12/blas.h>
#include <chemistry/qc/mbptr12/r12ia.h>
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

using namespace std;
using namespace sc;

#define INCLUDE_Q 1

void
R12IntEval::compute_BApp_()
{
  if (evaluated_)
    return;
  
  const bool vbs_eq_obs = r12info()->basis()->equiv(r12info()->basis_vir());
  const bool abs_eq_obs = r12info()->basis()->equiv(r12info()->basis_ri());
  const unsigned int maxnabs = r12info()->maxnabs();
  
  tim_enter("B(app. A'') intermediate");
  ExEnv::out0() << endl << indent
  << "Entered B(app. A'') intermediate evaluator" << endl;
  ExEnv::out0() << incindent;
  
  for(int s=0; s<nspincases2(); s++) {
    const SpinCase2 spincase2 = static_cast<SpinCase2>(s);
    const SpinCase1 spin1 = case1(spincase2);
    const SpinCase1 spin2 = case2(spincase2);
    Ref<MOIndexSpace> xspace1 = xspace(spin1);
    Ref<MOIndexSpace> xspace2 = xspace(spin2);

#if INCLUDE_Q

    std::string Qlabel = prepend_spincase(spincase2,"Q(A'') intermediate");
    tim_enter(Qlabel.c_str());
    ExEnv::out0() << endl << indent
                  << "Entered " << Qlabel << " evaluator" << endl;
    ExEnv::out0() << incindent;
    
    // compute Q = X_{xy}^{xy_{hj}}
    RefSCMatrix Q;
    if (maxnabs > 1) { // if can only use 2 RI index, h+J can be resolved by the RIBS
	Ref<MOIndexSpace> hj_x2 = hj_x_P(spin2);
	compute_X_(Q,spincase2,xspace1,xspace2,
		   xspace1,hj_x2);
    }
    else { // else do RI in orbital basis...
	if (vbs_eq_obs) { // which is just p if VBS == OBS
	    Ref<MOIndexSpace> hj_x2 = hj_x_p(spin2);
	    compute_X_(Q,spincase2,xspace1,xspace2,
		       xspace1,hj_x2);
	}
	else { // if VBS != OBS, p = m + a
	    Ref<MOIndexSpace> hj_x2 = hj_x_m(spin2);
	    compute_X_(Q,spincase2,xspace1,xspace2,
		       xspace1,hj_x2);
	    hj_x2 = hj_x_a(spin2);
	    compute_X_(Q,spincase2,xspace1,xspace2,
		       xspace1,hj_x2);
	}
    }

    if (xspace1 != xspace2) {
	if (maxnabs > 1) { // if can only use 2 RI index, h+J can be resolved by the RIBS
	    Ref<MOIndexSpace> hj_x1 = hj_x_P(spin1);
	    compute_X_(Q,spincase2,xspace1,xspace2,
		       hj_x1,xspace2);
	}
	else { // else do RI in orbital basis...
	    if (vbs_eq_obs) { // which is just p if VBS == OBS
		Ref<MOIndexSpace> hj_x1 = hj_x_p(spin1);
		compute_X_(Q,spincase2,xspace1,xspace2,
			   hj_x1,xspace2);
	    }
	    else { // if VBS != OBS, p = m + a
		Ref<MOIndexSpace> hj_x1 = hj_x_m(spin1);
		compute_X_(Q,spincase2,xspace1,xspace2,
			   hj_x1,xspace2);
		hj_x1 = hj_x_a(spin1);
		compute_X_(Q,spincase2,xspace1,xspace2,
			   hj_x1,xspace1);
	    }
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
    tim_exit(Qlabel.c_str());

    if (debug_ >= DefaultPrintThresholds::mostO4) {
      std::string label = prepend_spincase(spincase2,"Q(A'') contribution");
      Q.print(label.c_str());
    }
    B_[s].accumulate(Q); Q = 0;
#endif // INCLUDE_Q

    // Bra-Ket symmetrize the B(A'') contribution
    B_[s].scale(0.5);
    RefSCMatrix B_t = B_[s].t();
    B_[s].accumulate(B_t);
  }
  
  globally_sum_intermeds_();

  ExEnv::out0() << decindent;
  ExEnv::out0() << indent << "Exited B(app. A'') intermediate evaluator" << endl;

  tim_exit("B(app. A'') intermediate");
}

////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
