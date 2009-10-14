//
// gbc_contribs.cc
//
// Copyright (C) 2004 Edward Valeev
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
#include <chemistry/qc/mbptr12/print_scmat_norms.h>
#include <chemistry/qc/mbptr12/utils.h>
#include <chemistry/qc/mbptr12/utils.impl.h>
#include <chemistry/qc/mbptr12/print.h>

using namespace std;
using namespace sc;

#define INCLUDE_GBC1 1
#define INCLUDE_GBC2 1

void
R12IntEval::compute_B_gbc_()
{
  if (abs_method() == LinearR12::ABS_ABS || abs_method() == LinearR12::ABS_ABSPlus)
    throw std::runtime_error("R12IntEval::compute_B_gbc_1_() -- B(GBC1) term can only be computed using a CABS (or CABS+) approach");

  if (evaluated_)
    return;

  Timer tim_B_GBC("B(GBC) intermediate");
  ExEnv::out0() << endl << indent
  << "Entered B(GBC) intermediate evaluator" << endl;
  ExEnv::out0() << incindent;

  for(int s=0; s<nspincases2(); s++) {
    const SpinCase2 spincase2 = static_cast<SpinCase2>(s);
    const SpinCase1 spin1 = case1(spincase2);
    const SpinCase1 spin2 = case2(spincase2);

    Ref<SingleRefInfo> refinfo = r12info()->refinfo();
    Ref<OrbitalSpace> occ1 = occ(spin1);
    Ref<OrbitalSpace> occ2 = occ(spin2);
    Ref<OrbitalSpace> orbs1 = orbs(spin1);
    Ref<OrbitalSpace> orbs2 = orbs(spin2);
    Ref<OrbitalSpace> cabs1 = r12info()->ribs_space(spin1);
    Ref<OrbitalSpace> cabs2 = r12info()->ribs_space(spin2);
    Ref<OrbitalSpace> xspace1 = xspace(spin1);
    Ref<OrbitalSpace> xspace2 = xspace(spin2);
    Ref<OrbitalSpace> vir1 = vir(spin1);
    Ref<OrbitalSpace> vir2 = vir(spin2);
    Ref<OrbitalSpace> F_m_A_1 = F_m_A(spin1);
    Ref<OrbitalSpace> F_m_A_2 = F_m_A(spin2);
    Ref<OrbitalSpace> F_x_A_1 = F_x_A(spin1);
    Ref<OrbitalSpace> F_x_A_2 = F_x_A(spin2);

    RefSCMatrix B_gbc1 = B_[s].clone(); B_gbc1.assign(0.0);
    RefSCMatrix B_gbc2 = B_[s].clone(); B_gbc2.assign(0.0);

#if INCLUDE_GBC1
    if (!omit_P()) {

    // R_klAb F_Am R_mbij
    compute_FxF_(B_gbc1,spincase2,
                 xspace1,xspace2,
                 xspace1,xspace2,
                 vir1,vir2,
                 occ1,occ2,
                 F_m_A_1,F_m_A_2);

    if (r12info()->maxnabs() >= 2) {
      // R_klAB F_Am R_mBij
      compute_FxF_(B_gbc1,spincase2,
                   xspace1,xspace2,
                   xspace1,xspace2,
                   cabs1,cabs2,
                   occ1,occ2,
                   F_m_A_1,F_m_A_2);
    }

    B_gbc1.scale(-1.0);

    } // omit P ?
#endif // include GBC1

#if INCLUDE_GBC2

    //
    // GBC2 contribution
    //

    compute_X_(B_gbc2,spincase2,xspace1,xspace2,
               xspace1,F_x_A_2);
    if (xspace1 != xspace2) {
      compute_X_(B_gbc2,spincase2,xspace1,xspace2,
                 F_x_A_1,xspace2);
    }
    else {
      B_gbc2.scale(2.0);
      if (spincase2 == AlphaBeta) {
        symmetrize<false>(B_gbc2,B_gbc2,xspace1,xspace1);
      }
    }

#endif // include GBC2 ?

    std::string label = prepend_spincase(spincase2,"B(GBC1) contribution");
    if (debug_ >= DefaultPrintThresholds::mostO4) {
      B_gbc1.print(label.c_str());
    }
    if (debug_ >= DefaultPrintThresholds::mostN0) {
      print_scmat_norms(B_gbc1,label.c_str());
    }
    label = prepend_spincase(spincase2,"B(GBC2) contribution");
    if (debug_ >= DefaultPrintThresholds::mostO4) {
      B_gbc2.print(label.c_str());
    }
    if (debug_ >= DefaultPrintThresholds::mostN0) {
      print_scmat_norms(B_gbc2,label.c_str());
    }

    RefSCMatrix B_gbc;
    {
      B_gbc = B_gbc1;
      B_gbc.accumulate(B_gbc2);
      B_gbc1 = 0; B_gbc2 = 0;
    }
    // Symmetrize the B contribution
    B_gbc.scale(0.5);
    RefSCMatrix B_gbc_t = B_gbc.t();
    B_[s].accumulate(B_gbc); B_[s].accumulate(B_gbc_t);

  }

  globally_sum_intermeds_();

  ExEnv::out0() << decindent;
  ExEnv::out0() << indent << "Exited B(GBC) intermediate evaluator" << endl;
}

////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
