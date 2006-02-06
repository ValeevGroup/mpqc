//
// compute_a_gebc_vbs.cc
//
// Copyright (C) 2004 Edward Valeev
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
#include <math.h>
#include <limits.h>

#include <scconfig.h>
#include <util/misc/formio.h>
#include <util/misc/timer.h>
#include <util/class/class.h>
#include <util/state/state.h>
#include <util/state/state_text.h>
#include <util/state/state_bin.h>
#include <math/scmat/matrix.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/mbpt/bzerofast.h>
#include <chemistry/qc/mbptr12/r12ia.h>
#include <chemistry/qc/mbptr12/vxb_eval_info.h>
#include <chemistry/qc/mbptr12/pairiter.h>
#include <chemistry/qc/mbptr12/r12int_eval.h>

using namespace std;
using namespace sc;

void
R12IntEval::contrib_to_VXB_gebc_vbsneqobs_()
{
  if (evaluated_)
    return;

  // Compute VXB using new code
  using LinearR12::TwoParticleContraction;
  using LinearR12::Direct_Contraction;
  const LinearR12::ABSMethod absmethod = r12info()->abs_method();

  if (absmethod == LinearR12::ABS_ABS ||
      absmethod == LinearR12::ABS_ABSPlus)
    throw ProgrammingError("VXB != OBS only allowed with CABS/CABS+ RI method",__FILE__,__LINE__);
    
  for(int s=0; s<nspincases2(); s++) {
    const SpinCase2 spincase2 = static_cast<SpinCase2>(s);
    const SpinCase1 spin1 = case1(spincase2);
    const SpinCase1 spin2 = case2(spincase2);

    const Ref<MOIndexSpace>& occ1_act = r12info()->refinfo()->occ_act(spin1);
    const Ref<MOIndexSpace>& occ2_act = r12info()->refinfo()->occ_act(spin2);
    const Ref<MOIndexSpace>& occ1 = r12info()->refinfo()->occ(spin1);
    const Ref<MOIndexSpace>& occ2 = r12info()->refinfo()->occ(spin2);
    const Ref<MOIndexSpace>& vir1_act = r12info()->vir_act(spin1);
    const Ref<MOIndexSpace>& vir2_act = r12info()->vir_act(spin2);
    const Ref<MOIndexSpace>& ribs1 = r12info()->ribs_space(spin1);
    const Ref<MOIndexSpace>& ribs2 = r12info()->ribs_space(spin2);

    // (im|jn) contribution
    {
      Ref<TwoParticleContraction> tpcontract = new Direct_Contraction(occ1->rank(),occ2->rank(),-1.0);
      contrib_to_VXB_a_new_(occ1_act,occ1,occ2_act,occ2,
                            spincase2,tpcontract);
    }
    // (ia|jb) contribution
    {
      Ref<TwoParticleContraction> tpcontract = new Direct_Contraction(vir1_act->rank(),vir2_act->rank(),-1.0);
      contrib_to_VXB_a_new_(occ1_act,vir1_act,occ2_act,vir2_act,
                            spincase2,tpcontract);
    }
    // (im|ja) contribution
    {
      Ref<TwoParticleContraction> tpcontract = new Direct_Contraction(occ1->rank(),vir2_act->rank(),-1.0);
      contrib_to_VXB_a_new_(occ1_act,occ1,occ2_act,vir2_act,
                            spincase2,tpcontract);
      if (spincase2 == AlphaBeta && occ1_act != occ2_act) {
        Ref<TwoParticleContraction> tpcontract = new Direct_Contraction(vir1_act->rank(),occ2->rank(),-1.0);
        contrib_to_VXB_a_new_(occ1_act,vir1_act,occ2_act,occ2,
                              spincase2,tpcontract);
      }
    }
    // (im|jx) contribution
    {
      Ref<TwoParticleContraction> tpcontract = new Direct_Contraction(occ1->rank(),ribs2->rank(),-1.0);
      contrib_to_VXB_a_new_(occ1_act,occ1,occ2_act,ribs2,
                            spincase2,tpcontract);
      if (spincase2 == AlphaBeta && occ1_act != occ2_act) {
        Ref<TwoParticleContraction> tpcontract = new Direct_Contraction(ribs1->rank(),occ2->rank(),-1.0);
        contrib_to_VXB_a_new_(occ1_act,ribs1,occ2_act,occ2,
                              spincase2,tpcontract);
      }
    }
    
    if (debug_ > 1) {
      V_[s].print(prepend_spincase(static_cast<SpinCase2>(s),"V(diag+OBS+VBS+ABS) contribution").c_str());
      X_[s].print(prepend_spincase(static_cast<SpinCase2>(s),"X(diag+OBS+VBS+ABS) contribution").c_str());
      B_[s].print(prepend_spincase(static_cast<SpinCase2>(s),"B(diag+OBS+VBS+ABS) contribution").c_str());
    }
  }
  
  return;
}

////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
