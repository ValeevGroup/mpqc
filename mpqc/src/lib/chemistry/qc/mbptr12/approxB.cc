//
// approxB.cc
//
// Copyright (C) 2005 Edward Valeev
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

using namespace std;
using namespace sc;

#define SYMMETRIZE 1

void
R12IntEval::compute_BB_()
{
  if (evaluated_)
    return;
  
  const bool abs_eq_obs = r12info()->basis()->equiv(r12info()->basis_ri());
  
  tim_enter("B(app. B) intermediate");
  ExEnv::out0() << endl << indent
  << "Entered B(app. B) intermediate evaluator" << endl;
  ExEnv::out0() << incindent;
  
  for(int s=0; s<nspincases2(); s++) {
    const SpinCase2 spincase2 = static_cast<SpinCase2>(s);
    const SpinCase1 spin1 = case1(spincase2);
    const SpinCase1 spin2 = case2(spincase2);
    
    Ref<SingleRefInfo> refinfo = r12info()->refinfo();
    Ref<MOIndexSpace> occ1 = refinfo->occ(spin1);
    Ref<MOIndexSpace> occ2 = refinfo->occ(spin2);
    Ref<MOIndexSpace> orbs1 = refinfo->orbs(spin1);
    Ref<MOIndexSpace> orbs2 = refinfo->orbs(spin2);
    Ref<MOIndexSpace> occ1_act = occ_act(spin1);
    Ref<MOIndexSpace> occ2_act = occ_act(spin2);
    Ref<MOIndexSpace> vir1 = vir(spin1);
    Ref<MOIndexSpace> vir2 = vir(spin2);
    Ref<MOIndexSpace> kocc1_act_obs = kocc_act_obs(spin1);
    Ref<MOIndexSpace> kocc2_act_obs = kocc_act_obs(spin2);
    
    const LinearR12::ABSMethod absmethod = r12info()->abs_method();
    const bool include_Kp = (absmethod == LinearR12::ABS_CABS ||
                             absmethod == LinearR12::ABS_CABSPlus) || abs_eq_obs;
    
    // compute Q
    RefSCMatrix Q;
    if (include_Kp) {
      compute_X_(Q,spincase2,occ1_act,occ2_act,
                             occ1_act,kocc2_act_obs);
    }
    if (!abs_eq_obs) {
      Ref<MOIndexSpace> kocc2_act = kocc_act(spin2);
      compute_X_(Q,spincase2,occ1_act,occ2_act,
                             occ1_act,kocc2_act);
    }
    if (occ1_act != occ2_act) {
      if (include_Kp) {
        compute_X_(Q,spincase2,occ1_act,occ2_act,
                               kocc1_act_obs,occ2_act);
      }
      if (!abs_eq_obs) {
        Ref<MOIndexSpace> kocc1_act = kocc_act(spin1);
        compute_X_(Q,spincase2,occ1_act,occ2_act,
                               kocc1_act,occ2_act);
      }
    }
    else {
      Q.scale(2.0);
      symmetrize<false>(Q,Q,occ1_act,occ2_act);
    }
    if (debug_ > 1) {
      std::string label = prepend_spincase(spincase2,"B(Q) contribution");
      Q.print(label.c_str());
    }
    BB_[s].accumulate(Q); Q = 0;
    
    // compute P
    // WARNING implemented only using CABS/CABS+ approach
    if (!abs_eq_obs) {
      
      const LinearR12::ABSMethod absmethod = r12info()->abs_method();
      if (absmethod != LinearR12::ABS_CABS ||
          absmethod != LinearR12::ABS_CABSPlus) {
            throw FeatureNotImplemented("R12IntEval::compute_BB_() -- approximation B must be used with absmethod=cabs/cabs+ if OBS!=ABS",__FILE__,__LINE__);
      }
      
      Ref<MOIndexSpace> cabs1 = r12info()->ribs_space(spin1);
      Ref<MOIndexSpace> cabs2 = r12info()->ribs_space(spin2);
      
      RefSCMatrix P;
      if (r12info()->maxnabs() < 2) {

        Ref<MOIndexSpace> kvir1_obs = kvir_obs(spin1);
        Ref<MOIndexSpace> kvir2_obs = kvir_obs(spin2);

        // R_klpB K_pa R_aBij
        compute_FxF_(P,spincase2,
                     occ1_act,occ2_act,
                     occ1_act,occ2_act,
                     cabs1,cabs2,
                     vir1,vir2,
                     kvir1_obs,kvir2_obs);
      }
      else {
        
        Ref<MOIndexSpace> kcabs1 = kribs(spin1);
        Ref<MOIndexSpace> kcabs2 = kribs(spin2);
        Ref<MOIndexSpace> kvir1_ribs = kvir_ribs(spin1);
        Ref<MOIndexSpace> kvir2_ribs = kvir_ribs(spin2);
        
        // R_klPB K_PA R_ABij
        compute_FxF_(P,spincase2,
                     occ1_act,occ2_act,
                     occ1_act,occ2_act,
                     cabs1,cabs2,
                     cabs1,cabs2,
                     kcabs1,kcabs2);
        // R_klPb K_PA R_Abij
        compute_FxF_(P,spincase2,
                     occ1_act,occ2_act,
                     occ1_act,occ2_act,
                     vir1,vir2,
                     cabs1,cabs2,
                     kcabs1,kcabs2);
        // R_klPB K_Pa R_aBij
        compute_FxF_(P,spincase2,
                     occ1_act,occ2_act,
                     occ1_act,occ2_act,
                     cabs1,cabs2,
                     vir1,vir2,
                     kvir1_ribs,kvir2_ribs);
      }
      
      P.scale(-1.0);
      BB_[s].accumulate(P); P = 0;
    }
    
    // Bra-Ket symmetrize the B(B) contribution
    BB_[s].scale(0.5);
    RefSCMatrix BB_t = BB_[s].t();
    BB_[s].accumulate(BB_t);
  }
  
  globally_sum_intermeds_();

  ExEnv::out0() << decindent;
  ExEnv::out0() << indent << "Exited B(app. B) intermediate evaluator" << endl;

  tim_exit("B(app. B) intermediate");
}

////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
