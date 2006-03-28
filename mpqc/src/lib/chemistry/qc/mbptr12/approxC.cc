//
// approxC.cc
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

void
R12IntEval::compute_BC_()
{
  if (evaluated_)
    return;
  
  const bool abs_eq_obs = r12info()->basis()->equiv(r12info()->basis_ri());
  
  tim_enter("B(app. C) intermediate");
  ExEnv::out0() << endl << indent
  << "Entered B(app. C) intermediate evaluator" << endl;
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

#if INCLUDE_Q
    Ref<MOIndexSpace> hjocc1_act_ribs = hjactocc_ribs(spin1);
    Ref<MOIndexSpace> hjocc2_act_ribs = hjactocc_ribs(spin2);
    
    std::string Qlabel = prepend_spincase(spincase2,"Q(C) intermediate");
    tim_enter(Qlabel.c_str());
    ExEnv::out0() << endl << indent
                  << "Entered " << Qlabel << " evaluator" << endl;
    ExEnv::out0() << incindent;
    
    // compute Q
    RefSCMatrix Q;
    compute_X_(Q,spincase2,occ1_act,occ2_act,
               occ1_act,hjocc2_act_ribs);
    if (occ1_act != occ2_act) {
      compute_X_(Q,spincase2,occ1_act,occ2_act,
                 hjocc1_act_ribs,occ2_act);
    }
    else {
      Q.scale(2.0);
      if (spincase2 == AlphaBeta) {
        symmetrize<false>(Q,Q,occ1_act,occ2_act);
      }
    }

    ExEnv::out0() << decindent;
    ExEnv::out0() << indent << "Exited " << Qlabel << " evaluator" << endl;
    tim_exit(Qlabel.c_str());

    if (debug_ > 1) {
      std::string label = prepend_spincase(spincase2,"Q(C) contribution");
      Q.print(label.c_str());
    }
    BC_[s].accumulate(Q); Q = 0;
#endif // INCLUDE_Q

#if INCLUDE_P
    // compute P
    // WARNING implemented only using CABS/CABS+ approach
    if (!omit_P()) {
      
      const LinearR12::ABSMethod absmethod = r12info()->abs_method();
      if (absmethod != LinearR12::ABS_CABS &&
          absmethod != LinearR12::ABS_CABSPlus) {
            throw FeatureNotImplemented("R12IntEval::compute_BC_() -- approximation C must be used with absmethod=cabs/cabs+ if OBS!=ABS",__FILE__,__LINE__);
      }
      
      std::string Plabel = prepend_spincase(spincase2,"P(C) intermediate");
      tim_enter(Plabel.c_str());
      ExEnv::out0() << endl << indent
                    << "Entered " << Plabel << " evaluator" << endl;
      ExEnv::out0() << incindent;
      
      Ref<MOIndexSpace> ribs1, ribs2;
      if (abs_eq_obs) {
        ribs1 = orbs1;
        ribs2 = orbs2;
      }
      else {
        ribs1 = r12info()->abs_space();
        ribs2 = r12info()->abs_space();
      }
      RefSCMatrix P;
      
#if INCLUDE_P_PKP
      {
      Ref<MOIndexSpace> kribs1 = kribs_ribs(spin1);
      Ref<MOIndexSpace> kribs2 = kribs_ribs(spin2);
      // R_klPQ K_QR R_PRij
      compute_FxF_(P,spincase2,
                   occ1_act,occ2_act,
                   occ1_act,occ2_act,
                   ribs1,ribs2,
                   ribs1,ribs2,
                   kribs1,kribs2);
      }
#endif // INCLUDE_P_PKP
#if INCLUDE_P_PFP
      {
      Ref<MOIndexSpace> fribs1 = fribs_ribs(spin1);
      Ref<MOIndexSpace> fribs2 = fribs_ribs(spin2);
      // R_klPm F_PQ R_Qmij
      compute_FxF_(P,spincase2,
                   occ1_act,occ2_act,
                   occ1_act,occ2_act,
                   occ1,occ2,
                   ribs1,ribs2,
                   fribs1,fribs2);
      }
#endif // INCLUDE_P_PFP
#if INCLUDE_P_pFp
      {
      Ref<MOIndexSpace> forbs1 = fobs_obs(spin1);
      Ref<MOIndexSpace> forbs2 = fobs_obs(spin2);
      // R_klpa F_pq R_qaij
      compute_FxF_(P,spincase2,
                   occ1_act,occ2_act,
                   occ1_act,occ2_act,
                   occ1,occ2,
                   orbs1,orbs2,
                   forbs1,forbs2);
      }
#endif // INCLUDE_P_pFp

      if (!abs_eq_obs) {
        
        Ref<MOIndexSpace> cabs1 = r12info()->ribs_space(spin1);
        Ref<MOIndexSpace> cabs2 = r12info()->ribs_space(spin2);
        
#if INCLUDE_P_mFP
        {
          Ref<MOIndexSpace> focc1 = focc_ribs(spin1);
          Ref<MOIndexSpace> focc2 = focc_ribs(spin2);
          // R_klmA F_mP R_PAij
          compute_FxF_(P,spincase2,
                       occ1_act,occ2_act,
                       occ1_act,occ2_act,
                       cabs1,cabs2,
                       occ1,occ2,
                       focc1,focc2);
        }
#endif // INCLUDE_P_mFP
#if INCLUDE_P_pFA
        {
          Ref<MOIndexSpace> forbs1 = fobs_cabs(spin1);
          Ref<MOIndexSpace> forbs2 = fobs_cabs(spin2);
          // R_klpa F_pA R_Aaij
          compute_FxF_(P,spincase2,
                       occ1_act,occ2_act,
                       occ1_act,occ2_act,
                       vir1,vir2,
                       orbs1,orbs2,
                       forbs1,forbs2);
        }
#endif // INCLUDE_P_pFA
        P.scale(-1.0);
#if INCLUDE_P_mFm
        {
          Ref<MOIndexSpace> focc1 = focc_occ(spin1);
          Ref<MOIndexSpace> focc2 = focc_occ(spin2);
          // R_klmA F_mn R_nAij
          compute_FxF_(P,spincase2,
                       occ1_act,occ2_act,
                       occ1_act,occ2_act,
                       cabs1,cabs2,
                       occ1,occ2,
                       focc1,focc2);
        }
#endif // INCLUDE_P_mFm
        
      }
      else {
        P.scale(-1.0);
      }
      
      ExEnv::out0() << decindent;
      ExEnv::out0() << indent << "Exited " << Plabel << " evaluator" << endl;
      tim_exit(Plabel.c_str());

      BC_[s].accumulate(P); P = 0;
    }
#endif // INCLUDE_P

    // Bra-Ket symmetrize the B(C) contribution
    BC_[s].scale(0.5);
    RefSCMatrix BC_t = BC_[s].t();
    BC_[s].accumulate(BC_t);
  }
  
  globally_sum_intermeds_();

  ExEnv::out0() << decindent;
  ExEnv::out0() << indent << "Exited B(app. C) intermediate evaluator" << endl;

  tim_exit("B(app. C) intermediate");
}

////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
