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
#define INCLUDE_P 1
#define INCLUDE_P_AKB 1
#define INCLUDE_P_AKb 1
#define INCLUDE_P_aKB 1
#define TEST_P_AKB 0
#define TEST_P_AKb 0
#define TEST_P_aKB 0
#define TEST_P_aKb 0
#define TEST_P_AKb_EXPL 0

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
    Ref<MOIndexSpace> vir1 = vir(spin1);
    Ref<MOIndexSpace> vir2 = vir(spin2);
    const Ref<MOIndexSpace>& xspace1 = xspace(spin1);
    const Ref<MOIndexSpace>& xspace2 = xspace(spin2);

#if INCLUDE_Q
    Ref<MOIndexSpace> Kxp1 = K_x_p(spin1);
    Ref<MOIndexSpace> Kxp2 = K_x_p(spin2);
    
    const LinearR12::ABSMethod absmethod = r12info()->abs_method();
    const bool include_Kp = (absmethod == LinearR12::ABS_CABS ||
                             absmethod == LinearR12::ABS_CABSPlus) || abs_eq_obs;
    
    std::string Qlabel = prepend_spincase(spincase2,"Q intermediate");
    tim_enter(Qlabel.c_str());
    ExEnv::out0() << endl << indent
                  << "Entered " << Qlabel << " evaluator" << endl;
    ExEnv::out0() << incindent;
    
    // compute Q
    RefSCMatrix Q;
    if (include_Kp) {
	compute_X_(Q,spincase2,xspace1,xspace2,
 		               xspace1,Kxp2);
    }
    if (!abs_eq_obs) {
	Ref<MOIndexSpace> KxA2 = K_x_A(spin2);
	compute_X_(Q,spincase2,xspace1,xspace2,
                               xspace1,KxA2);
    }
    if (xspace1 != xspace2) {
      if (include_Kp) {
	compute_X_(Q,spincase2,xspace1,xspace2,
		               Kxp1,xspace2);
      }
      if (!abs_eq_obs) {
	Ref<MOIndexSpace> KxA1 = K_x_A(spin1);
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
    tim_exit(Qlabel.c_str());

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
      
      const LinearR12::ABSMethod absmethod = r12info()->abs_method();
      if (absmethod != LinearR12::ABS_CABS &&
          absmethod != LinearR12::ABS_CABSPlus) {
            throw FeatureNotImplemented("R12IntEval::compute_BB_() -- approximation B must be used with absmethod=cabs/cabs+ if OBS!=ABS",__FILE__,__LINE__);
      }
      
      std::string Plabel = prepend_spincase(spincase2,"P intermediate");
      tim_enter(Plabel.c_str());
      ExEnv::out0() << endl << indent
                    << "Entered " << Plabel << " evaluator" << endl;
      ExEnv::out0() << incindent;
      
      Ref<MOIndexSpace> cabs1 = r12info()->ribs_space(spin1);
      Ref<MOIndexSpace> cabs2 = r12info()->ribs_space(spin2);
      
      RefSCMatrix P;
      if (r12info()->maxnabs() < 2) {

        Ref<MOIndexSpace> kvir1_obs = kvir_obs(spin1);
        Ref<MOIndexSpace> kvir2_obs = kvir_obs(spin2);

        // R_klpB K_pa R_aBij
        compute_FxF_(P,spincase2,
                     xspace1,xspace2,
                     xspace1,xspace2,
                     cabs1,cabs2,
                     vir1,vir2,
                     kvir1_obs,kvir2_obs);
      }
      else {

#if INCLUDE_P_AKB || INCLUDE_P_AKb
        Ref<MOIndexSpace> kcabs1 = kribs(spin1);
        Ref<MOIndexSpace> kcabs2 = kribs(spin2);
#endif

#if INCLUDE_P_AKB
        // R_klPB K_PA R_ABij
        compute_FxF_(P,spincase2,
                     xspace1,xspace2,
                     xspace1,xspace2,
                     cabs1,cabs2,
                     cabs1,cabs2,
                     kcabs1,kcabs2);

#if TEST_P_AKB
        P.print("R_klPB K_PA R_ABij");

        RefSCMatrix P_AKB = P.clone(); P_AKB.assign(0.0);
        Ref<MOIndexSpace> kabs1 = kribs_T(spin1);
        Ref<MOIndexSpace> kabs2 = kribs_T(spin2);
        Ref<MOIndexSpace> abs1 = r12info()->abs_space();
        Ref<MOIndexSpace> abs2 = r12info()->abs_space();
        // R_klAB K_AP R_PBij
        compute_FxF_(P_AKB,spincase2,
                     xspace1,xspace2,
                     xspace1,xspace2,
                     cabs1,cabs2,
                     abs1,abs2,
                     kabs1,kabs2);
        P_AKB.print("R_klAB K_AP R_PBij");
#endif
#endif // INCLUDE_P_AKB

#if INCLUDE_P_AKb
        // R_klPb K_PA R_Abij
        compute_FxF_(P,spincase2,
                     xspace1,xspace2,
                     xspace1,xspace2,
                     vir1,vir2,
                     cabs1,cabs2,
                     kcabs1,kcabs2);

#if TEST_P_AKb
        P.print("R_klPb K_PA R_Abij");

	RefSCMatrix P_AKb = P.clone(); P_AKb.assign(0.0);
	Ref<MOIndexSpace> kabs1 = kribs_T(spin1);
        Ref<MOIndexSpace> kabs2 = kribs_T(spin2);
        Ref<MOIndexSpace> abs1 = r12info()->abs_space();
        Ref<MOIndexSpace> abs2 = r12info()->abs_space();
        // R_klAb K_AP R_Pbij
        compute_FxF_(P_AKb,spincase2,
                     xspace1,xspace2,
                     xspace1,xspace2,
                     vir1,vir2,
                     abs1,abs2,
                     kabs1,kabs2);
        P_AKb.print("R_klAb K_AP R_Pbij");
#endif

#if TEST_P_AKb_EXPL
        Ref<R12IntEval> thisref(this);
        // (i a |i a') tforms
        std::vector<  Ref<TwoBodyMOIntsTransform> > tforms_iaiA;
        {
          NewTransformCreator tform_creator(thisref,xspace1,vir1,xspace2,cabs2,true);
          fill_container(tform_creator,tforms_iaiA);
        }
        // (i a |i p') tforms
        std::vector<  Ref<TwoBodyMOIntsTransform> > tforms_iaiP;
        {
          NewTransformCreator tform_creator(thisref,xspace1,vir1,xspace2,abs2,true);
          fill_container(tform_creator,tforms_iaiP);
        }

        Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
        RefSCDimension dim_aA = new SCDimension(vir1->rank() * cabs2->rank());
        RefSCDimension dim_aP = new SCDimension(vir1->rank() * abs2->rank());
        RefSCMatrix F12_ij_aA = localkit->matrix(dim_oo(spincase2),dim_aA);  F12_ij_aA.assign(0.0);
        RefSCMatrix F12_ij_aP = localkit->matrix(dim_oo(spincase2),dim_aP);  F12_ij_aP.assign(0.0);
        compute_tbint_tensor<ManyBodyTensors::I_to_T,true,false>(
          F12_ij_aA, corrfactor()->tbint_type_f12(),
          xspace1, vir1,
          xspace2, cabs2,
          spincase2!=AlphaBeta, tforms_iaiA
        );
        F12_ij_aA.print("F12 (aA)");
        compute_tbint_tensor<ManyBodyTensors::I_to_T,true,false>(
          F12_ij_aP, corrfactor()->tbint_type_f12(),
          xspace1, vir1,
          xspace2, abs2,
          spincase2!=AlphaBeta, tforms_iaiP
        );
        F12_ij_aP.print("F12 (aP)");

        // Test these too
        // (i a'|i a) tforms
        std::vector<  Ref<TwoBodyMOIntsTransform> > tforms_iAia;
        {
          NewTransformCreator tform_creator(thisref,xspace1,cabs1,xspace2,vir2,true);
          fill_container(tform_creator,tforms_iAia);
        }
        // (i p'|i a) tforms
        std::vector<  Ref<TwoBodyMOIntsTransform> > tforms_iPia;
        {
          NewTransformCreator tform_creator(thisref,xspace1,abs1,xspace2,vir2,true);
          fill_container(tform_creator,tforms_iPia);
        }
        RefSCMatrix F12_ij_Aa = localkit->matrix(dim_oo(spincase2),dim_aA);  F12_ij_Aa.assign(0.0);
        RefSCMatrix F12_ij_Pa = localkit->matrix(dim_oo(spincase2),dim_aP);  F12_ij_Pa.assign(0.0);
        compute_tbint_tensor<ManyBodyTensors::I_to_T,true,false>(
          F12_ij_Aa, corrfactor()->tbint_type_f12(),
          xspace1, cabs1,
          xspace2, vir2,
          spincase2!=AlphaBeta, tforms_iAia
        );
        F12_ij_Aa.print("F12 (Aa)");
        compute_tbint_tensor<ManyBodyTensors::I_to_T,true,false>(
          F12_ij_Pa, corrfactor()->tbint_type_f12(),
          xspace1, abs1,
          xspace2, vir2,
          spincase2!=AlphaBeta, tforms_iPia
        );
        F12_ij_Pa.print("F12 (Pa)");

        RefSCMatrix K_PA = exchange_(occ2,abs2,cabs2);
        K_PA.print("Exchange (ABS/RI-BS)");

        RefSCMatrix K_PA_local = localkit->matrix(K_PA.rowdim(),K_PA.coldim());
        double* tmparr = new double[K_PA.rowdim().n() * K_PA.coldim().n()];
        K_PA.convert(tmparr);
        K_PA_local.assign(tmparr);
        delete[] tmparr;

        const unsigned int nij = dim_oo(spincase2).n();
        RefSCMatrix F12_ij_aAK = F12_ij_aA.clone();  F12_ij_aAK.assign(0.0);
        double* aP_block = new double[dim_aP.n()];
        double* aA_block = new double[dim_aA.n()];
        RefSCMatrix tmp_aP = localkit->matrix(new SCDimension(vir1->rank()), K_PA.rowdim());
        for(unsigned int ij=0; ij<nij; ij++) {
          RefSCVector rij_aP = F12_ij_aP.get_row(ij);
          rij_aP.convert(aP_block);
          tmp_aP.assign(aP_block);
          RefSCMatrix tmp_aA = tmp_aP * K_PA_local;
          tmp_aA.convert(aA_block);
          RefSCVector rij_aA = F12_ij_aAK.get_row(ij);
          rij_aA.assign(aA_block);
          F12_ij_aAK.assign_row(rij_aA,ij);
        }
        F12_ij_aAK.print("F12 (aAK)");
        RefSCMatrix P_AKb_expl = F12_ij_aA * F12_ij_aAK.t();
        P_AKb_expl.scale(2.0);
        symmetrize<false>(P_AKb_expl,P_AKb_expl,xspace1,xspace2);
        P_AKb_expl.print("R_klPb K_PA R_Abij (expl)");
#endif
        
#endif // INCLUDE_P_AKb

#if INCLUDE_P_aKB
        Ref<MOIndexSpace> kvir1_ribs = kvir_ribs(spin1);
        Ref<MOIndexSpace> kvir2_ribs = kvir_ribs(spin2);
        // R_klPB K_Pa R_aBij
        compute_FxF_(P,spincase2,
                     xspace1,xspace2,
                     xspace1,xspace2,
                     cabs1,cabs2,
                     vir1,vir2,
                     kvir1_ribs,kvir2_ribs);
#if TEST_P_aKB
        P.print("R_klPB K_Pa R_aBij");

        RefSCMatrix P_aKB = P.clone(); P_aKB.assign(0.0);
        Ref<MOIndexSpace> kabs1 = kvir_ribs_T(spin1);
        Ref<MOIndexSpace> kabs2 = kvir_ribs_T(spin2);
        Ref<MOIndexSpace> abs1 = r12info()->abs_space();
        Ref<MOIndexSpace> abs2 = r12info()->abs_space();
        // R_klaB K_aP R_PBij
        compute_FxF_(P_aKB,spincase2,
                     xspace1,xspace2,
                     xspace1,xspace2,
                     cabs1,cabs2,
                     abs1,abs2,
                     kabs1,kabs2);
        P_aKB.print("R_klaB K_aP R_PBij");
#endif

#endif // INCLUDE_P_aKB
      }

#if TEST_P_aKb
      RefSCMatrix P_aKb = BB_[s].clone(); P_aKb.assign(0.0);
      Ref<MOIndexSpace> kvir1_ribs = kvir_ribs(spin1);
      Ref<MOIndexSpace> kvir2_ribs = kvir_ribs(spin2);
      // R_klPb K_Pa R_abij
      compute_FxF_(P_aKb,spincase2,
		   xspace1,xspace2,
		   xspace1,xspace2,
                   vir1,vir2,
                   vir1,vir2,
                   kvir1_ribs,kvir2_ribs);
      P_aKb.print("P_aKb test");
#endif

      P.scale(-1.0);
      
      ExEnv::out0() << decindent;
      ExEnv::out0() << indent << "Exited " << Plabel << " evaluator" << endl;
      tim_exit(Plabel.c_str());

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

  tim_exit("B(app. B) intermediate");
}

////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
