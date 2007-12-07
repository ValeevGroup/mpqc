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
#include <chemistry/qc/mbptr12/print.h>

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
  
    const bool vbs_eq_obs = r12info()->basis()->equiv(r12info()->basis_vir());
    const bool abs_eq_obs = r12info()->basis()->equiv(r12info()->basis_ri());
    const unsigned int maxnabs = r12info()->maxnabs();

    const unsigned int nf12 = corrfactor()->nfunctions();

    // some combinations are not implemented yet or are not sane
    if (!vbs_eq_obs && ansatz()->projector() == LinearR12::Projector_3)
	throw FeatureNotImplemented("B(C) cannot be evaluated yet when using ansatz 3 and VBS!=OBS",__FILE__,__LINE__);
  
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
	Ref<MOIndexSpace> xspace1 = xspace(spin1);
	Ref<MOIndexSpace> xspace2 = xspace(spin2);
	Ref<MOIndexSpace> vir1 = vir(spin1);
	Ref<MOIndexSpace> vir2 = vir(spin2);
	bool empty_vir_space = vir1->rank()==0 || vir2->rank()==0;

#if INCLUDE_Q
	// if can only use 1 RI index, h+J can be resolved by the OBS
	Ref<MOIndexSpace> hj_x1, hj_x2;
	if (maxnabs > 1) {
	    hj_x1 = hj_x_P(spin1);
	    hj_x2 = hj_x_P(spin2);
	}
	else {
	    hj_x1 = hj_x_p(spin1);
	    hj_x2 = hj_x_p(spin2);
	}    
	std::string Qlabel = prepend_spincase(spincase2,"Q(C) intermediate");
	tim_enter(Qlabel.c_str());
	ExEnv::out0() << endl << indent
		      << "Entered " << Qlabel << " evaluator" << endl;
	ExEnv::out0() << incindent;
    
	// compute Q = F12^2 (note F2_only = true in compute_X_ calls)
	RefSCMatrix Q;
	compute_X_(Q,spincase2,xspace1,xspace2,
		   xspace1,hj_x2,true);
	if (xspace1 != xspace2) {
	    compute_X_(Q,spincase2,xspace1,xspace2,
		       hj_x1,xspace2,true);
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

#if 0
	if (debug_ >= DefaultPrintThresholds::mostO4) {
#else
	    {
#endif
		std::string label = prepend_spincase(spincase2,"Q(C) contribution");
		Q.print(label.c_str());
	    }
	    B_[s].accumulate(Q); Q = 0;
#endif // INCLUDE_Q

#if INCLUDE_P
	    // compute P
	    // WARNING implemented only using CABS/CABS+ approach when OBS!=ABS
      
	    const LinearR12::ABSMethod absmethod = r12info()->abs_method();
	    if (!abs_eq_obs && 
		absmethod != LinearR12::ABS_CABS &&
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
		ribs1 = r12info()->ribs_space();
		ribs2 = r12info()->ribs_space();
	    }
	    RefSCMatrix P;
      
#if INCLUDE_P_PKP
	    // R_klPQ K_QR R_PRij is included with projectors 2 and 3
	    {
		Ref<MOIndexSpace> kribs1, kribs2;
		if (abs_eq_obs) {
		    kribs1 = K_p_p(spin1);
		    kribs2 = K_p_p(spin2);
		}
		else {
		    kribs1 = K_P_P(spin1);
		    kribs2 = K_P_P(spin2);
		}
		compute_FxF_(P,spincase2,
			     xspace1,xspace2,
			     xspace1,xspace2,
			     ribs1,ribs2,
			     ribs1,ribs2,
			     kribs1,kribs2);

                if (debug_ >= DefaultPrintThresholds::allO4)
                  P.print("P(incl R_klPQ K_QR R_PRij)");
	    }
#endif // INCLUDE_P_PKP

	    if (ansatz()->projector() == LinearR12::Projector_2) {
#if INCLUDE_P_PFP
		{
		    Ref<MOIndexSpace> fribs1,fribs2;
		    if (abs_eq_obs) {
			fribs1 = F_p_p(spin1);
			fribs2 = F_p_p(spin2);
		    }
		    else {
			fribs1 = F_P_P(spin1);
			fribs2 = F_P_P(spin2);
		    }
		    // R_klPm F_PQ R_Qmij
		    compute_FxF_(P,spincase2,
				 xspace1,xspace2,
				 xspace1,xspace2,
				 occ1,occ2,
				 ribs1,ribs2,
				 fribs1,fribs2);
                    if (debug_ >= DefaultPrintThresholds::allO4)
                      P.print("P(incl R_klPm F_PQ R_Qmij)");
		}
#endif // INCLUDE_P_PFP
	    }

	    {
#if INCLUDE_P_pFp
		// the form of the pFp term depends on the projector
		Ref<MOIndexSpace> z1 = (ansatz()->projector() == LinearR12::Projector_2 ? vir1 : orbs1);
		Ref<MOIndexSpace> z2 = (ansatz()->projector() == LinearR12::Projector_2 ? vir2 : orbs2);

		if (!empty_vir_space || ansatz()->projector() == LinearR12::Projector_3) {
		    // and on whether VBS==OBS
		    if (vbs_eq_obs) {
			Ref<MOIndexSpace> forbs1, forbs2;
			forbs1 = F_p_p(spin1);
			forbs2 = F_p_p(spin2);
			// R_klpa F_pq R_qaij
			compute_FxF_(P,spincase2,
				     xspace1,xspace2,
				     xspace1,xspace2,
				     z1,z2,
				     orbs1,orbs2,
				     forbs1,forbs2);
                        if (debug_ >= DefaultPrintThresholds::allO4)
                          P.print("P(incl R_klpa F_pq R_qaij)");
		    }
		    else {
			// else it's a sum of
			// R_klma F_mn R_naij
			{
			    Ref<MOIndexSpace> x1 = occ1;
			    Ref<MOIndexSpace> x2 = occ2;
			    Ref<MOIndexSpace> fx1 = F_m_m(spin1);
			    Ref<MOIndexSpace> fx2 = F_m_m(spin2);
			    compute_FxF_(P,spincase2,
					 xspace1,xspace2,
					 xspace1,xspace2,
					 z1,z2,
					 x1,x2,
					 fx1,fx2);
                            if (debug_ >= DefaultPrintThresholds::allO4)
                              P.print("P(incl R_klma F_mn R_naij)");
			}
			// R_klca F_cd R_daij
			{
			    Ref<MOIndexSpace> x1 = vir1;
			    Ref<MOIndexSpace> x2 = vir2;
			    Ref<MOIndexSpace> fx1 = F_a_a(spin1);
			    Ref<MOIndexSpace> fx2 = F_a_a(spin2);
			    compute_FxF_(P,spincase2,
					 xspace1,xspace2,
					 xspace1,xspace2,
					 z1,z2,
					 x1,x2,
					 fx1,fx2);
                            if (debug_ >= DefaultPrintThresholds::allO4)
                              P.print("P(incl R_klca F_cd R_daij)");
			}
			// R_klca F_cm R_maij + R_klma F_mc R_caij
			{
			    Ref<MOIndexSpace> x1 = occ1;
			    Ref<MOIndexSpace> x2 = occ2;
			    Ref<MOIndexSpace> fx1 = F_m_a(spin1);
			    Ref<MOIndexSpace> fx2 = F_m_a(spin2);
			    RefSCMatrix Ptmp;
			    compute_FxF_(Ptmp,spincase2,
					 xspace1,xspace2,
					 xspace1,xspace2,
					 z1,z2,
					 x1,x2,
					 fx1,fx2);
                            // bra-ket symmetrization will take care of
                            // the R_klma F_mc R_caij term
			    Ptmp.scale(2.0);
			    P.accumulate(Ptmp);
                            if (debug_ >= DefaultPrintThresholds::allO4)
                              P.print("P(incl 2 R_klca F_cm R_maij");
			}
		    } // VBS != OBS
		} // pFp contributes
#endif // INCLUDE_P_pFp
	    }

      
	    if (!abs_eq_obs) {
        
		Ref<MOIndexSpace> cabs1 = r12info()->ribs_space(spin1);
		Ref<MOIndexSpace> cabs2 = r12info()->ribs_space(spin2);
        
		if (ansatz()->projector() == LinearR12::Projector_2) {
#if INCLUDE_P_mFP
                    // R_klmA F_mP R_PAij + R_klPA F_Pm R_mAij
		    {
			Ref<MOIndexSpace> focc1, focc2;
			focc1 = F_m_P(spin1);
			focc2 = F_m_P(spin2);
			RefSCMatrix Ptmp;
			compute_FxF_(Ptmp,spincase2,
				     xspace1,xspace2,
				     xspace1,xspace2,
				     cabs1,cabs2,
				     occ1,occ2,
				     focc1,focc2);
                        // bra-ket symmetrization will take care of
                        // the R_klma F_mc R_caij term
                        Ptmp.scale(2.0);
			P.accumulate(Ptmp);
                        if (debug_ >= DefaultPrintThresholds::allO4)
                          P.print("P(incl 2 R_klmA F_mP R_PAij)");
		    }
#endif // INCLUDE_P_mFP
		}

		{
#if INCLUDE_P_pFA
		    // the form of the pFp term depends on the projector
		    Ref<MOIndexSpace> z1 = (ansatz()->projector() == LinearR12::Projector_2 ? vir1 : orbs1);
		    Ref<MOIndexSpace> z2 = (ansatz()->projector() == LinearR12::Projector_2 ? vir2 : orbs2);
		    
		    // R_klpz F_pA R_Azij
		    if (!empty_vir_space || ansatz()->projector() == LinearR12::Projector_3) {
			// and on whether VBS==OBS
			if (vbs_eq_obs) {
			    Ref<MOIndexSpace> forbs1 = F_p_A(spin1);
			    Ref<MOIndexSpace> forbs2 = F_p_A(spin2);
			    RefSCMatrix Ptmp;
			    compute_FxF_(Ptmp,spincase2,
					 xspace1,xspace2,
					 xspace1,xspace2,
					 z1,z2,
					 orbs1,orbs2,
					 forbs1,forbs2);
			    Ptmp.scale(2.0);
			    P.accumulate(Ptmp);
			} // VBS == OBS
			else { // VBS != OBS

			    RefSCMatrix Ptmp;
			    {
				Ref<MOIndexSpace> x1 = occ1;
				Ref<MOIndexSpace> x2 = occ2;
				Ref<MOIndexSpace> fx1 = F_m_A(spin1);
				Ref<MOIndexSpace> fx2 = F_m_A(spin2);
				compute_FxF_(Ptmp,spincase2,
					     xspace1,xspace2,
					     xspace1,xspace2,
					     z1,z2,
					     x1,x2,
					     fx1,fx2);
                                if (debug_ >= DefaultPrintThresholds::allO4)
                                  Ptmp.print("P(R_klmz F_mA R_Azij)");
			    }
			    {
				Ref<MOIndexSpace> x1 = vir1;
				Ref<MOIndexSpace> x2 = vir2;
				Ref<MOIndexSpace> fx1 = F_a_A(spin1);
				Ref<MOIndexSpace> fx2 = F_a_A(spin2);
				compute_FxF_(Ptmp,spincase2,
					     xspace1,xspace2,
					     xspace1,xspace2,
					     z1,z2,
					     x1,x2,
					     fx1,fx2);
                                if (debug_ >= DefaultPrintThresholds::allO4)
                                  Ptmp.print("P(R_klaz F_aA R_Azij)");
			    }
			    Ptmp.scale(2.0);
			    P.accumulate(Ptmp);
				
			} // VBS != OBS
                        if (debug_ >= DefaultPrintThresholds::allO4)
                          P.print("P(incl R_klpz F_pA R_Azij)");
		    } // pFA contributes
#endif // INCLUDE_P_pFA
		}

		P.scale(-1.0);

		if (ansatz()->projector() == LinearR12::Projector_2) {
#if INCLUDE_P_mFm
		    {
			Ref<MOIndexSpace> focc1 = F_m_m(spin1);
			Ref<MOIndexSpace> focc2 = F_m_m(spin2);
			// R_klmA F_mn R_nAij
			compute_FxF_(P,spincase2,
				     xspace1,xspace2,
				     xspace1,xspace2,
				     cabs1,cabs2,
				     occ1,occ2,
				     focc1,focc2);
                        if (debug_ >= DefaultPrintThresholds::allO4)
                          P.print("P(incl R_klmA F_mn R_nAij)");
		    }
#endif // INCLUDE_P_mFm
		}
        
	    } // ABS != OBS
	    else {
		P.scale(-1.0);
	    }
      
	    ExEnv::out0() << decindent;
	    ExEnv::out0() << indent << "Exited " << Plabel << " evaluator" << endl;
	    tim_exit(Plabel.c_str());

	    B_[s].accumulate(P); P = 0;
#endif // INCLUDE_P

	    // Bra-Ket symmetrize the B(C) contribution
	    B_[s].scale(0.5);
	    RefSCMatrix B_t = B_[s].t();
	    B_[s].accumulate(B_t);
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
