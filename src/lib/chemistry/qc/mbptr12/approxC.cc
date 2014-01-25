//
// approxC.cc
//
// Copyright (C) 2006 Edward Valeev
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

#include <mpqc_config.h>
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
#include <chemistry/qc/mbptr12/debug.h>

using namespace std;
using namespace sc;

#define DEBUG_PRINT_ALL_B_CONTRIBUTIONS 0
#define INCLUDE_Q 1
#define INCLUDE_P 1
#define INCLUDE_P_PKP 1
#if !EVALUATE_MV_VIA_RI
# define INCLUDE_P_PFP 1
# define INCLUDE_P_pFp 1
# define INCLUDE_P_mFP 1
# define INCLUDE_P_pFA 1
# define INCLUDE_P_mFm 1
#else
# define INCLUDE_P_PFP 0
# define INCLUDE_P_pFp 0
# define INCLUDE_P_mFP 0
# define INCLUDE_P_pFA 0
# define INCLUDE_P_mFm 0
#endif

void
R12IntEval::compute_BC_()
{
    if (evaluated_)
	return;

    const bool vbs_eq_obs = r12world()->obs_eq_vbs();
    const bool abs_eq_obs = r12world()->obs_eq_ribs();
    const unsigned int maxnabs = r12world()->r12tech()->maxnabs();

    const unsigned int nf12 = corrfactor()->nfunctions();

    // some combinations are not implemented yet or are not sane
    if (!vbs_eq_obs && ansatz()->projector() == R12Technology::Projector_3)
	throw FeatureNotImplemented("B(C) cannot be evaluated yet when using ansatz 3 and VBS!=OBS",__FILE__,__LINE__);

    Timer tim_B_app_C("B(app. C) intermediate");
    ExEnv::out0() << endl << indent
		  << "Entered B(app. C) intermediate evaluator" << endl;
    ExEnv::out0() << incindent;

    //
    // if I'm relativistic, all relativistic terms will be treated as K, hence
    // 1) make new space Kr = K - (H(r) - H(nr)) = K - dH
    // 2) make new space hJnr = h(nr) + J = (h(r) + J) - (H(r) - H(nr)) = hJ - dH
    //
    Ref<OrbitalSpace> Kr[NSpinCases1];
    Ref<OrbitalSpace> hJnr[NSpinCases1];
    if (this->dk() > 0) {

      const R12Technology::H0_dk_approx_pauli H0_dk_approx_pauli = r12world()->r12tech()->H0_dk_approx();

      const int nspins1 = this->nspincases1();
      for (int s = 0; s < nspins1; ++s) {
        const SpinCase1 spin = static_cast<SpinCase1> (s);
        // which space is used as RIBS?
        Ref<OrbitalSpace> ribs = (!abs_eq_obs) ? r12world()->ribs_space()
                                               : this->orbs(spin);
        // get AO space for RIBS
        const Ref<OrbitalSpace>& aoribs =
          this->ao_registry()->value(ribs->basis());

        Ref<OrbitalSpace> kribs = (!abs_eq_obs) ? K_P_P(spin) : K_p_p(spin);
        if (H0_dk_approx_pauli == R12Technology::H0_dk_approx_pauli_false) { // correct exchange
          // compute dH = H(rel) - H(nonrel) in AO basis
          RefSCMatrix Hr = this->fock(aoribs, aoribs, spin, 0.0, 0.0, 1.0);
          const std::string nonrel_hkey =
            ParsedOneBodyIntKey::key(aoribs->id(), aoribs->id(),
                                     std::string("H"));
          RefSCMatrix Hnr = fockbuild_runtime()->get(nonrel_hkey);
          RefSCMatrix dH = Hnr.clone();
          dH->convert(Hr);
          Hnr.scale(-1.0);
          dH.accumulate(Hnr);
          Hr = 0;
          Hnr = 0;
          // transform to MO basis
          RefSCMatrix dH_mo = ribs->coefs().t() * dH * ribs->coefs();

          // Make new space K - dH
          std::string id = ribs->id();
          id += "_Kr(";
          id += ribs->id();
          id += ")";
          ExEnv::out0() << "id = " << id << endl;
          std::string name = "(Kr)-weighted space";
#if 1
          Kr[spin] = new OrbitalSpace(id, name, kribs, kribs->coefs() - ribs->coefs() * dH_mo, ribs->basis());
#else
          Kr[spin] = new OrbitalSpace(id, name, kribs, kribs->coefs(),
                                    ribs->basis());
#endif

          this->r12world()->world()->tfactory()->orbital_registry()->add(make_keyspace_pair(Kr[s]));
        }
        else { // use pure exchange
          Kr[spin] = kribs;
        }

        Ref<OrbitalSpace> hJ_x_P = (!abs_eq_obs) ? hj_x_P(spin) : hj_x_p(spin);
        if (H0_dk_approx_pauli == R12Technology::H0_dk_approx_pauli_false) { // use nonrelativistic hamiltonian in h+J
          const Ref<OrbitalSpace>& x = GGspace(spin);
          const Ref<OrbitalSpace>& aox = this->ao_registry()->value(x->basis());
          // compute dH = H(rel) - H(nonrel) in AO basis
          RefSCMatrix Hr = this->fock(aox, aoribs, spin, 0.0, 0.0, 1.0);
          const std::string nonrel_hkey =
            ParsedOneBodyIntKey::key(aox->id(), aoribs->id(),
                                     std::string("H"));
          RefSCMatrix Hnr = fockbuild_runtime()->get(nonrel_hkey);
          RefSCMatrix dH = Hnr.clone();
          dH->convert(Hr);
          Hnr.scale(-1.0);
          dH.accumulate(Hnr);
          Hr = 0;
          Hnr = 0;
          // transform to MO basis, transpose so that x dimension is coldim
          RefSCMatrix dH_mo = x->coefs().t() * dH * ribs->coefs();
          dH_mo = dH_mo.t();

          // Make new space hJ - dH
          std::string id = x->id();
          id += "_hJnr(";
          id += ribs->id();
          id += ")";
          ExEnv::out0() << "id = " << id << endl;
          std::string name = "(hJnr)-weighted space";
#if 1
          hJnr[spin] = new OrbitalSpace(id, name, hJ_x_P, hJ_x_P->coefs() - ribs->coefs() * dH_mo, hJ_x_P->basis());
#else
          hJnr[spin] = new OrbitalSpace(id, name, hJ_x_P, hJ_x_P->coefs(), hJ_x_P->basis());
#endif

          this->r12world()->world()->tfactory()->orbital_registry()->add(make_keyspace_pair(hJnr[s]));
        }
        else if (H0_dk_approx_pauli == R12Technology::H0_dk_approx_pauli_fHf) { // use pauli hamitonian in h+J
          const Ref<OrbitalSpace>& x = GGspace(spin);
          const Ref<OrbitalSpace>& aox = this->ao_registry()->value(x->basis());
          // compute dH = H(rel) - H(pauli) in AO basis
          RefSCMatrix Hr = this->fock(aox, aoribs, spin, 0.0, 0.0, 1.0);
          const int override_pauli = 1;
          RefSCMatrix Hpauli = this->fock(aox, aoribs, spin, 0.0, 0.0, 1.0, override_pauli);
          Hpauli.scale(-1.0);  Hpauli.accumulate(Hr);
          RefSCMatrix dH = x->coefs().kit()->matrix(x->coefs().rowdim(), ribs->coefs().rowdim());
          dH->convert(Hpauli);
          Hr = 0;
          Hpauli = 0;
          // transform to MO basis, transpose so that x dimension is coldim
          RefSCMatrix dH_mo = x->coefs().t() * dH * ribs->coefs();
          dH_mo = dH_mo.t();

          // Make new space hJ - dH
          std::string id = x->id();
          id += "_hJp(";
          id += ribs->id();
          id += ")";
          ExEnv::out0() << "id = " << id << endl;
          std::string name = "(hJp)-weighted space";
#if 1
          hJnr[spin] = new OrbitalSpace(id, name, hJ_x_P, hJ_x_P->coefs() - ribs->coefs() * dH_mo, hJ_x_P->basis());
#else
          hJnr[spin] = new OrbitalSpace(id, name, hJ_x_P, hJ_x_P->coefs(), hJ_x_P->basis());
#endif

          this->r12world()->world()->tfactory()->orbital_registry()->add(make_keyspace_pair(hJnr[s]));
        }
        else { // use pure h+J
          hJnr[spin] = hJ_x_P;
        }

      }
      if (Kr[Beta] == 0) Kr[Beta] = Kr[Alpha];
      if (hJnr[Beta] == 0) hJnr[Beta] = hJnr[Alpha];
    }


    for(int s=0; s<nspincases2(); s++) {
	const SpinCase2 spincase2 = static_cast<SpinCase2>(s);
	const SpinCase1 spin1 = case1(spincase2);
	const SpinCase1 spin2 = case2(spincase2);

	Ref<OrbitalSpace> occ1 = occ(spin1);
	Ref<OrbitalSpace> occ2 = occ(spin2);
	Ref<OrbitalSpace> orbs1 = orbs(spin1);
	Ref<OrbitalSpace> orbs2 = orbs(spin2);
	Ref<OrbitalSpace> GGspace1 = GGspace(spin1);
	Ref<OrbitalSpace> GGspace2 = GGspace(spin2);
	Ref<OrbitalSpace> vir1 = vir(spin1);
	Ref<OrbitalSpace> vir2 = vir(spin2);
	bool empty_vir_space = vir1->rank()==0 || vir2->rank()==0;

	// make sure that I have electrons for both spins
	if (GGspace1->rank() == 0 || GGspace2->rank() == 0)
	  continue;


#if INCLUDE_Q
	// if can only use 1 RI index, h+J can be resolved by the OBS
	Ref<OrbitalSpace> hj_x1, hj_x2;
	if (maxnabs > 0) {
	    hj_x1 = hj_x_P(spin1);
	    hj_x2 = hj_x_P(spin2);
	}
	else {
	    hj_x1 = hj_x_p(spin1);
	    hj_x2 = hj_x_p(spin2);
	}
    // if I'm relativistic, all relativistic terms are treated just as K, use hJ here that uses nonrelativistic core
    if (this->dk() > 0) {
      hj_x1 = hJnr[spin1];
      hj_x2 = hJnr[spin2];
    }
	std::string Qlabel = prepend_spincase(spincase2,"Q(C) intermediate");
	Timer tim_Q(Qlabel);
	ExEnv::out0() << endl << indent
		      << "Entered " << Qlabel << " evaluator" << endl;
	ExEnv::out0() << incindent;

	// compute Q = F12^2 (note F2_only = true in compute_X_ calls)
	RefSCMatrix Q;
	compute_X_(Q,spincase2,GGspace1,GGspace2,
		   GGspace1,hj_x2,true);
	if (GGspace1 != GGspace2) {
	    compute_X_(Q,spincase2,GGspace1,GGspace2,
		       hj_x1,GGspace2,true);
	}
	else {
	    Q.scale(2.0);
	    if (spincase2 == AlphaBeta) {
		symmetrize<false>(Q,Q,GGspace1,GGspace1);
	    }
	}

	ExEnv::out0() << decindent;
	ExEnv::out0() << indent << "Exited " << Qlabel << " evaluator" << endl;
	tim_Q.exit();

	if (debug_ >= DefaultPrintThresholds::mostO4 || DEBUG_PRINT_ALL_B_CONTRIBUTIONS) {
	  globally_sum_intermeds_();
	  std::string label = prepend_spincase(spincase2,"Q(C) contribution");
	  ExEnv::out0() << indent << __FILE__ << ": "<<__LINE__<<"\n";
      Q.print(label.c_str());
	}
	B_[s].accumulate(Q); Q = 0;
#endif // INCLUDE_Q

#if INCLUDE_P
	    // compute P

	    std::string Plabel = prepend_spincase(spincase2,"P(C) intermediate");
	    Timer tim_P(Plabel);
	    ExEnv::out0() << endl << indent
			  << "Entered " << Plabel << " evaluator" << endl;
	    ExEnv::out0() << incindent;

	    Ref<OrbitalSpace> ribs1, ribs2;
	    if (abs_eq_obs) {
		ribs1 = orbs1;
		ribs2 = orbs2;
	    }
	    else {
		ribs1 = r12world()->ribs_space();
		ribs2 = r12world()->ribs_space();
	    }
	    RefSCMatrix P;

#if INCLUDE_P_PKP
	    // R_klPQ K_QR R_PRij is included with projectors 2 and 3
	    {
		Ref<OrbitalSpace> kribs1, kribs2;
		if (abs_eq_obs) {
		    kribs1 = K_p_p(spin1);
		    kribs2 = K_p_p(spin2);
		}
		else {
		    kribs1 = K_P_P(spin1);
		    kribs2 = K_P_P(spin2);
		}

		// if I'm relativistic, treat all relativistic terms just as K, reassign kribs1 and kribs2
		if (this->dk() > 0) {
		  kribs1 = Kr[spin1];
          kribs2 = Kr[spin2];
		}

		compute_FxF_(P,spincase2,
			     GGspace1,GGspace2,
			     GGspace1,GGspace2,
			     ribs1,ribs2,
			     ribs1,ribs2,
			     kribs1,kribs2);

                if (debug_ >= DefaultPrintThresholds::allO4 || DEBUG_PRINT_ALL_B_CONTRIBUTIONS)
                  P.print("P(incl R_klPQ K_QR R_PRij)");
	    }
#endif // INCLUDE_P_PKP

	    // in this special case PFP and pFp terms can be combined
	    if (ansatz()->projector() == R12Technology::Projector_2 && abs_eq_obs && vbs_eq_obs) {
#if INCLUDE_P_PFP && INCLUDE_P_pFp
	      {
	        const Ref<OrbitalSpace> forbs1 = F_p_p(spin1);
	        const Ref<OrbitalSpace> forbs2 = F_p_p(spin2);
	        // R_klpr F_pq R_qrij
	        compute_FxF_(P,spincase2,
	                     GGspace1,GGspace2,
	                     GGspace1,GGspace2,
	                     orbs1,orbs2,
	                     orbs1,orbs2,
	                     forbs1,forbs2);
	        if (debug_ >= DefaultPrintThresholds::allO4 || DEBUG_PRINT_ALL_B_CONTRIBUTIONS)
	          P.print("P(incl R_klpr F_pq R_qrij)");
	      }
#endif // INCLUDE_P_PFP
	    }
	    else {

	    if (ansatz()->projector() == R12Technology::Projector_2) {
#if INCLUDE_P_PFP
		{
		    Ref<OrbitalSpace> fribs1,fribs2;
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
				 GGspace1,GGspace2,
				 GGspace1,GGspace2,
				 occ1,occ2,
				 ribs1,ribs2,
				 fribs1,fribs2);
                    if (debug_ >= DefaultPrintThresholds::allO4 || DEBUG_PRINT_ALL_B_CONTRIBUTIONS)
                      P.print("P(incl R_klPm F_PQ R_Qmij)");
		}
#endif // INCLUDE_P_PFP
	    }

	    {
#if INCLUDE_P_pFp
		// the form of the pFp term depends on the projector
		Ref<OrbitalSpace> z1 = (ansatz()->projector() == R12Technology::Projector_2 ? vir1 : orbs1);
		Ref<OrbitalSpace> z2 = (ansatz()->projector() == R12Technology::Projector_2 ? vir2 : orbs2);

		if (!empty_vir_space || ansatz()->projector() == R12Technology::Projector_3) {
		    // and on whether VBS==OBS
		    if (vbs_eq_obs) {
			Ref<OrbitalSpace> forbs1, forbs2;
			forbs1 = F_p_p(spin1);
			forbs2 = F_p_p(spin2);
			// R_klpa F_pq R_qaij
			compute_FxF_(P,spincase2,
				     GGspace1,GGspace2,
				     GGspace1,GGspace2,
				     z1,z2,
				     orbs1,orbs2,
				     forbs1,forbs2);
                        if (debug_ >= DefaultPrintThresholds::allO4 || DEBUG_PRINT_ALL_B_CONTRIBUTIONS)
                          P.print("P(incl R_klpa F_pq R_qaij)");
		    }
		    else {
			// else it's a sum of
			// R_klma F_mn R_naij
			{
			    Ref<OrbitalSpace> x1 = occ1;
			    Ref<OrbitalSpace> x2 = occ2;
			    Ref<OrbitalSpace> fx1 = F_m_m(spin1);
			    Ref<OrbitalSpace> fx2 = F_m_m(spin2);
			    compute_FxF_(P,spincase2,
					 GGspace1,GGspace2,
					 GGspace1,GGspace2,
					 z1,z2,
					 x1,x2,
					 fx1,fx2);
                            if (debug_ >= DefaultPrintThresholds::allO4 || DEBUG_PRINT_ALL_B_CONTRIBUTIONS)
                              P.print("P(incl R_klma F_mn R_naij)");
			}
			// R_klca F_cd R_daij
			{
			    Ref<OrbitalSpace> x1 = vir1;
			    Ref<OrbitalSpace> x2 = vir2;
			    Ref<OrbitalSpace> fx1 = F_a_a(spin1);
			    Ref<OrbitalSpace> fx2 = F_a_a(spin2);
			    compute_FxF_(P,spincase2,
					 GGspace1,GGspace2,
					 GGspace1,GGspace2,
					 z1,z2,
					 x1,x2,
					 fx1,fx2);
                            if (debug_ >= DefaultPrintThresholds::allO4 || DEBUG_PRINT_ALL_B_CONTRIBUTIONS)
                              P.print("P(incl R_klca F_cd R_daij)");
			}
			// R_klca F_cm R_maij + R_klma F_mc R_caij
			{
			    Ref<OrbitalSpace> x1 = occ1;
			    Ref<OrbitalSpace> x2 = occ2;
			    Ref<OrbitalSpace> fx1 = F_m_a(spin1);
			    Ref<OrbitalSpace> fx2 = F_m_a(spin2);
			    RefSCMatrix Ptmp;
			    compute_FxF_(Ptmp,spincase2,
					 GGspace1,GGspace2,
					 GGspace1,GGspace2,
					 z1,z2,
					 x1,x2,
					 fx1,fx2);
                            // bra-ket symmetrization will take care of
                            // the R_klma F_mc R_caij term
			    Ptmp.scale(2.0);
			    P.accumulate(Ptmp);
                            if (debug_ >= DefaultPrintThresholds::allO4 || DEBUG_PRINT_ALL_B_CONTRIBUTIONS)
                              P.print("P(incl 2 R_klca F_cm R_maij");
			}
		    } // VBS != OBS
		} // pFp contributes
#endif // INCLUDE_P_pFp
	    }
	    } // end of PFP + pFp code for OBS != ABS


	    if (!abs_eq_obs) {

		Ref<OrbitalSpace> cabs1 = r12world()->cabs_space(spin1);
		Ref<OrbitalSpace> cabs2 = r12world()->cabs_space(spin2);

		if (ansatz()->projector() == R12Technology::Projector_2) {
#if INCLUDE_P_mFP
                    // R_klmA F_mP R_PAij + R_klPA F_Pm R_mAij
		    {
			Ref<OrbitalSpace> focc1, focc2;
			focc1 = F_m_P(spin1);
			focc2 = F_m_P(spin2);
			RefSCMatrix Ptmp;
			compute_FxF_(Ptmp,spincase2,
				     GGspace1,GGspace2,
				     GGspace1,GGspace2,
				     cabs1,cabs2,
				     occ1,occ2,
				     focc1,focc2);
                        // bra-ket symmetrization will take care of
                        // the R_klma F_mc R_caij term
                        Ptmp.scale(2.0);
			P.accumulate(Ptmp);
                        if (debug_ >= DefaultPrintThresholds::allO4 || DEBUG_PRINT_ALL_B_CONTRIBUTIONS)
                          P.print("P(incl 2 R_klmA F_mP R_PAij)");
		    }
#endif // INCLUDE_P_mFP
		}

		{
#if INCLUDE_P_pFA
		    // the form of the pFp term depends on the projector
		    Ref<OrbitalSpace> z1 = (ansatz()->projector() == R12Technology::Projector_2 ? vir1 : orbs1);
		    Ref<OrbitalSpace> z2 = (ansatz()->projector() == R12Technology::Projector_2 ? vir2 : orbs2);

		    // R_klpz F_pA R_Azij
		    if (!empty_vir_space || ansatz()->projector() == R12Technology::Projector_3) {
			// and on whether VBS==OBS
			if (vbs_eq_obs) {
			    Ref<OrbitalSpace> forbs1 = F_p_A(spin1);
			    Ref<OrbitalSpace> forbs2 = F_p_A(spin2);
			    RefSCMatrix Ptmp;
			    compute_FxF_(Ptmp,spincase2,
					 GGspace1,GGspace2,
					 GGspace1,GGspace2,
					 z1,z2,
					 orbs1,orbs2,
					 forbs1,forbs2);
			    Ptmp.scale(2.0);
			    P.accumulate(Ptmp);
			} // VBS == OBS
			else { // VBS != OBS

			    RefSCMatrix Ptmp;
			    {
				Ref<OrbitalSpace> x1 = occ1;
				Ref<OrbitalSpace> x2 = occ2;
				Ref<OrbitalSpace> fx1 = F_m_A(spin1);
				Ref<OrbitalSpace> fx2 = F_m_A(spin2);
				compute_FxF_(Ptmp,spincase2,
					     GGspace1,GGspace2,
					     GGspace1,GGspace2,
					     z1,z2,
					     x1,x2,
					     fx1,fx2);
                                if (debug_ >= DefaultPrintThresholds::allO4 || DEBUG_PRINT_ALL_B_CONTRIBUTIONS)
                                  Ptmp.print("P(R_klmz F_mA R_Azij)");
			    }
			    {
				Ref<OrbitalSpace> x1 = vir1;
				Ref<OrbitalSpace> x2 = vir2;
				Ref<OrbitalSpace> fx1 = F_a_A(spin1);
				Ref<OrbitalSpace> fx2 = F_a_A(spin2);
				compute_FxF_(Ptmp,spincase2,
					     GGspace1,GGspace2,
					     GGspace1,GGspace2,
					     z1,z2,
					     x1,x2,
					     fx1,fx2);
                                if (debug_ >= DefaultPrintThresholds::allO4 || DEBUG_PRINT_ALL_B_CONTRIBUTIONS)
                                  Ptmp.print("P(R_klaz F_aA R_Azij)");
			    }
			    Ptmp.scale(2.0);
			    P.accumulate(Ptmp);

			} // VBS != OBS
                        if (debug_ >= DefaultPrintThresholds::allO4 || DEBUG_PRINT_ALL_B_CONTRIBUTIONS)
                          P.print("P(incl R_klpz F_pA R_Azij)");
		    } // pFA contributes
#endif // INCLUDE_P_pFA
		}

		P.scale(-1.0);

		if (ansatz()->projector() == R12Technology::Projector_2) {
#if INCLUDE_P_mFm
		    {
			Ref<OrbitalSpace> focc1 = F_m_m(spin1);
			Ref<OrbitalSpace> focc2 = F_m_m(spin2);
			// R_klmA F_mn R_nAij
			compute_FxF_(P,spincase2,
				     GGspace1,GGspace2,
				     GGspace1,GGspace2,
				     cabs1,cabs2,
				     occ1,occ2,
				     focc1,focc2);
                        if (debug_ >= DefaultPrintThresholds::allO4 || DEBUG_PRINT_ALL_B_CONTRIBUTIONS)
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
	    tim_P.exit();

	    if (debug_ >= DefaultPrintThresholds::mostO4 || DEBUG_PRINT_ALL_B_CONTRIBUTIONS) {
	      std::string label = prepend_spincase(spincase2,"P(C) contribution");
	      ExEnv::out0() << indent << __FILE__ << ": "<<__LINE__<<"\n";
          P.print(label.c_str());
        }

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

	tim_B_app_C.exit();
    }

////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
