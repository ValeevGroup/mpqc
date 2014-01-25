//
// Xgen.cc
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
#include <sstream>
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
#include <chemistry/qc/mbptr12/contract_tbint_tensor.h>
#include <chemistry/qc/mbptr12/twoparticlecontraction.h>
#include <chemistry/qc/lcao/utils.h>
#include <chemistry/qc/lcao/utils.impl.h>
#include <util/misc/print.h>

using namespace std;
using namespace sc;

#define SYMMETRIZE 1
#define PRINT_COMPONENTS 0

namespace {
  void print_component(const RefSCMatrix& X, const std::string& label) {
    std::ostringstream oss;
    oss << "Component of X: " << label << std::endl;
    X.print(oss.str().c_str());
    X.assign(0.0);
  }
};


void
R12IntEval::compute_X_(RefSCMatrix& X,
                       SpinCase2 spincase2,
                       const Ref<OrbitalSpace>& bra1,
                       const Ref<OrbitalSpace>& bra2,
                       const Ref<OrbitalSpace>& ket1,
                       const Ref<OrbitalSpace>& ket2,
                       bool F2_only)
{
    using namespace sc::mbptr12;

    // equations depend on whether VBS == OBS ..
    const bool vbs_eq_obs = r12world()->obs_eq_vbs();
    // .. and ABS == OBS
    const bool abs_eq_obs = r12world()->obs_eq_ribs();
    // if particle 1 and 2 are equivalent, can use permutational symmetry
    const bool part1_equiv_part2 = (bra1 == bra2 && ket1 == ket2);

    //
    // check number of ABS indices, make sure user-imposed maxnabs not exceeded
    //
    const Ref<GaussianBasisSet> abs = r12world()->basis_ri();
    const unsigned int nabs_in_bra1 = abs_eq_obs ? 0 : (bra1->basis() == abs);
    const unsigned int nabs_in_bra2 = abs_eq_obs ? 0 : (bra2->basis() == abs);
    const unsigned int nabs_in_ket1 = abs_eq_obs ? 0 : (ket1->basis() == abs);
    const unsigned int nabs_in_ket2 = abs_eq_obs ? 0 : (ket2->basis() == abs);
    const unsigned int nabs_in_bra = nabs_in_bra1 + nabs_in_bra2;
    const unsigned int nabs_in_ket = nabs_in_ket1 + nabs_in_ket2;
    const unsigned int maxnabs = r12world()->r12tech()->maxnabs();
    if (!F2_only && (nabs_in_bra > maxnabs ||
	             nabs_in_ket > maxnabs)) {
	throw ProgrammingError("R12IntEval::compute_X_() -- maxnabs is exceeded",__FILE__,__LINE__);
    }
    const unsigned int nabs = max(nabs_in_bra,nabs_in_ket);
    // check if RI needs to be done in ABS
    const bool do_ri_in_abs = !abs_eq_obs && (maxnabs - nabs > 0);

    ////////////////////////////////
    // Game begins
    ////////////////////////////////
    Timer tim("generic X intermediate");
    ExEnv::out0() << indent << "Entered generic X intermediate evaluator" << endl;
    ExEnv::out0() << incindent;

    // geminal dimensions
    const unsigned int nf12 = corrfactor()->nfunctions();
    SpinMOPairIter braiter(bra1->rank(), (spincase2==AlphaBeta ? bra2 : bra1)->rank(), spincase2);
    SpinMOPairIter ketiter(ket1->rank(), (spincase2==AlphaBeta ? ket2 : ket1)->rank(), spincase2);
    const unsigned int nbra = nf12 * braiter.nij();
    const unsigned int nket = nf12 * ketiter.nij();

    // init the target X matrix
    if (X == 0) {
	// use the same matrix kit as the intermediates
	X = B_[AlphaBeta].kit()->matrix(new SCDimension(nbra),
					new SCDimension(nket));
	X.assign(0.0);
    }
    else {
	if (X.rowdim().n() != nbra)
	    throw ProgrammingError("R12IntEval::compute_X_() -- row dimension of the given X doesn't match given bra dimensions",__FILE__,__LINE__);
	if (X.coldim().n() != nket)
	    throw ProgrammingError("R12IntEval::compute_X_() -- column dimension of the given X doesn't match given ket dimensions",__FILE__,__LINE__);
    }

    // get orbital spaces and verify their sanity
    const SpinCase1 spin1 = case1(spincase2);
    const SpinCase1 spin2 = case2(spincase2);
    Ref<OrbitalSpace> occ1 = occ(spin1);
    Ref<OrbitalSpace> occ2 = occ(spin2);
    Ref<OrbitalSpace> orbs1 = orbs(spin1);
    Ref<OrbitalSpace> orbs2 = orbs(spin2);
    // if orbs1 and orbs2 have different rank -- something is TERRIBLY wrong
    if (orbs1->rank() != orbs2->rank())
	throw ProgrammingError("R12IntEval::compute_X_() -- orbs1 and orbs2 have different ranks",__FILE__,__LINE__);
    const unsigned int nobs = orbs1->rank();


    // *** Working equations for ansatz 2 ***
    // if VBS==OBS: X_{ij}^{kl} = RR_{ij}^{kl} - 1/2 R_{ij}^{pq} * R_{pq}^{kl} - R_{ij}^{ma'} * R_{ma'}^{kl}
    // if VBS!=OBS: X_{ij}^{kl} = RR_{ij}^{kl} - 1/2 R_{ij}^{mn} * R_{mn}^{kl} - 1/2 R_{ij}^{ab} * R_{ab}^{kl}
    //                            - R_{ij}^{ma} * R_{ma}^{kl} - R_{ij}^{ma'} * R_{ma'}^{kl}
    // hence the diagonal (RR) and the RI parts of X do not depend on whether VBS==OBS.
    // in ansatz 3 the RI terms do not contribute at all

    //
    // F12^2 contribution depends on the type of correlation factor
    //
    enum {r12corrfactor, g12corrfactor, gg12corrfactor} corrfac;
    Ref<R12Technology::R12CorrelationFactor> r12corrptr; r12corrptr << corrfactor();
    Ref<R12Technology::G12CorrelationFactor> g12corrptr; g12corrptr << corrfactor();
    Ref<R12Technology::G12NCCorrelationFactor> g12nccorrptr; g12nccorrptr << corrfactor();
    if (r12corrptr) corrfac = r12corrfactor;
    if (g12corrptr) corrfac = g12corrfactor;
    if (g12nccorrptr) corrfac = g12corrfactor;

    switch (corrfac) {
    case r12corrfactor:  // R12^2 reduces to one-electron integrals
    {
        RefSCMatrix R2_ijkl = compute_r2_(bra1,bra2,ket1,ket2);
        if (spincase2 == AlphaBeta) {
	    X.accumulate(R2_ijkl);
        }
        else {
	    if (!part1_equiv_part2 /* && spincase2 != AlphaBeta */ )
		symmetrize<false>(R2_ijkl,R2_ijkl,bra1,ket1);
	    antisymmetrize(X,R2_ijkl,bra1,ket1,true);
        }
    }
    break;

    // G12^2 involves two-electron integrals
    case g12corrfactor:
    case gg12corrfactor:
    {
        // (i k |j l) tforms
        std::vector<std::string> tforms_ikjl;
        {
	    R12TwoBodyIntKeyCreator tformkey_creator(moints_runtime4(),bra1,ket1,bra2,ket2,
	                                             corrfactor(),true,true);
	    fill_container(tformkey_creator,tforms_ikjl);
        }
        compute_tbint_tensor<ManyBodyTensors::I_to_T,true,true>(
	    X, corrfactor()->tbint_type_f12f12(),
	    bra1, ket1, bra2, ket2, spincase2!=AlphaBeta,
	    tforms_ikjl
	    );
    }
    break;

    default:
	throw ProgrammingError("R12IntEval::compute_X_() -- unrecognized type of correlation factor",__FILE__,__LINE__);
    }

#if PRINT_COMPONENTS
    print_component(X,"F12^2");
#endif

    if (!F2_only) {

	if (vbs_eq_obs) {

	    Ref<TwoParticleContraction> contract_pp = new CABS_OBS_Contraction(nobs);

	    // (i p |j p) tforms
	    std::vector<std::string> tforms_ipjp;
	    {
		R12TwoBodyIntKeyCreator tformkey_creator(moints_runtime4(),bra1,orbs1,bra2,orbs2,
		                                         corrfactor(),true);
		fill_container(tformkey_creator,tforms_ipjp);
	    }
	    // (k p |l p) tforms
	    std::vector<std::string> tforms_kplp;
	    {
		R12TwoBodyIntKeyCreator tformkey_creator(moints_runtime4(),ket1,orbs1,ket2,orbs2,
                                                 corrfactor(),true);
		fill_container(tformkey_creator,tforms_kplp);
	    }


	    // compute ABS/CABS contraction for <ij|F12|pp> . <kl|F12|pp>
	    contract_tbint_tensor<ManyBodyTensors::I_to_T,ManyBodyTensors::I_to_T,ManyBodyTensors::I_to_T,true,true,false>
		(X, corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_f12(),
		 bra1, bra2, orbs1, orbs2,
		 ket1, ket2, orbs1, orbs2,
		 contract_pp, spincase2!=AlphaBeta, tforms_ipjp, tforms_kplp
		    );
#if PRINT_COMPONENTS
	    print_component(X,"<ij|pp>");
#endif
	}
	// VBS != OBS
	else {

	    const double asymm_contr_pfac = part1_equiv_part2 ? -2.0 : -1.0;
	    // (im|jn) contribution
	    {
		Ref<OrbitalSpace> cs1 = occ1;
		Ref<OrbitalSpace> cs2 = occ2;
		Ref<TwoParticleContraction> tpcontract = new Direct_Contraction(cs1->rank(),cs2->rank(),-1.0);
		std::vector<std::string> tforms_f12_ij;
		{
		    R12TwoBodyIntKeyCreator tformkey_creator(moints_runtime4(),bra1,cs1,bra2,cs2,
	                                                 corrfactor(),true);
		    fill_container(tformkey_creator,tforms_f12_ij);
		}
		std::vector<std::string> tforms_f12_kl;
		{
		    R12TwoBodyIntKeyCreator tformkey_creator(moints_runtime4(),ket1,cs1,ket2,cs2,
	                                                 corrfactor(),true);
		    fill_container(tformkey_creator,tforms_f12_kl);
		}
		contract_tbint_tensor<ManyBodyTensors::I_to_T,ManyBodyTensors::I_to_T,ManyBodyTensors::I_to_T,true,true,false>
		    (X, corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_f12(),
		     bra1, bra2, cs1, cs2,
		     ket1, ket2, cs1, cs2,
		     tpcontract, spincase2!=AlphaBeta, tforms_f12_ij, tforms_f12_kl);
	    } // (im|jn)
	    // (ia|jb) contribution
	    {
		Ref<OrbitalSpace> cs1 = vir_act(spin1);
		Ref<OrbitalSpace> cs2 = vir_act(spin2);
		Ref<TwoParticleContraction> tpcontract = new Direct_Contraction(cs1->rank(),cs2->rank(),-1.0);
		std::vector<std::string> tforms_f12_ij;
		{
		    R12TwoBodyIntKeyCreator tformkey_creator(moints_runtime4(),bra1,cs1,bra2,cs2,
	                                                 corrfactor(),true);
		    fill_container(tformkey_creator,tforms_f12_ij);
		}
		std::vector<std::string> tforms_f12_kl;
		{
		    R12TwoBodyIntKeyCreator tformkey_creator(moints_runtime4(),ket1,cs1,ket2,cs2,
	                                                 corrfactor(),true);
		    fill_container(tformkey_creator,tforms_f12_kl);
		}
		contract_tbint_tensor<ManyBodyTensors::I_to_T,ManyBodyTensors::I_to_T,ManyBodyTensors::I_to_T,true,true,false>
		    (X, corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_f12(),
		     bra1, bra2, cs1, cs2,
		     ket1, ket2, cs1, cs2,
		     tpcontract, spincase2!=AlphaBeta, tforms_f12_ij, tforms_f12_kl);
	    } // (ia|jb)
	    // (im|ja) contribution
	    {
		Ref<OrbitalSpace> cs1 = occ1;
		Ref<OrbitalSpace> cs2 = vir_act(spin2);
		Ref<TwoParticleContraction> tpcontract = new Direct_Contraction(cs1->rank(),cs2->rank(),asymm_contr_pfac);
		std::vector<std::string> tforms_f12_ij;
		{
		    R12TwoBodyIntKeyCreator tformkey_creator(moints_runtime4(),bra1,cs1,bra2,cs2,
	                                                 corrfactor(),true);
		    fill_container(tformkey_creator,tforms_f12_ij);
		}
		std::vector<std::string> tforms_f12_kl;
		{
		    R12TwoBodyIntKeyCreator tformkey_creator(moints_runtime4(),ket1,cs1,ket2,cs2,
	                                                 corrfactor(),true);
		    fill_container(tformkey_creator,tforms_f12_kl);
		}
		contract_tbint_tensor<ManyBodyTensors::I_to_T,ManyBodyTensors::I_to_T,ManyBodyTensors::I_to_T,true,true,false>
		    (X, corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_f12(),
		     bra1, bra2, cs1, cs2,
		     ket1, ket2, cs1, cs2,
		     tpcontract, spincase2!=AlphaBeta, tforms_f12_ij, tforms_f12_kl);
	    } // (im|ja)
	    // (ia|jm) contribution
	    if (!part1_equiv_part2) {
		Ref<OrbitalSpace> cs1 = vir_act(spin1);
		Ref<OrbitalSpace> cs2 = occ2;
		Ref<TwoParticleContraction> tpcontract = new Direct_Contraction(cs1->rank(),cs2->rank(),-1.0);
		std::vector<std::string> tforms_f12_ij;
		{
		    R12TwoBodyIntKeyCreator tformkey_creator(moints_runtime4(),bra1,cs1,bra2,cs2,
	                                                 corrfactor(),true);
		    fill_container(tformkey_creator,tforms_f12_ij);
		}
		std::vector<std::string> tforms_f12_kl;
		{
		    R12TwoBodyIntKeyCreator tformkey_creator(moints_runtime4(),ket1,cs1,ket2,cs2,
	                                                 corrfactor(),true);
		    fill_container(tformkey_creator,tforms_f12_kl);
		}
		contract_tbint_tensor<ManyBodyTensors::I_to_T,ManyBodyTensors::I_to_T,ManyBodyTensors::I_to_T,true,true,false>
		    (X, corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_f12(),
		     bra1, bra2, cs1, cs2,
		     ket1, ket2, cs1, cs2,
		     tpcontract, spincase2!=AlphaBeta, tforms_f12_ij, tforms_f12_kl);
	    } // (ia|jm)

	} // VBS != OBS

	// These are only needed in ansatz 2
	if (ansatz()->projector() == R12Technology::Projector_2 && do_ri_in_abs) {
	    Ref<OrbitalSpace> ribs2 = r12world()->cabs_space(spin2);

	    // (i m |j a') tforms
	    std::vector<std::string> tforms_imjA;
	    {
		R12TwoBodyIntKeyCreator tformkey_creator(moints_runtime4(),bra1,occ1,bra2,ribs2,
                                                 corrfactor(),true);
		fill_container(tformkey_creator,tforms_imjA);
	    }
	    // (k m |l a') tforms
	    std::vector<std::string> tforms_kmlA;
	    {
		R12TwoBodyIntKeyCreator tformkey_creator(moints_runtime4(),ket1,occ1,ket2,ribs2,
                                                 corrfactor(),true);
		fill_container(tformkey_creator,tforms_kmlA);
	    }

	    const double perm_pfac = (part1_equiv_part2 ? 2.0 : 1.0);
	    Ref<TwoParticleContraction> dircontract_mA =
		new Direct_Contraction(occ1->rank(),ribs2->rank(),perm_pfac);

	    // compute contraction -1 * <ij|F12|m a'> . <kl|F12|m a'>
	    contract_tbint_tensor<ManyBodyTensors::I_to_T,ManyBodyTensors::I_to_T,ManyBodyTensors::I_to_mT,true,true,false>
		(
		    X, corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_f12(),
		    bra1, bra2, occ1, ribs2,
		    ket1, ket2, occ1, ribs2,
		    dircontract_mA, spincase2!=AlphaBeta, tforms_imjA, tforms_kmlA
		    );
#if PRINT_COMPONENTS
	    print_component(X,"<ij|ma'>");
#endif

	    if (!part1_equiv_part2) {

		Ref<OrbitalSpace> ribs1 = r12world()->cabs_space(spin1);

		// (i a' |j m) tforms
		std::vector<std::string> tforms_iAjm;
		{
		    R12TwoBodyIntKeyCreator tformkey_creator(moints_runtime4(),bra1,ribs1,bra2,occ2,
	                                                 corrfactor(),true);
		    fill_container(tformkey_creator,tforms_iAjm);
		}
		// (k a' |l m) tforms
		std::vector<std::string> tforms_kAlm;
		{
		    R12TwoBodyIntKeyCreator tformkey_creator(moints_runtime4(),ket1,ribs1,ket2,occ2,
	                                                 corrfactor(),true);
		    fill_container(tformkey_creator,tforms_kAlm);
		}

		Ref<TwoParticleContraction> dircontract_Am =
		    new Direct_Contraction(ribs1->rank(),occ2->rank(),1.0);

		// compute contraction -1 * <ij|F12|a'm> . <kl|F12|a'm>
		contract_tbint_tensor<ManyBodyTensors::I_to_T,ManyBodyTensors::I_to_T,ManyBodyTensors::I_to_mT,true,true,false>
		    (
			X, corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_f12(),
			bra1, bra2, ribs1, occ2,
			ket1, ket2, ribs1, occ2,
			dircontract_Am, spincase2!=AlphaBeta, tforms_iAjm, tforms_kAlm
			);
#if PRINT_COMPONENTS
		print_component(X,"<ij|a'm>");
#endif
	    }
	}

	if (debug_ >= DefaultPrintThresholds::mostO4) {
	    std::string label = prepend_spincase(spincase2,"generic X");
	    X.print(label.c_str());
	}

	// Bra-Ket symmetrize
	X.scale(0.5);
	RefSCMatrix X_t = X.t();
	X.accumulate(X_t);  X_t = 0;
    } // if (!F2_only)

    ////////////////////////////////
    // Game over
    ////////////////////////////////
    globally_sum_scmatrix_(X);
    ExEnv::out0() << decindent;
    ExEnv::out0() << indent << "Exited generic X intermediate evaluator" << endl;
}

////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
