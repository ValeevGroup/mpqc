//
// Xgen.cc
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
R12IntEval::compute_X_(RefSCMatrix& X,
                       SpinCase2 spincase2,
                       const Ref<MOIndexSpace>& bra1,
                       const Ref<MOIndexSpace>& bra2,
                       const Ref<MOIndexSpace>& ket1,
                       const Ref<MOIndexSpace>& ket2)
{
  const bool abs_eq_obs = r12info()->basis()->equiv(r12info()->basis_ri());
  const bool part1_equiv_part2 = (bra1 == bra2 && ket1 == ket2);
  
  // Check semantics
  bool correct_semantics = (spincase2 != AlphaBeta && part1_equiv_part2) ||
                           (spincase2 == AlphaBeta);
  if (!correct_semantics)
    throw ProgrammingError("R12IntEval::compute_X_() -- incorrect call semantics",__FILE__,__LINE__);
  
  // check number of ABS indices
  const Ref<GaussianBasisSet> abs = r12info()->basis_ri();
  const unsigned int nabs_in_bra1 = (bra1->basis() == abs);
  const unsigned int nabs_in_bra2 = (bra2->basis() == abs);
  const unsigned int nabs_in_ket1 = (ket1->basis() == abs);
  const unsigned int nabs_in_ket2 = (ket2->basis() == abs);
  const unsigned int nabs_in_bra = nabs_in_bra1 + nabs_in_bra2;
  const unsigned int nabs_in_ket = nabs_in_ket1 + nabs_in_ket2;
  const unsigned int maxnabs = r12info()->maxnabs();
  if (nabs_in_bra > maxnabs ||
      nabs_in_ket > maxnabs) {
    throw ProgrammingError("R12IntEval::compute_X_() -- maxnabs is exceeded",__FILE__,__LINE__);
  }
  const unsigned int nabs = max(nabs_in_bra,nabs_in_ket);
  // check if RI needs to be done in ABS
  const bool do_ri_in_abs = !abs_eq_obs && (maxnabs - nabs > 0);
  
  tim_enter("generic X intermediate");
  ExEnv::out0() << indent << "Entered generic X intermediate evaluator" << endl;
  ExEnv::out0() << incindent;
  
  SpinMOPairIter braiter(bra1,bra2,spincase2);
  SpinMOPairIter ketiter(ket1,ket2,spincase2);
  const unsigned int nbra = braiter.nij();
  const unsigned int nket = ketiter.nij();
  
  if (X.null()) {
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
    
    const SpinCase1 spin1 = case1(spincase2);
    const SpinCase1 spin2 = case2(spincase2);
    
    Ref<SingleRefInfo> refinfo = r12info()->refinfo();
    Ref<MOIndexSpace> occ1 = refinfo->occ(spin1);
    Ref<MOIndexSpace> occ2 = refinfo->occ(spin2);
    Ref<MOIndexSpace> vir1 = vir(spin1);
    Ref<MOIndexSpace> vir2 = vir(spin2);
    Ref<MOIndexSpace> orbs1 = refinfo->orbs(spin1);
    Ref<MOIndexSpace> orbs2 = refinfo->orbs(spin2);
    // if orbs1 and orbs2 have different rank -- something is TERRIBLY wrong
    if (orbs1->rank() != orbs2->rank())
      throw ProgrammingError("R12IntEval::compute_X_() -- orbs1 and orbs2 have different ranks",__FILE__,__LINE__);
    const unsigned int nobs = orbs1->rank();
    
    Ref<R12IntEval> thisref(this);
    // (i p |j p) tforms
    std::vector<  Ref<TwoBodyMOIntsTransform> > tforms_ipjp;
    {
      NewTransformCreator tform_creator(thisref,bra1,orbs1,bra2,orbs2,true);
      fill_container(tform_creator,tforms_ipjp);
    }
    // (k p |l p) tforms
    std::vector<  Ref<TwoBodyMOIntsTransform> > tforms_kplp;
    {
      NewTransformCreator tform_creator(thisref,ket1,orbs1,ket2,orbs2,true);
      fill_container(tform_creator,tforms_kplp);
    }
    
    //
    // F12^2 contribution depends on the type of correlation factor
    //
    enum {r12corrfactor, g12corrfactor} corrfac;
    Ref<LinearR12::R12CorrelationFactor> r12corrptr; r12corrptr << r12info()->corrfactor();
    Ref<LinearR12::G12CorrelationFactor> g12corrptr; g12corrptr << r12info()->corrfactor();
    if (r12corrptr.nonnull()) corrfac = r12corrfactor;
    if (g12corrptr.nonnull()) corrfac = g12corrfactor;
    
    switch (corrfac) {
      case r12corrfactor:  // R12^2 reduces to one-electron integrals
      {
        RefSCMatrix R2_ijkl = compute_r2_(bra1,bra2,ket1,ket2);
        if (spincase2 == AlphaBeta) {
          X.accumulate(R2_ijkl);
        }
        else {
          antisymmetrize(X,R2_ijkl,bra1,ket1,true);
        }
        R2_ijkl = 0;
      }
      break;
      
      case g12corrfactor: // G12^2 involves two-electron integrals
      {
        // (i k |j l) tforms
        std::vector<  Ref<TwoBodyMOIntsTransform> > tforms_ikjl;
        {
          NewTransformCreator tform_creator(thisref,bra1,ket1,bra2,ket2,true,true);
          fill_container(tform_creator,tforms_ikjl);
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
    
    // ABS and CABS method differ by the TwoParticleContraction
    using LinearR12::TwoParticleContraction;
    using LinearR12::ABS_OBS_Contraction;
    using LinearR12::CABS_OBS_Contraction;
    using LinearR12::Direct_Contraction;
    const LinearR12::ABSMethod absmethod = r12info()->abs_method();
    Ref<TwoParticleContraction> contract_pp;
    if ((absmethod == LinearR12::ABS_ABS ||
         absmethod == LinearR12::ABS_ABSPlus) && do_ri_in_abs)
      contract_pp = new ABS_OBS_Contraction(nobs,
                                            occ1->rank(),
                                            occ2->rank());
    else
      contract_pp = new CABS_OBS_Contraction(nobs);
    
    // compute ABS/CABS contraction for <ij|F12|pp> . <kl|F12|pp>
    contract_tbint_tensor<
      ManyBodyTensors::I_to_T,
      ManyBodyTensors::I_to_T,
      ManyBodyTensors::I_to_T,
      true,true,false>
      (
        X, corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_f12(),
        bra1, bra2,
        orbs1, orbs2,
        ket1, ket2,
        orbs1, orbs2,
        contract_pp,
        spincase2!=AlphaBeta, tforms_ipjp, tforms_kplp
      );
    
    if (do_ri_in_abs) {
      Ref<MOIndexSpace> ribs1 = r12info()->ribs_space(spin1);      Ref<MOIndexSpace> ribs2 = r12info()->ribs_space(spin2);
      
      // (i m |j a') tforms
      std::vector<  Ref<TwoBodyMOIntsTransform> > tforms_imjA;
      {
        NewTransformCreator tform_creator(thisref,bra1,occ1,bra2,ribs2,true);
        fill_container(tform_creator,tforms_imjA);
      }
      // (k m |l a') tforms
      std::vector<  Ref<TwoBodyMOIntsTransform> > tforms_kmlA;
      {
        NewTransformCreator tform_creator(thisref,ket1,occ1,ket2,ribs2,true);
        fill_container(tform_creator,tforms_kmlA);
      }
      
      Ref<TwoParticleContraction> dircontract_mA =
        new Direct_Contraction(occ1->rank(),ribs2->rank(),1.0);
      
      // compute contraction -1 * <ij|F12|m a'> . <kl|F12|m a'>
      contract_tbint_tensor<
        ManyBodyTensors::I_to_T,
        ManyBodyTensors::I_to_T,
        ManyBodyTensors::I_to_mT,
        true,true,false>
        (
          X, corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_f12(),
          bra1, bra2,
          occ1, ribs2,
          ket1, ket2,
          occ1, ribs2,
          dircontract_mA,
          spincase2!=AlphaBeta, tforms_imjA, tforms_kmlA
        );
      
      if (!part1_equiv_part2) {

        Ref<MOIndexSpace> ribs1 = r12info()->ribs_space(spin1);

        // (i a' |j m) tforms
        std::vector<  Ref<TwoBodyMOIntsTransform> > tforms_iAjm;
        {
          NewTransformCreator tform_creator(thisref,bra1,ribs1,bra2,occ2,true);
          fill_container(tform_creator,tforms_iAjm);
        }
        // (k a' |l m) tforms
        std::vector<  Ref<TwoBodyMOIntsTransform> > tforms_kAlm;
        {
          NewTransformCreator tform_creator(thisref,ket1,ribs1,ket2,occ2,true);
          fill_container(tform_creator,tforms_kAlm);
        }
        
        Ref<TwoParticleContraction> dircontract_Am =
        new Direct_Contraction(ribs1->rank(),occ2->rank(),1.0);
        
        // compute contraction -1 * <ij|F12|a'm> . <kl|F12|a'm>
        contract_tbint_tensor<
          ManyBodyTensors::I_to_T,
          ManyBodyTensors::I_to_T,
          ManyBodyTensors::I_to_mT,
          true,true,false>
          (
            X, corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_f12(),
            bra1, bra2,
            ribs1, occ2,
            ket1, ket2,
            ribs1, occ2,
            dircontract_Am,
            spincase2!=AlphaBeta, tforms_iAjm, tforms_kAlm
          );
      }
    }
    
    // make particles equivalent, if necessary
    if (part1_equiv_part2) {
      symmetrize<false>(X,X,bra1,ket1);
    }

    if (debug_ > 1) {
      std::string label = prepend_spincase(spincase2,"Generic X");
      X.print(label.c_str());
    }
    
    // Bra-Ket symmetrize
    X.scale(0.5);
    RefSCMatrix X_t = X.t();
    X.accumulate(X_t);
  
  globally_sum_scmatrix_(X);
  
  ExEnv::out0() << decindent;
  ExEnv::out0() << indent << "Exited generic X intermediate evaluator" << endl;

  tim_exit("generic X intermediate");
}

////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
