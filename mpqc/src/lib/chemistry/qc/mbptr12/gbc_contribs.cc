//
// gbc_contribs.cc
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

#define INCLUDE_GBC1 1
#define INCLUDE_GBC2 1
#define COMPUTE_GBC1_AS_FXF 1
#define COMPUTE_GBC2_AS_X 1

void
R12IntEval::compute_B_gbc_()
{
  if (abs_method_ == LinearR12::ABS_ABS || abs_method_ == LinearR12::ABS_ABSPlus)
    throw std::runtime_error("R12IntEval::compute_B_gbc_1_() -- B(GBC1) term can only be computed using a CABS (or CABS+) approach");
  
  if (evaluated_)
    return;
  
  tim_enter("B(GBC) intermediate");
  ExEnv::out0() << endl << indent
  << "Entered B(GBC) intermediate evaluator" << endl;
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
    Ref<MOIndexSpace> ribs1 = r12info()->ribs_space(spin1);
    Ref<MOIndexSpace> ribs2 = r12info()->ribs_space(spin2);
    Ref<MOIndexSpace> occ1_act = occ_act(spin1);
    Ref<MOIndexSpace> occ2_act = occ_act(spin2);
    Ref<MOIndexSpace> vir1 = vir(spin1);
    Ref<MOIndexSpace> vir2 = vir(spin2);
    Ref<MOIndexSpace> focc1 = focc(spin1);
    Ref<MOIndexSpace> focc2 = focc(spin2);
    Ref<MOIndexSpace> focc1_act = focc_act(spin1);
    Ref<MOIndexSpace> focc2_act = focc_act(spin2);
    
    RefSCMatrix B_gbc1 = B_[s].clone(); B_gbc1.assign(0.0);
    RefSCMatrix B_gbc2 = B_[s].clone(); B_gbc2.assign(0.0);
    
#if !COMPUTE_GBC1_AS_FXF || !COMPUTE_GBC2_AS_X
    using namespace sc::LinearR12;
    Ref<TwoParticleContraction> dircontract_mA =
      new Direct_Contraction(occ1->rank(),ribs2->rank(),1.0);
    Ref<TwoParticleContraction> dircontract_ma =
      new Direct_Contraction(occ1->rank(),vir2->rank(),1.0);
    Ref<TwoParticleContraction> dircontract_pp =
      new Direct_Contraction(orbs1->rank(),orbs2->rank(),1.0);
    
    // (i p |i p) tforms are ready
    std::vector<  Ref<TwoBodyMOIntsTransform> > tforms_ipip;
    Ref<R12IntEval> thisref(this);
    {
      NamedTransformCreator tform_creator(thisref,occ1_act,orbs1,occ2_act,orbs2,true);
      fill_container(tform_creator,tforms_ipip);
    }
    
    // (i m |i a') tforms are ready
    std::vector<  Ref<TwoBodyMOIntsTransform> > tforms_imiA;
    {
      NewTransformCreator tform_creator(thisref,occ1_act,occ1,occ2_act,ribs2,true);
      fill_container(tform_creator,tforms_imiA);
    }
    
    // (i m_F|i a) tforms are new
    std::vector<  Ref<TwoBodyMOIntsTransform> > tforms_iMia;
    {
      NewTransformCreator tform_creator(thisref,occ1_act,focc1,occ2_act,vir2,true);
      fill_container(tform_creator,tforms_iMia);
    }
    
    // (i m_F|i a') tforms are new
    std::vector<  Ref<TwoBodyMOIntsTransform> > tforms_iMiA;
    {
      NewTransformCreator tform_creator(thisref,occ1_act,focc1,occ2_act,ribs2,true);
      fill_container(tform_creator,tforms_iMiA);
    }
    
    // (i m|i_F a') tforms are new
    std::vector<  Ref<TwoBodyMOIntsTransform> > tforms_imIA;
    {
      NewTransformCreator tform_creator(thisref,occ1_act,occ1,focc2_act,ribs2,true);
      fill_container(tform_creator,tforms_imIA);
    }
    
    // (i_F m|i a') tforms are new
    std::vector<  Ref<TwoBodyMOIntsTransform> > tforms_ImiA;
    {
      NewTransformCreator tform_creator(thisref,focc1_act,occ1,occ2_act,ribs2,true);
      fill_container(tform_creator,tforms_ImiA);
    }
    
    // (i p|i_F p) tforms are new
    std::vector<  Ref<TwoBodyMOIntsTransform> > tforms_ipIp;
    {
      NewTransformCreator tform_creator(thisref,occ1_act,orbs1,focc2_act,orbs2,true);
      fill_container(tform_creator,tforms_ipIp);
    }
#endif
    
#if INCLUDE_GBC1

#if !COMPUTE_GBC1_AS_FXF
    
    //
    // GBC1 contribution
    //
    
    // compute contraction <ii|F12|m a'> . <ii|F12|m_F a'>
    contract_tbint_tensor<
      ManyBodyTensors::I_to_T,
      ManyBodyTensors::I_to_T,
      ManyBodyTensors::I_to_T,
      true,true,false>
      (
        B_gbc1, corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_f12(),
        occ1_act, occ2_act,
        occ1, ribs2,
        occ1_act, occ2_act,
        focc1, ribs2,
        dircontract_mA,
        spincase2!=AlphaBeta, tforms_imiA, tforms_iMiA
      );
    
    // compute contraction <ii|F12|m a> . <ii|F12|m_F a>
    contract_tbint_tensor<
      ManyBodyTensors::I_to_T,
      ManyBodyTensors::I_to_T,
      ManyBodyTensors::I_to_T,
      true,true,false>
      (
        B_gbc1, corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_f12(),
        occ1_act, occ2_act,
        occ1, vir2,
        occ1_act, occ2_act,
        focc1, vir2,
        dircontract_ma,
        spincase2!=AlphaBeta, tforms_ipip, tforms_iMia
      );
#else // use compute_FxF_
    
    // R_klAb F_Am R_mbij
    compute_FxF_(B_gbc1,spincase2,
                 occ1_act,occ2_act,
                 occ1_act,occ2_act,
                 vir1,vir2,
                 occ1,occ2,
                 focc1,focc2);
    
    if (r12info()->maxnabs() >= 2) {
      // R_klAB F_Am R_mBij
      compute_FxF_(B_gbc1,spincase2,
                   occ1_act,occ2_act,
                   occ1_act,occ2_act,
                   ribs1,ribs2,
                   occ1,occ2,
                   focc1,focc2);
    }
    
    B_gbc1.scale(-1.0);
    
#endif // if use compute FxF
    
#endif // include GBC1
    
#if INCLUDE_GBC2
    
    //
    // GBC2 contribution
    //
    
#if !COMPUTE_GBC2_AS_X
    
    RefSCMatrix R2_ijkL = compute_r2_(occ1_act,occ2_act,occ1_act,focc2_act);
    B_gbc2.accumulate(R2_ijkL);
    
    //
    // Compute contribution X -= r_{ij}^{\alpha'm} r_{m\alpha'}^{k l_f}
    //                         + r_{ji}^{\alpha'm} r_{\alpha'm}^{k l_f}
    //
    
    // compute contraction -1 * <ii|F12|m a'> . <iI|F12|m a'>
    contract_tbint_tensor<
      ManyBodyTensors::I_to_T,
      ManyBodyTensors::I_to_T,
      ManyBodyTensors::I_to_mT,
      true,true,false>
      (
        B_gbc2, corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_f12(),
        occ1_act, occ2_act,
        occ1, ribs2,
        occ1_act, focc2_act,
        occ1, ribs2,
        dircontract_mA,
        spincase2!=AlphaBeta, tforms_imiA, tforms_imIA
      );
    
    // compute contraction -1 * <ii|F12|m a'> . <Ii|F12|m a'>
    contract_tbint_tensor<
      ManyBodyTensors::I_to_T,
      ManyBodyTensors::I_to_T,
      ManyBodyTensors::I_to_mT,
      true,true,false>
      (
        B_gbc2, corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_f12(),
        occ1_act, occ2_act,
        occ1, ribs2,
        focc1_act, occ2_act,
        occ1, ribs2,
        dircontract_mA,
        spincase2!=AlphaBeta, tforms_imiA, tforms_ImiA
      );
    
    //
    // Compute contribution X -= r_{ij}^{pq} r_{pq}^{k l_f}
    //
    
    // compute contraction -1 * <ii|F12|pp> . <iI|F12|pp>
    contract_tbint_tensor<
      ManyBodyTensors::I_to_T,
      ManyBodyTensors::I_to_T,
      ManyBodyTensors::I_to_mT,
      true,true,false>
      (
        B_gbc2, corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_f12(),
        occ1_act, occ2_act,
        orbs1, orbs2,
        occ1_act, focc2_act,
        orbs2, orbs2,
        dircontract_pp,
        spincase2!=AlphaBeta, tforms_ipip, tforms_ipIp
      );
    
    
    // make particles equivalent, if necessary
    if (!spin_polarized()) {
      symmetrize<false>(B_gbc1,B_gbc1,occ1_act,occ1_act);
      symmetrize<false>(B_gbc2,B_gbc2,occ1_act,occ1_act);
    }
#else  // use compute_X_
    
    compute_X_(B_gbc2,spincase2,occ1_act,occ2_act,
               occ1_act,focc2_act);
    if (occ1_act != occ2_act) {
      compute_X_(B_gbc2,spincase2,occ1_act,occ2_act,
                 focc1_act,occ2_act);
    }
    else {
      B_gbc2.scale(2.0);
      if (spincase2 == AlphaBeta) {
        symmetrize<false>(B_gbc2,B_gbc2,occ1_act,occ2_act);
      }
    }
    
#endif // use compute_X_ ?
    
#endif // include GBC2 ?
    
    if (debug_ > 1) {
      std::string label = prepend_spincase(spincase2,"B(GBC1) contribution");
      B_gbc1.print(label.c_str());
      label = prepend_spincase(spincase2,"B(GBC2) contribution");
      B_gbc2.print(label.c_str());
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

  tim_exit("B(GBC) intermediate");
}

////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
