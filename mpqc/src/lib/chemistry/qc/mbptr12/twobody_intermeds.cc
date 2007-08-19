//
// twobody_intermeds.cc
//
// Copyright (C) 2007 Edward Valeev
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

#ifdef __GNUG__
#pragma implementation
#endif

#include <math/scmat/local.h>
#include <chemistry/qc/mbptr12/r12int_eval.h>
#include <chemistry/qc/mbptr12/creator.h>
#include <chemistry/qc/mbptr12/container.h>
#include <chemistry/qc/mbptr12/compute_tbint_tensor.h>
#include <chemistry/qc/mbptr12/contract_tbint_tensor.h>

using namespace std;
using namespace sc;

RefSCMatrix
R12IntEval::V(SpinCase2 spincase2,
	      const Ref<MOIndexSpace>& p1,
	      const Ref<MOIndexSpace>& p2)
{
  using namespace sc::LinearR12;

  Ref<LocalSCMatrixKit> local_matrix_kit = new LocalSCMatrixKit();

  const bool obs_eq_vbs = r12info_->basis_vir()->equiv(r12info_->basis());
  const bool obs_eq_ribs = r12info()->basis_ri()->equiv(r12info()->basis());

  const bool p1_eq_p2 = (p1 == p2);
  //if (!p1_eq_p2)
  //  throw FeatureNotImplemented("R12IntEval::V() -- p1 == p2 must be true",__FILE__,__LINE__);

  const SpinCase1 spin1 = case1(spincase2);
  const SpinCase1 spin2 = case2(spincase2);
  Ref<SingleRefInfo> refinfo = r12info()->refinfo();

  RefSCMatrix V;
  if (!spin_polarized() && (spincase2 == AlphaAlpha || spincase2 == BetaBeta)) {
    const unsigned int nx1 = xspace(spin1)->rank();
    const unsigned int np1 = p1->rank();
    V = local_matrix_kit->matrix(new SCDimension(nx1*(nx1-1)/2),
				 new SCDimension(np1*(np1-1)/2));
    RefSCMatrix Vab = this->V(AlphaBeta,p1,p1);
    antisymmetrize(V,Vab,xspace(Alpha),p1);
    return V;
  }


  tim_enter("R12 intermeds (tensor contract): Vpqxy");
  
  const Ref<MOIndexSpace>& xspace1 = xspace(spin1);
  const Ref<MOIndexSpace>& xspace2 = xspace(spin2);
  const Ref<MOIndexSpace>& orbs1 = refinfo->orbs(spin1);
  const Ref<MOIndexSpace>& orbs2 = refinfo->orbs(spin2);

  // The diagonal contribution
  Ref<LinearR12::G12CorrelationFactor> g12ptr; g12ptr << corrfactor();
  Ref<LinearR12::G12NCCorrelationFactor> g12ncptr; g12ncptr << corrfactor();
  Ref<LinearR12::GenG12CorrelationFactor> gg12ptr; gg12ptr << corrfactor();
  Ref<LinearR12::R12CorrelationFactor> r12ptr; r12ptr << corrfactor();
  if (r12ptr.nonnull()) {
    RefSCMatrix I = compute_I_(xspace1,xspace2,p1,p2);
    V = I.clone();
    V.accumulate(I);
  }
  else if (g12ptr.nonnull() || g12ncptr.nonnull() || gg12ptr.nonnull()) {
    std::vector<  Ref<TwoBodyMOIntsTransform> > tforms_f12_xmyn;
    {
      Ref<R12IntEval> thisref(this);
      NewTransformCreator tform_creator(thisref,
					xspace1,p1,
					xspace2,p2,
					true
					);
      fill_container(tform_creator,tforms_f12_xmyn);
    }
    compute_tbint_tensor<ManyBodyTensors::I_to_T,true,false>(
      V, corrfactor()->tbint_type_f12eri(),
      xspace1, p1,
      xspace2, p2,
      antisymmetrize,
      tforms_f12_xmyn);
  }
  
  // are particles 1 and 2 equivalent?
  const bool part1_equiv_part2 =  spincase2 != AlphaBeta || p1_eq_p2;
  // Need to antisymmetrize 1 and 2
  const bool antisymmetrize = spincase2 != AlphaBeta;

  // some transforms can be skipped if p1/p2 is a subset of x1/x2
  const bool p1p2_in_x1x2 = in(*p1,*xspace1) && in(*p2,*xspace2);

  Ref<TwoParticleContraction> tpcontract;
  const ABSMethod absmethod = r12info()->abs_method();
  // "ABS"-type contraction is used for projector 2 ABS/ABS+ method when OBS != RIBS
  // it involves this term +O1O2-V1V2
  if ((absmethod == LinearR12::ABS_ABS ||
       absmethod == LinearR12::ABS_ABSPlus) && !obs_eq_ribs &&
      ansatz()->projector() == LinearR12::Projector_2)
    tpcontract = new ABS_OBS_Contraction(refinfo->orbs(spin1)->rank(),
					 refinfo->occ(spin1)->rank(),
					 refinfo->occ(spin2)->rank());
  else
    // involves this term -P1P2
    tpcontract = new CABS_OBS_Contraction(refinfo->orbs(spin1)->rank());

  std::vector<  Ref<TwoBodyMOIntsTransform> > tforms;
  std::vector<  Ref<TwoBodyMOIntsTransform> > tforms_f12;
  {
    Ref<R12IntEval> thisref(this);
    NewTransformCreator tform_creator(
				      thisref,
				      xspace1,
				      orbs1,
				      xspace2,
				      orbs2,true
				      );
    fill_container(tform_creator,tforms_f12);
  }
  if (!p1p2_in_x1x2) {
    Ref<R12IntEval> thisref(this);
    NewTransformCreator tform_creator(
				      thisref,
				      p1,
				      orbs1,
				      p2,
				      orbs2
				      );
    fill_container(tform_creator,tforms);
  }
  else
    tforms.push_back(tforms_f12[0]);

  contract_tbint_tensor<ManyBodyTensors::I_to_T,
    ManyBodyTensors::I_to_T,
    ManyBodyTensors::I_to_T,
    true,false,false>
    (
     V, corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_eri(),
     xspace1, xspace2,
     orbs1, orbs2,
     p1, p2,
     orbs1, orbs2,
     tpcontract,
     spincase2!=AlphaBeta, tforms_f12, tforms
     );

  if (debug_ >= DefaultPrintThresholds::O4) {
    V.print(prepend_spincase(spincase2,"Vpqxy: diag+OBS contribution").c_str());
  }

  // These terms only contribute if Projector=2
  if (!obs_eq_ribs && ansatz()->projector() == LinearR12::Projector_2) {
    
    const bool cabs_method = (absmethod ==  LinearR12::ABS_CABS ||
			      absmethod == LinearR12::ABS_CABSPlus);
    const Ref<MOIndexSpace>& occ1 = occ(spin1);
    const Ref<MOIndexSpace>& occ2 = occ(spin2);
    Ref<MOIndexSpace> rispace1, rispace2;
    if (cabs_method) {
      rispace1 = r12info()->ribs_space(spin1);
      rispace2 = r12info()->ribs_space(spin2);
    }
    else {
      rispace1 = r12info()->ribs_space();
      rispace2 = r12info()->ribs_space();
    }
    // If particles are equivalent, <ij|Pm> = <ji|mP>, hence in the same set of integrals.
    // Can then skip <ij|Pm>, simply add 2<ij|mP> and (anti)symmetrize
    Ref<TwoParticleContraction> tpcontract =
      new Direct_Contraction(
			     occ1->rank(),
			     rispace2->rank(),part1_equiv_part2 ? -2.0 : -1.0
			     );

    std::vector<  Ref<TwoBodyMOIntsTransform> > tforms_imjP;
    std::vector<  Ref<TwoBodyMOIntsTransform> > tforms_f12_xmyP;
    {
      Ref<R12IntEval> thisref(this);
      NewTransformCreator tform_creator(
					thisref,
					xspace1,
					occ1,
					xspace2,
					rispace2,true
					);
      fill_container(tform_creator,tforms_f12_xmyP);
    }
    if (!p1p2_in_x1x2) {
      Ref<R12IntEval> thisref(this);
      NewTransformCreator tform_creator(
					thisref,
					p1,
					occ1,
					p2,
					rispace2
					);
      fill_container(tform_creator,tforms_imjP);
    }
    else
      tforms_imjP.push_back(tforms_f12_xmyP[0]);
    
    contract_tbint_tensor<ManyBodyTensors::I_to_T,
      ManyBodyTensors::I_to_T,
      ManyBodyTensors::I_to_T,
      true,false,false>
      (
       V, corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_eri(),
       xspace1, xspace2,
       occ1, rispace2,
       p1, p2,
       occ1, rispace2,
       tpcontract,
       antisymmetrize, tforms_f12_xmyP, tforms_imjP
       );

    // If particles 1 and 2 are not equivalent, also need another set of terms
    if (!part1_equiv_part2) {
      
      std::vector<  Ref<TwoBodyMOIntsTransform> > tforms_iPjm;
      std::vector<  Ref<TwoBodyMOIntsTransform> > tforms_f12_xPym;
      {
	Ref<R12IntEval> thisref(this);
	NewTransformCreator tform_creator(
					  thisref,
					  xspace1,
					  rispace1,
					  xspace2,
					  occ2,true
					  );
	fill_container(tform_creator,tforms_f12_xPym);
      }
      if (!p1p2_in_x1x2) {
	Ref<R12IntEval> thisref(this);
	NewTransformCreator tform_creator(
					  thisref,
					  p1,
					  rispace1,
					  p2,
					  occ2
					  );
	fill_container(tform_creator,tforms_iPjm);
      }
      else
	tforms_iPjm.push_back(tforms_f12_xPym[0]);
      
      Ref<TwoParticleContraction> tpcontract =
	new Direct_Contraction(
			       rispace1->rank(),
			       occ2->rank(),-1.0
			       );
      
      contract_tbint_tensor<ManyBodyTensors::I_to_T,
	ManyBodyTensors::I_to_T,
	ManyBodyTensors::I_to_T,
	true,false,false>
	(
	 V, corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_eri(),
	 xspace1, xspace2,
	 rispace1, occ2,
	 p1, p2,
	 rispace1, occ2,
	 tpcontract,
	 antisymmetrize, tforms_f12_xPym, tforms_iPjm
	 );
    } // if part1_equiv_part2
  } // ABS != OBS

  if (!antisymmetrize && part1_equiv_part2) {
    symmetrize<false>(V,V,xspace1,p1);
  }

  if (debug_ >= DefaultPrintThresholds::O4) {
    V.print(prepend_spincase(spincase2,"Vpqxy: diag+OBS+ABS contribution").c_str());
  }
  
  tim_exit("R12 intermeds (tensor contract): Vpqxy");
}
