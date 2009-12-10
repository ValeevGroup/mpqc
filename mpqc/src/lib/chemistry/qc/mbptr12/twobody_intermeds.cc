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
              const Ref<OrbitalSpace>& p1,
              const Ref<OrbitalSpace>& p2)
{
  const LinearR12::ABSMethod absmethod = r12world()->r12tech()->abs_method();
  const bool obs_eq_ribs = r12world()->obs_eq_ribs();
  // "ABS"-type contraction is used for projector 2 ABS/ABS+ method when OBS != RIBS
  // it involves this term +O1O2-V1V2
  if ((absmethod == LinearR12::ABS_ABS ||
       absmethod == LinearR12::ABS_ABSPlus) && !obs_eq_ribs &&
      ansatz()->projector() == LinearR12::Projector_2)
    return this->V_abs(spincase2, p1, p2);
  else
    return this->V_cabs(spincase2, p1, p2);
}

RefSCMatrix
R12IntEval::V_abs(SpinCase2 spincase2,
                  const Ref<OrbitalSpace>& p1,
                  const Ref<OrbitalSpace>& p2)
{
  using namespace sc::LinearR12;

  Ref<LocalSCMatrixKit> local_matrix_kit = new LocalSCMatrixKit();

  const bool obs_eq_vbs = r12world()->obs_eq_vbs();
  const bool obs_eq_ribs = r12world()->obs_eq_ribs();

  const bool p1_eq_p2 = (p1 == p2);
  // are particles 1 and 2 equivalent?
  const bool part1_equiv_part2 =  spincase2 != AlphaBeta || p1_eq_p2;
  // Need to antisymmetrize 1 and 2
  const bool antisymmetrize = spincase2 != AlphaBeta;
  //if (!p1_eq_p2 && antisymmetrize)
  //  throw FeatureNotImplemented("R12IntEval::V() -- p1 == p2 must be true if AA or BB spin case",__FILE__,__LINE__);

  const SpinCase1 spin1 = case1(spincase2);
  const SpinCase1 spin2 = case2(spincase2);

  RefSCMatrix V;
  if (!spin_polarized() && (spincase2 == AlphaAlpha || spincase2 == BetaBeta)) {
    const unsigned int nx1 = xspace(spin1)->rank();
    const unsigned int np1 = p1->rank();
    V = local_matrix_kit->matrix(dim_f12(spincase2),
				 new SCDimension(np1*(np1-1)/2));
    RefSCMatrix Vab = this->V(AlphaBeta,p1,p1);
    sc::antisymmetrize(V,Vab,xspace(Alpha),p1);
    return V;
  }

  Timer tim("R12 intermeds (tensor contract): Vpqxy");

  const Ref<OrbitalSpace>& xspace1 = xspace(spin1);
  const Ref<OrbitalSpace>& xspace2 = xspace(spin2);
  const Ref<OrbitalSpace>& orbs1 = orbs(spin1);
  const Ref<OrbitalSpace>& orbs2 = orbs(spin2);

  // some transforms can be skipped if p1/p2 is a subset of x1/x2
  const bool p1p2_in_x1x2 = in(*p1,*xspace1) && in(*p2,*xspace2);

  // allocate V
  const unsigned int np12 = p1_eq_p2 && spincase2 != AlphaBeta ? p1->rank()*(p1->rank()-1)/2 : p1->rank()*p2->rank();
  RefSCDimension dim_p12 = new SCDimension(np12);
  V = local_matrix_kit->matrix(dim_f12(spincase2), dim_p12);
  V.assign(0.0);

  // The diagonal contribution
  Ref<LinearR12::G12CorrelationFactor> g12ptr; g12ptr << corrfactor();
  Ref<LinearR12::G12NCCorrelationFactor> g12ncptr; g12ncptr << corrfactor();
  Ref<LinearR12::GenG12CorrelationFactor> gg12ptr; gg12ptr << corrfactor();
  Ref<LinearR12::R12CorrelationFactor> r12ptr; r12ptr << corrfactor();
  if (r12ptr.nonnull()) {
    RefSCMatrix I = compute_I_(xspace1,xspace2,p1,p2);
    if (!antisymmetrize)
      V.accumulate(I);
    else
      sc::antisymmetrize<true>(V,I,xspace1,xspace2,p1,p2);
  }
  else if (g12ptr.nonnull() || g12ncptr.nonnull() || gg12ptr.nonnull()) {
    std::vector<std::string> tforms_f12_xmyn;
    {
      R12TwoBodyIntKeyCreator tformkey_creator(moints_runtime4(),
					xspace1,p1,
					xspace2,p2,
					corrfactor(),true
					);
      fill_container(tformkey_creator,tforms_f12_xmyn);
    }
    compute_tbint_tensor<ManyBodyTensors::I_to_T,true,false>(
      V, corrfactor()->tbint_type_f12eri(),
      xspace1, p1,
      xspace2, p2,
      antisymmetrize,
      tforms_f12_xmyn);
  }
  if (debug_ >= DefaultPrintThresholds::O4) {
    V.print(prepend_spincase(spincase2,"Vpqxy: diag contribution").c_str());
  }

  Ref<TwoParticleContraction> tpcontract = new ABS_OBS_Contraction(orbs1->rank(),
                                                                   occ(spin1)->rank(),
                                                                   occ(spin2)->rank());

  std::vector<std::string> tforms;
  std::vector<std::string> tforms_f12;
  {
    R12TwoBodyIntKeyCreator tformkey_creator(
				      moints_runtime4(),
				      xspace1,
				      orbs1,
				      xspace2,
				      orbs2,
				      corrfactor(),true
				      );
    fill_container(tformkey_creator,tforms_f12);
  }
  if (!p1p2_in_x1x2) {
    const std::string tform_key = ParsedTwoBodyFourCenterIntKey::key(p1->id(),p2->id(),
                                                                     orbs1->id(),orbs2->id(),
                                                                     std::string("ERI"),
                                                                     std::string(TwoBodyIntLayout::b1b2_k1k2));
    tforms.push_back(tform_key);
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

    const LinearR12::ABSMethod absmethod = r12world()->r12tech()->abs_method();
    const bool cabs_method = (absmethod ==  LinearR12::ABS_CABS ||
			      absmethod == LinearR12::ABS_CABSPlus);
    const Ref<OrbitalSpace>& occ1 = occ(spin1);
    const Ref<OrbitalSpace>& occ2 = occ(spin2);
    Ref<OrbitalSpace> rispace1, rispace2;
    if (cabs_method) {
      rispace1 = r12world()->cabs_space(spin1);
      rispace2 = r12world()->cabs_space(spin2);
    }
    else {
      rispace1 = r12world()->ribs_space();
      rispace2 = r12world()->ribs_space();
    }
    // If particles are equivalent, <ij|Pm> = <ji|mP>, hence in the same set of integrals.
    // Can then skip <ij|Pm>, simply add 2<ij|mP> and (anti)symmetrize
    Ref<TwoParticleContraction> tpcontract =
      new Direct_Contraction(
			     occ1->rank(),
			     rispace2->rank(),part1_equiv_part2 ? -2.0 : -1.0
			     );

    std::vector<std::string> tforms_imjP;
    std::vector<std::string> tforms_f12_xmyP;
    {
      R12TwoBodyIntKeyCreator tformkey_creator(
					moints_runtime4(),
					xspace1,
					occ1,
					xspace2,
					rispace2,
                    corrfactor(),true
					);
      fill_container(tformkey_creator,tforms_f12_xmyP);
    }
    if (!p1p2_in_x1x2) {
      const std::string tform_key = ParsedTwoBodyFourCenterIntKey::key(p1->id(),p2->id(),
                                                             occ1->id(),rispace2->id(),
                                                             std::string("ERI"),
                                                             std::string(TwoBodyIntLayout::b1b2_k1k2));
      tforms_imjP.push_back(tform_key);
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

      std::vector<std::string> tforms_iPjm;
      std::vector<std::string> tforms_f12_xPym;
      {
	R12TwoBodyIntKeyCreator tformkey_creator(
					  moints_runtime4(),
					  xspace1,
					  rispace1,
					  xspace2,
					  occ2,
                      corrfactor(),true
					  );
	fill_container(tformkey_creator,tforms_f12_xPym);
      }
      if (!p1p2_in_x1x2) {
        const std::string tform_key = ParsedTwoBodyFourCenterIntKey::key(p1->id(),p2->id(),
                                                               rispace1->id(),occ2->id(),
                                                               std::string("ERI"),
                                                               std::string(TwoBodyIntLayout::b1b2_k1k2));
        tforms_iPjm.push_back(tform_key);
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

  return V;
}


RefSCMatrix
R12IntEval::V_cabs(SpinCase2 spincase2,
                   const Ref<OrbitalSpace>& p1,
                   const Ref<OrbitalSpace>& p2)
{
  using namespace sc::LinearR12;

  Ref<LocalSCMatrixKit> local_matrix_kit = new LocalSCMatrixKit();

  const bool obs_eq_vbs = r12world()->obs_eq_vbs();
  const bool obs_eq_ribs = r12world()->obs_eq_ribs();

  const bool p1_eq_p2 = (p1 == p2);
  // are particles 1 and 2 equivalent?
  const bool part1_equiv_part2 =  spincase2 != AlphaBeta || p1_eq_p2;
  // Need to antisymmetrize 1 and 2
  const bool antisymmetrize = spincase2 != AlphaBeta;
  //if (!p1_eq_p2 && antisymmetrize)
  //  throw FeatureNotImplemented("R12IntEval::V() -- p1 == p2 must be true if AA or BB spin case",__FILE__,__LINE__);

  const SpinCase1 spin1 = case1(spincase2);
  const SpinCase1 spin2 = case2(spincase2);

  RefSCMatrix V;
  if (!spin_polarized() && (spincase2 == AlphaAlpha || spincase2 == BetaBeta)) {
    const unsigned int nx1 = xspace(spin1)->rank();
    const unsigned int np1 = p1->rank();
    V = local_matrix_kit->matrix(dim_f12(spincase2),
                 new SCDimension(np1*(np1-1)/2));
    RefSCMatrix Vab = this->V(AlphaBeta,p1,p1);
    sc::antisymmetrize(V,Vab,xspace(Alpha),p1);
    return V;
  }

  Timer tim("R12 intermeds (tensor contract): Vpqxy");

  const Ref<OrbitalSpace>& xspace1 = xspace(spin1);
  const Ref<OrbitalSpace>& xspace2 = xspace(spin2);
  const Ref<OrbitalSpace>& orbs1 = orbs(spin1);
  const Ref<OrbitalSpace>& orbs2 = orbs(spin2);

  // some transforms can be skipped if p1/p2 equals x1/x2
  const bool p1p2_eq_x1x2 = (p1 == xspace1) && (p2 == xspace2);

  // allocate V
  const unsigned int np12 = p1_eq_p2 && spincase2 != AlphaBeta ? p1->rank()*(p1->rank()-1)/2 : p1->rank()*p2->rank();
  RefSCDimension dim_p12 = new SCDimension(np12);
  V = local_matrix_kit->matrix(dim_f12(spincase2), dim_p12);
  V.assign(0.0);

  // The diagonal contribution
  Ref<LinearR12::G12CorrelationFactor> g12ptr; g12ptr << corrfactor();
  Ref<LinearR12::G12NCCorrelationFactor> g12ncptr; g12ncptr << corrfactor();
  Ref<LinearR12::GenG12CorrelationFactor> gg12ptr; gg12ptr << corrfactor();
  Ref<LinearR12::R12CorrelationFactor> r12ptr; r12ptr << corrfactor();
  if (r12ptr.nonnull()) {
    RefSCMatrix I = compute_I_(xspace1,xspace2,p1,p2);
    if (!antisymmetrize)
      V.accumulate(I);
    else
      sc::antisymmetrize<true>(V,I,xspace1,xspace2,p1,p2);
  }
  else if (g12ptr.nonnull() || g12ncptr.nonnull() || gg12ptr.nonnull()) {
    std::vector<std::string> tforms_f12_xmyn;
    {
      R12TwoBodyIntKeyCreator tformkey_creator(moints_runtime4(),
                    xspace1,p1,
                    xspace2,p2,
                    corrfactor(),true
                    );
      fill_container(tformkey_creator,tforms_f12_xmyn);
    }
    compute_tbint_tensor<ManyBodyTensors::I_to_T,true,false>(
      V, corrfactor()->tbint_type_f12eri(),
      xspace1, p1,
      xspace2, p2,
      antisymmetrize,
      tforms_f12_xmyn);
  }
  if (debug_ >= DefaultPrintThresholds::O4) {
    V.print(prepend_spincase(spincase2,"Vpqxy: diag contribution").c_str());
  }

  std::vector<std::string> tforms;
  std::vector<std::string> tforms_f12;
  {
    R12TwoBodyIntKeyCreator tformkey_creator(
                      moints_runtime4(),
                      xspace1,
                      orbs1,
                      xspace2,
                      orbs2,
                      corrfactor(),true
                      );
    fill_container(tformkey_creator,tforms_f12);
  }
  if (!p1p2_eq_x1x2) {
    const std::string tform_key = ParsedTwoBodyFourCenterIntKey::key(p1->id(),p2->id(),
                                                                     orbs1->id(),orbs2->id(),
                                                                     std::string("ERI"),
                                                                     std::string(TwoBodyIntLayout::b1b2_k1k2));
    tforms.push_back(tform_key);
  }
  else
    tforms.push_back(tforms_f12[0]);

  contract_tbint_tensor<true,false>
    (
     V, corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_eri(),
     -1.0,
     xspace1, xspace2,
     orbs1, orbs2,
     p1, p2,
     orbs1, orbs2,
     spincase2!=AlphaBeta, tforms_f12, tforms
     );

  if (debug_ >= DefaultPrintThresholds::O4) {
    V.print(prepend_spincase(spincase2,"Vpqxy: diag+OBS contribution").c_str());
  }

  // These terms only contribute if Projector=2
  if (!obs_eq_ribs && ansatz()->projector() == LinearR12::Projector_2) {

    const Ref<OrbitalSpace>& occ1 = occ(spin1);
    const Ref<OrbitalSpace>& occ2 = occ(spin2);
    Ref<OrbitalSpace> rispace1, rispace2;
    rispace1 = r12world()->cabs_space(spin1);
    rispace2 = r12world()->cabs_space(spin2);
    // If particles are equivalent, <ij|Pm> = <ji|mP>, hence in the same set of integrals.
    // Can then skip <ij|Pm>, simply add 2<ij|mP> and (anti)symmetrize
    const double perm_factor = part1_equiv_part2 ? -2.0 : -1.0;

    std::vector<std::string> tforms_imjP;
    std::vector<std::string> tforms_f12_xmyP;
    {
      R12TwoBodyIntKeyCreator tformkey_creator(
                    moints_runtime4(),
                    xspace1,
                    occ1,
                    xspace2,
                    rispace2,
                    corrfactor(),true
                    );
      fill_container(tformkey_creator,tforms_f12_xmyP);
    }
    if (!p1p2_eq_x1x2) {
      const std::string tform_key = ParsedTwoBodyFourCenterIntKey::key(p1->id(),p2->id(),
                                                             occ1->id(),rispace2->id(),
                                                             std::string("ERI"),
                                                             std::string(TwoBodyIntLayout::b1b2_k1k2));
      tforms_imjP.push_back(tform_key);
    }
    else
      tforms_imjP.push_back(tforms_f12_xmyP[0]);

    contract_tbint_tensor<true,false>
      (
       V, corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_eri(),
       perm_factor,
       xspace1, xspace2,
       occ1, rispace2,
       p1, p2,
       occ1, rispace2,
       antisymmetrize, tforms_f12_xmyP, tforms_imjP
       );

    // If particles 1 and 2 are not equivalent, also need another set of terms
    if (!part1_equiv_part2) {

      std::vector<std::string> tforms_iPjm;
      std::vector<std::string> tforms_f12_xPym;
      {
    R12TwoBodyIntKeyCreator tformkey_creator(
                      moints_runtime4(),
                      xspace1,
                      rispace1,
                      xspace2,
                      occ2,
                      corrfactor(),true
                      );
    fill_container(tformkey_creator,tforms_f12_xPym);
      }
      if (!p1p2_eq_x1x2) {
        const std::string tform_key = ParsedTwoBodyFourCenterIntKey::key(p1->id(),p2->id(),
                                                               rispace1->id(),occ2->id(),
                                                               std::string("ERI"),
                                                               std::string(TwoBodyIntLayout::b1b2_k1k2));
        tforms_iPjm.push_back(tform_key);
      }
      else
    tforms_iPjm.push_back(tforms_f12_xPym[0]);

      contract_tbint_tensor<true,false>
    (
     V, corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_eri(),
     -1.0,
     xspace1, xspace2,
     rispace1, occ2,
     p1, p2,
     rispace1, occ2,
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

  return V;
}

RefSymmSCMatrix
R12IntEval::P(SpinCase2 spincase2)
{
  using namespace sc::LinearR12;

  Ref<LocalSCMatrixKit> local_matrix_kit = new LocalSCMatrixKit();

  const bool obs_eq_vbs = r12world()->obs_eq_vbs();
  const bool obs_eq_ribs = r12world()->obs_eq_ribs();

  const SpinCase1 spin1 = case1(spincase2);
  const SpinCase1 spin2 = case2(spincase2);
  const Ref<OrbitalSpace>& xspace1 = xspace(spin1);
  const Ref<OrbitalSpace>& xspace2 = xspace(spin2);
  const Ref<OrbitalSpace>& orbs1 = orbs(spin1);
  const Ref<OrbitalSpace>& orbs2 = orbs(spin2);
  const Ref<OrbitalSpace>& occ1 = occ(spin1);
  const Ref<OrbitalSpace>& occ2 = occ(spin2);
  const LinearR12::ABSMethod absmethod = r12world()->r12tech()->abs_method();
  if (absmethod == LinearR12::ABS_ABS ||
      absmethod == LinearR12::ABS_ABSPlus)
    throw InputError("R12IntEval::P() -- cannot use if abs_method = abs/abs+. Try cabs or cabs+.");
  const Ref<OrbitalSpace>& cabs1 = r12world()->cabs_space(spin1);
  const Ref<OrbitalSpace>& cabs2 = r12world()->cabs_space(spin2);
  const RefSCDimension dimf12 = dim_f12(spincase2);

  // are particles 1 and 2 equivalent?
  const bool part1_equiv_part2 =  spincase2 != AlphaBeta || (xspace1 == xspace2);
  // Need to antisymmetrize 1 and 2
  const bool antisymmetrize = spincase2 != AlphaBeta;
  // For RHF compute alpha-alpha and beta-beta P from alpha-beta P
  if (!spin_polarized() && (spincase2 == AlphaAlpha || spincase2 == BetaBeta)) {
    RefSymmSCMatrix Paa = local_matrix_kit->symmmatrix(dimf12);
    RefSymmSCMatrix Pab = this->P(AlphaBeta);
    sc::antisymmetrize<false>(Paa,Pab,xspace1);
    return Paa;
  }

  Timer tim("R12 intermeds (tensor contract): Pxyow");

  // allocate P
  RefSCMatrix P = local_matrix_kit->matrix(dimf12,dimf12);
  P.assign(0.0);

  //
  // Several contributions depend on the form of the correlation factor:
  // 1) the diagonal contribution = RGR, i.e. f(r12) * f(r12) / r12
  // 2) RG = f(r12) / r12
  //
  Ref<LinearR12::G12CorrelationFactor> g12ptr; g12ptr << corrfactor();
  Ref<LinearR12::G12NCCorrelationFactor> g12ncptr; g12ncptr << corrfactor();
  Ref<LinearR12::GenG12CorrelationFactor> gg12ptr; gg12ptr << corrfactor();
  Ref<LinearR12::R12CorrelationFactor> r12ptr; r12ptr << corrfactor();

  //
  // Diagonal contribution: P_{xy}^{wz} = (RGR)_{xy}^{wz}
  //
  if (r12ptr.nonnull()) {
    std::vector<std::string> tforms_f12_xoyw;
    {
      R12TwoBodyIntKeyCreator tformkey_creator(moints_runtime4(),
                    xspace1,orbs1,
                    xspace2,orbs2,
                    corrfactor(),
                    true
                    );
      fill_container(tformkey_creator,tforms_f12_xoyw);
    }
    compute_tbint_tensor<ManyBodyTensors::I_to_T,true,false>(
      P, corrfactor()->tbint_type_f12(),
      xspace1, xspace1,
      xspace2, xspace2,
      antisymmetrize,
      tforms_f12_xoyw);
  }
  else if (g12ptr.nonnull() || g12ncptr.nonnull() || gg12ptr.nonnull()) {
    std::vector<std::string> tforms_f12f12_xoyw;
    {
      R12TwoBodyIntKeyCreator tformkey_creator(
        moints_runtime4(),
        xspace1,
        xspace1,
        xspace2,
        xspace2,
        corrfactor(),true,true
        );
      fill_container(tformkey_creator,tforms_f12f12_xoyw);
    }
    compute_tbint_tensor<ManyBodyTensors::I_to_T,true,true>(
      P, corrfactor()->tbint_type_f12eri(),
      xspace1, xspace1,
      xspace2, xspace2,
      antisymmetrize,
      tforms_f12f12_xoyw);
  }
  if (debug_ >= DefaultPrintThresholds::O4) {
    P.print(prepend_spincase(spincase2,"Pxyow: diag contribution").c_str());
  }

  //
  // OBS contribution: P_{xy}^{wz} -= 1/2 ( V_{xy}^{pq} + (RG)_{xy}^{pq} ) r_{pq}^{wz}
  //
  RefSCMatrix V_pp = V(spincase2,orbs1,orbs2);
  V_pp.print(prepend_spincase(spincase2,"Pxyow: V_pp").c_str());
  // accumulate RG into V
  if (r12ptr.nonnull()) {
    RefSCMatrix I = compute_I_(xspace1,xspace2,orbs1,orbs2);
    if (!antisymmetrize)
      V_pp.accumulate(I);
    else
      sc::antisymmetrize<true>(V_pp,I,xspace1,xspace2,orbs1,orbs2);
  }
  else if (g12ptr.nonnull() || g12ncptr.nonnull() || gg12ptr.nonnull()) {
    std::vector<std::string> tforms_f12_xpyq;
    {
      R12TwoBodyIntKeyCreator tformkey_creator(
        moints_runtime4(),
        xspace1,
        orbs1,
        xspace2,
        orbs2,
        corrfactor(),true,false
        );
      fill_container(tformkey_creator,tforms_f12_xpyq);
    }
    compute_tbint_tensor<ManyBodyTensors::I_to_T,true,false>(
      V_pp, corrfactor()->tbint_type_f12eri(),
      xspace1, orbs1,
      xspace2, orbs2,
      antisymmetrize,
      tforms_f12_xpyq);
  }
  V_pp.print(prepend_spincase(spincase2,"Pxyow: V_pp + RG_pp").c_str());
  // get R_pp
  RefSCMatrix R_pp = local_matrix_kit->matrix(dimf12,dim_aa(spincase2));  R_pp.assign(0.0);
  {
    std::vector<std::string> tforms_f12_xpyq;
    {
      R12TwoBodyIntKeyCreator tformkey_creator(
        moints_runtime4(),
        xspace1,
        orbs1,
        xspace2,
        orbs2,
        corrfactor(),true,false
        );
      fill_container(tformkey_creator,tforms_f12_xpyq);
    }
    compute_tbint_tensor<ManyBodyTensors::I_to_T,true,false>(
      R_pp, corrfactor()->tbint_type_f12(),
      xspace1, orbs1,
      xspace2, orbs2,
      antisymmetrize,
      tforms_f12_xpyq);
  }
  V_pp.scale(-1.0);
  P.accumulate(V_pp * R_pp.t());
  V_pp = 0;  R_pp = 0;

  if (debug_ >= DefaultPrintThresholds::O4) {
    P.print(prepend_spincase(spincase2,"Pxyow: diag+OBS contribution").c_str());
  }

  //
  // ABS contribution: P_{xy}^{wz} -= (V_{xy}^{ma'} + (RG)_{xy}^{ma'} ) r_{ma'}^{wz}
  //
  RefSCMatrix V_iA = V(spincase2,occ1,cabs2);
  V_iA.print(prepend_spincase(spincase2,"Pxyow: V_iA").c_str());
  // accumulate RG into V
  if (r12ptr.nonnull()) {
    RefSCMatrix I = compute_I_(xspace1,xspace2,occ1,cabs2);
    if (!antisymmetrize)
      V_iA.accumulate(I);
    else
      sc::antisymmetrize<true>(V_iA,I,xspace1,xspace2,occ1,cabs2);
  }
  else if (g12ptr.nonnull() || g12ncptr.nonnull() || gg12ptr.nonnull()) {
    std::vector<std::string> tforms_f12_xiyA;
    {
      R12TwoBodyIntKeyCreator tformkey_creator(
        moints_runtime4(),
        xspace1,
        occ1,
        xspace2,
        cabs2,
        corrfactor(),true,false
        );
      fill_container(tformkey_creator,tforms_f12_xiyA);
    }
    compute_tbint_tensor<ManyBodyTensors::I_to_T,true,false>(
      V_iA, corrfactor()->tbint_type_f12eri(),
      xspace1, occ1,
      xspace2, cabs2,
      antisymmetrize,
      tforms_f12_xiyA);
  }
  V_iA.print(prepend_spincase(spincase2,"Pxyow: V_iA + RG_iA").c_str());
  // get R_iA
  const unsigned int niA = occ1->rank() * cabs2->rank();
  RefSCMatrix R_iA = local_matrix_kit->matrix(dimf12,new SCDimension(niA));  R_iA.assign(0.0);
  {
    std::vector<std::string> tforms_f12_xiyA;
    {
      R12TwoBodyIntKeyCreator tformkey_creator(
        moints_runtime4(),
        xspace1,
        occ1,
        xspace2,
        cabs2,
        corrfactor(),true,false
        );
      fill_container(tformkey_creator,tforms_f12_xiyA);
    }
    compute_tbint_tensor<ManyBodyTensors::I_to_T,true,false>(
      R_iA, corrfactor()->tbint_type_f12(),
      xspace1, occ1,
      xspace2, cabs2,
      antisymmetrize,
      tforms_f12_xiyA);
  }
  // if particles 1 and 2 are equivalent, iA and Ai contributions are identical, hence just scale iA contribution by 2
  V_iA.scale(part1_equiv_part2 ? -2.0 : -1.0);
  P.accumulate(V_iA * R_iA.t());
  V_iA = 0;  R_iA = 0;
  if (!part1_equiv_part2) {
    RefSCMatrix V_Ai = V(spincase2,cabs1,occ2);
    V_Ai.print(prepend_spincase(spincase2,"Pxyow: V_Ai").c_str());
    // accumulate RG into V
    if (r12ptr.nonnull()) {
      RefSCMatrix I = compute_I_(xspace1,xspace2,cabs1,occ2);
      if (!antisymmetrize)
        V_Ai.accumulate(I);
      else
        sc::antisymmetrize<true>(V_Ai,I,xspace1,xspace2,cabs1,occ2);
    }
    else if (g12ptr.nonnull() || g12ncptr.nonnull() || gg12ptr.nonnull()) {
      std::vector<std::string> tforms_f12_xAyi;
      {
        R12TwoBodyIntKeyCreator tformkey_creator(
          moints_runtime4(),
          xspace1,
          cabs1,
          xspace2,
          occ2,
          corrfactor(),true,false
          );
        fill_container(tformkey_creator,tforms_f12_xAyi);
      }
      compute_tbint_tensor<ManyBodyTensors::I_to_T,true,false>(
        V_Ai, corrfactor()->tbint_type_f12eri(),
        xspace1, cabs1,
        xspace2, occ2,
        antisymmetrize,
        tforms_f12_xAyi);
    }
    V_Ai.print(prepend_spincase(spincase2,"Pxyow: V_Ai + RG_Ai").c_str());
    // get R_Ai
    const unsigned int nAi = cabs1->rank() * occ2->rank();
    RefSCMatrix R_Ai = local_matrix_kit->matrix(dimf12,new SCDimension(nAi));  R_Ai.assign(0.0);
    {
      std::vector<std::string> tforms_f12_xAyi;
      {
        R12TwoBodyIntKeyCreator tformkey_creator(
          moints_runtime4(),
          xspace1,
          cabs1,
          xspace2,
          occ2,
          corrfactor(),true,false
          );
        fill_container(tformkey_creator,tforms_f12_xAyi);
      }
      compute_tbint_tensor<ManyBodyTensors::I_to_T,true,false>(
        R_Ai, corrfactor()->tbint_type_f12(),
        xspace1, cabs1,
        xspace2, occ2,
        antisymmetrize,
        tforms_f12_xAyi);
    }
    V_Ai.scale(-1.0);
    P.accumulate(V_Ai * R_Ai.t());
    V_Ai = 0;  R_Ai = 0;
  }

  // Because we skipped Ai term if particles 1 and 2 are equivalent, must make particles equivalent explicitly
  if (part1_equiv_part2)
    symmetrize<false>(P,P,xspace1,xspace1);

  if (debug_ >= DefaultPrintThresholds::O4) {
    P.print(prepend_spincase(spincase2,"Pxyow: diag+OBS+ABS contribution").c_str());
  }

  return to_lower_triangle(P);
}
