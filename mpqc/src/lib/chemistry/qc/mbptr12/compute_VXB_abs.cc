//
// compute_VXB_abs.cc
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

#include <chemistry/qc/mbptr12/r12int_eval.h>
#include <chemistry/qc/mbptr12/creator.h>
#include <chemistry/qc/mbptr12/container.h>
#include <chemistry/qc/mbptr12/compute_tbint_tensor.h>
#include <chemistry/qc/mbptr12/contract_tbint_tensor.h>

using namespace std;
using namespace sc;

/**
   R12IntEval::contrib_to_VXB_a_() computes V, X, and B intermediates in standard approximation A using tensor contract functions
*/
void
R12IntEval::contrib_to_VXB_abs_()
{
  if (evaluated_)
    return;

  const bool obs_eq_vbs = r12info()->obs_eq_vbs();
  const bool obs_eq_ribs = r12info()->obs_eq_ribs();
  // commutators only appear in A', A'', and B
  const bool compute_B = (stdapprox() == LinearR12::StdApprox_Ap ||
      stdapprox() == LinearR12::StdApprox_App || stdapprox() == LinearR12::StdApprox_B);

  if (!obs_eq_vbs)
      throw ProgrammingError("R12IntEval::contrib_to_VXB_a_() -- can't use this builder if OBS != VBS",__FILE__,__LINE__);

  Timer tim("mp2-f12a intermeds (tensor contract)");

  // Test new tensor compute function
  for(int s=0; s<nspincases2(); s++) {
      using namespace sc::LinearR12;
      const SpinCase2 spincase2 = static_cast<SpinCase2>(s);
      const SpinCase1 spin1 = case1(spincase2);
      const SpinCase1 spin2 = case2(spincase2);
      Ref<SingleRefInfo> refinfo = r12info()->refinfo();

      if (dim_oo(spincase2).n() == 0)
        continue;

      const Ref<OrbitalSpace>& occ1 = occ(spin1);
      const Ref<OrbitalSpace>& occ2 = occ(spin2);
      const Ref<OrbitalSpace>& occ1_act = occ_act(spin1);
      const Ref<OrbitalSpace>& occ2_act = occ_act(spin2);
      const Ref<OrbitalSpace>& orbs1 = orbs(spin1);
      const Ref<OrbitalSpace>& orbs2 = orbs(spin2);
      const Ref<OrbitalSpace>& xspace1 = xspace(spin1);
      const Ref<OrbitalSpace>& xspace2 = xspace(spin2);

      // for now geminal-generating products must have same equivalence as the occupied orbitals
      const bool occ1_eq_occ2 = (occ1_act == occ2_act);
      const bool x1_eq_x2 = (xspace1 == xspace2);
      if (occ1_eq_occ2 ^ x1_eq_x2) {
      throw ProgrammingError("R12IntEval::contrib_to_VXB_a_() -- this orbital_product cannot be handled yet",__FILE__,__LINE__);
      }

      // are particles 1 and 2 equivalent?
      const bool part1_equiv_part2 =  spincase2 != AlphaBeta || occ1_eq_occ2;
      // Need to antisymmetrize 1 and 2
      const bool antisymmetrize = spincase2 != AlphaBeta;

      // some transforms can be skipped if occ1/occ2 is a subset of x1/x2
      // for now it's always true since can only use ij and pq products to generate geminals
      const bool occ12_in_x12 = true;

      Ref<TwoParticleContraction> tpcontract;
      const ABSMethod absmethod = r12info()->abs_method();
      // "ABS"-type contraction is used for projector 2 ABS/ABS+ method when OBS != RIBS
      // it involves this term +O1O2-V1V2
      if ((absmethod == LinearR12::ABS_ABS ||
       absmethod == LinearR12::ABS_ABSPlus) && !obs_eq_ribs &&
      ansatz()->projector() == LinearR12::Projector_2)
          tpcontract = new ABS_OBS_Contraction(orbs1->rank(),
                                               occ1->rank(),
                                               occ2->rank());
      else
      // involves this term -P1P2
          tpcontract = new CABS_OBS_Contraction(orbs1->rank());

      std::vector<std::string> tforms;
      std::vector<std::string> tforms_f12;
      {
      R12TwoBodyIntKeyCreator tformkey_creator(
          r12info()->moints_runtime4(),
          xspace1,
          orbs1,
          xspace2,
          orbs2,
          r12info()->corrfactor(),true
          );
      fill_container(tformkey_creator,tforms_f12);
      }
      if (!occ12_in_x12) {
      R12TwoBodyIntKeyCreator tformkey_creator(
          r12info()->moints_runtime4(),
          occ1_act,
          orbs1,
          occ2_act,
          orbs2,
          r12info()->corrfactor()
          );
      fill_container(tformkey_creator,tforms);
      }
      else
      tforms.push_back(tforms_f12[0]);

      contract_tbint_tensor<ManyBodyTensors::I_to_T,
      ManyBodyTensors::I_to_T,
      ManyBodyTensors::I_to_T,
      true,false,false>
          (
          V_[s],
          corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_eri(),
          xspace1, xspace2,
          orbs1, orbs2,
          occ1_act, occ2_act,
          orbs1, orbs2,
          tpcontract,
          spincase2!=AlphaBeta, tforms_f12, tforms
      );

      contract_tbint_tensor<ManyBodyTensors::I_to_T,
      ManyBodyTensors::I_to_T,
      ManyBodyTensors::I_to_T,
      true,true,false>
          (
          X_[s], corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_f12(),
          xspace1, xspace2,
          orbs1, orbs2,
          xspace1, xspace2,
          orbs1, orbs2,
          tpcontract,
          spincase2!=AlphaBeta, tforms_f12, tforms_f12
          );
      if (compute_B) {
      contract_tbint_tensor<ManyBodyTensors::I_to_T,
          ManyBodyTensors::I_to_T,
          ManyBodyTensors::I_to_T,
          true,true,false>
          (
          B_[s], corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_t1f12(),
          xspace1, xspace2,
          orbs1, orbs2,
          xspace1, xspace2,
          orbs1, orbs2,
          tpcontract,
          spincase2!=AlphaBeta, tforms_f12, tforms_f12
          );
      contract_tbint_tensor<ManyBodyTensors::I_to_T,
          ManyBodyTensors::I_to_T,
          ManyBodyTensors::I_to_T,
          true,true,false>
          (
          B_[s], corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_t2f12(),
          xspace1, xspace2,
          orbs1, orbs2,
          xspace1, xspace2,
          orbs1, orbs2,
          tpcontract,
          spincase2!=AlphaBeta, tforms_f12, tforms_f12
          );
      B_[s].scale(0.5); RefSCMatrix Bt = B_[s].t(); B_[s].accumulate(Bt);
      }

      if (debug_ >= DefaultPrintThresholds::O4) {
          globally_sum_intermeds_();
          V_[s].print(prepend_spincase(static_cast<SpinCase2>(s),"V(diag+OBS) contribution").c_str());
          X_[s].print(prepend_spincase(static_cast<SpinCase2>(s),"X(diag+OBS) contribution").c_str());
      if (compute_B)
          B_[s].print(prepend_spincase(static_cast<SpinCase2>(s),"B(diag+OBS) contribution").c_str());
      }

#define DEBUG_SKIP_CABS_TERMS 0
#if !DEBUG_SKIP_CABS_TERMS
      // These terms only contribute if Projector=2
      if (!obs_eq_ribs && ansatz()->projector() == LinearR12::Projector_2) {

      const bool cabs_method = (absmethod ==  LinearR12::ABS_CABS ||
                    absmethod == LinearR12::ABS_CABSPlus);
      Ref<OrbitalSpace> rispace1, rispace2;
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

      std::vector<std::string> tforms_imjP;
      std::vector<std::string> tforms_f12_xmyP;
      {
          R12TwoBodyIntKeyCreator tformkey_creator(
          r12info()->moints_runtime4(),
          xspace1,
          occ1,
          xspace2,
          rispace2,
          r12info()->corrfactor(),true
          );
          fill_container(tformkey_creator,tforms_f12_xmyP);
      }
      if (!occ12_in_x12) {
          R12TwoBodyIntKeyCreator tformkey_creator(
          r12info()->moints_runtime4(),
          occ1_act,
          occ1,
          occ2_act,
          rispace2,
          r12info()->corrfactor()
          );
          fill_container(tformkey_creator,tforms_imjP);
      }
      else
          tforms_imjP.push_back(tforms_f12_xmyP[0]);

      contract_tbint_tensor<ManyBodyTensors::I_to_T,
          ManyBodyTensors::I_to_T,
          ManyBodyTensors::I_to_T,
          true,false,false>
          (
          V_[s], corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_eri(),
          xspace1, xspace2,
          occ1, rispace2,
          occ1_act, occ2_act,
          occ1, rispace2,
          tpcontract,
          antisymmetrize, tforms_f12_xmyP, tforms_imjP
          );
      contract_tbint_tensor<ManyBodyTensors::I_to_T,
          ManyBodyTensors::I_to_T,
          ManyBodyTensors::I_to_T,
          true,true,false>
          (
          X_[s], corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_f12(),
          xspace1, xspace2,
          occ1, rispace2,
          xspace1, xspace2,
          occ1, rispace2,
          tpcontract,
          antisymmetrize, tforms_f12_xmyP, tforms_f12_xmyP
          );

      if (compute_B) {
          contract_tbint_tensor<ManyBodyTensors::I_to_T,
          ManyBodyTensors::I_to_T,
          ManyBodyTensors::I_to_T,
          true,true,false>
          (
              B_[s], corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_t1f12(),
              xspace1, xspace2,
              occ1, rispace2,
              xspace1, xspace2,
              occ1, rispace2,
              tpcontract,
              antisymmetrize, tforms_f12_xmyP, tforms_f12_xmyP
          );
          contract_tbint_tensor<ManyBodyTensors::I_to_T,
          ManyBodyTensors::I_to_T,
          ManyBodyTensors::I_to_T,
          true,true,false>
          (
              B_[s], corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_t2f12(),
              xspace1, xspace2,
              occ1, rispace2,
              xspace1, xspace2,
              occ1, rispace2,
              tpcontract,
              antisymmetrize, tforms_f12_xmyP, tforms_f12_xmyP
          );
          B_[s].scale(0.5); RefSCMatrix Bt = B_[s].t(); B_[s].accumulate(Bt);
      }

      // If particles 1 and 2 are not equivalent, also need another set of terms
      if (!part1_equiv_part2) {

          std::vector<std::string> tforms_iPjm;
          std::vector<std::string> tforms_f12_xPym;
          {
          R12TwoBodyIntKeyCreator tformkey_creator(
              r12info()->moints_runtime4(),
              xspace1,
              rispace1,
              xspace2,
              occ2,
              r12info()->corrfactor(),true
              );
          fill_container(tformkey_creator,tforms_f12_xPym);
          }
          if (!occ12_in_x12) {
          R12TwoBodyIntKeyCreator tformkey_creator(
              r12info()->moints_runtime4(),
              occ1_act,
              rispace1,
              occ2_act,
              occ2,
              r12info()->corrfactor()
              );
          fill_container(tformkey_creator,tforms_iPjm);
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
              V_[s], corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_eri(),
              xspace1, xspace2,
              rispace1, occ2,
              occ1_act, occ2_act,
              rispace1, occ2,
              tpcontract,
              antisymmetrize, tforms_f12_xPym, tforms_iPjm
              );
          contract_tbint_tensor<ManyBodyTensors::I_to_T,
          ManyBodyTensors::I_to_T,
          ManyBodyTensors::I_to_T,
          true,true,false>
          (
              X_[s], corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_f12(),
              xspace1, xspace2,
              rispace1, occ2,
              xspace1, xspace2,
              rispace1, occ2,
              tpcontract,
              antisymmetrize, tforms_f12_xPym, tforms_f12_xPym
              );

          if (compute_B) {
          contract_tbint_tensor<ManyBodyTensors::I_to_T,
              ManyBodyTensors::I_to_T,
              ManyBodyTensors::I_to_T,
              true,true,false>
              (
              B_[s], corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_t1f12(),
              xspace1, xspace2,
              rispace1, occ2,
              xspace1, xspace2,
              rispace1, occ2,
              tpcontract,
              antisymmetrize, tforms_f12_xPym, tforms_f12_xPym
              );
          contract_tbint_tensor<ManyBodyTensors::I_to_T,
              ManyBodyTensors::I_to_T,
              ManyBodyTensors::I_to_T,
              true,true,false>
              (
              B_[s], corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_t2f12(),
              xspace1, xspace2,
              rispace1, occ2,
              xspace1, xspace2,
              rispace1, occ2,
              tpcontract,
              antisymmetrize, tforms_f12_xPym, tforms_f12_xPym
              );

          B_[s].scale(0.5); RefSCMatrix Bt = B_[s].t(); B_[s].accumulate(Bt);
          }
      }

      if (!antisymmetrize && part1_equiv_part2) {
          symmetrize<false>(V_[s],V_[s],xspace1,occ1_act);
          symmetrize<false>(X_[s],X_[s],xspace1,xspace1);
          if (compute_B)
          symmetrize<false>(B_[s],B_[s],xspace1,xspace1);
      }

      if (debug_ >= DefaultPrintThresholds::O4) {
          globally_sum_intermeds_();
          V_[s].print(prepend_spincase(static_cast<SpinCase2>(s),"V(diag+OBS+ABS) contribution").c_str());
          X_[s].print(prepend_spincase(static_cast<SpinCase2>(s),"X(diag+OBS+ABS) contribution").c_str());
          if (compute_B)
          B_[s].print(prepend_spincase(static_cast<SpinCase2>(s),"B(diag+OBS+ABS) contribution").c_str());
      }
      }

#endif

  }
}

////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
