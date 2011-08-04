//
// compute_VXB_c_ansatz1.cc
//
// Copyright (C) 2008 Martin Torheyden
//
// Author: Martin Torheyden <mtorhey@vt.edu>
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
   R12IntEval::contrib_to_VXB_c_ansatz1_() computes V, X, and B intermediates in standard approximation C, Ansatz 1 using tensor contract functions
*/

void R12IntEval::contrib_to_VXB_c_ansatz1_() {
  if (evaluated_)
    return;

  Ref<R12IntEval> thisref(this);

  Timer timer("mp2-f12c_ansatz1 intermeds (tensor contract)");

  //if((r12world()->ansatz()->orbital_product_gg() != R12Technology::OrbProdgg_pq) ||
  //    (r12world()->ansatz()->orbital_product_GG() != R12Technology::OrbProdGG_pq)) {
  //  throw InputError("R12IntEval::contrib_to_VXB_c_ansatz1() -- gg space and GG space must be both spanned by a pq orbital product Ansatz.",__FILE__,__LINE__);
  //}

  if(r12world()->r12tech()->ansatz()->projector() != R12Technology::Projector_1) {
    throw InputError("R12IntEval::contrib_to_VXB_c_ansatz1() -- this routine only works with an Ansatz 1 projector.",__FILE__,__LINE__);
  }

  for(int s=0; s<nspincases2(); s++) {
    using namespace sc::mbptr12;
    const SpinCase2 spincase2 = static_cast<SpinCase2>(s);
    const SpinCase1 spin1 = case1(spincase2);
    const SpinCase1 spin2 = case2(spincase2);

    if (dim_oo(spincase2).n() == 0)
      continue;

    const Ref<OrbitalSpace>& gg1_space = ggspace(spin1);
    const Ref<OrbitalSpace>& gg2_space = ggspace(spin2);
    const Ref<OrbitalSpace>& orbs1 = this->orbs(spin1);
    const Ref<OrbitalSpace>& orbs2 = this->orbs(spin2);
    const Ref<OrbitalSpace>& GG1_space = GGspace(spin1);
    const Ref<OrbitalSpace>& GG2_space = GGspace(spin2);
    const Ref<OrbitalSpace>& abs1 = r12world()->abs_space();
    const Ref<OrbitalSpace>& abs2 = r12world()->abs_space();
    const Ref<OrbitalSpace>& cabs1 = r12world()->cabs_space(spin1);
    const Ref<OrbitalSpace>& cabs2 = r12world()->cabs_space(spin2);

    const bool gg1_eq_gg2 = (gg1_space==gg2_space);
    const bool GG1_eq_GG2 = (GG1_space==GG2_space);
    const bool orbs1_eq_orbs2 = (orbs1==orbs2);

    if(gg1_eq_gg2 ^ GG1_eq_GG2) {
      throw ProgrammingError("R12IntEval::contrib_to_VXB_c_ansatz1 -- gg1 and gg2 space must be of the same structure as GG1 and GG2 space",__FILE__,__LINE__);
    }

    // are particles 1 and 2 equivalent?
    const bool part1_equiv_part2 =  spincase2 != AlphaBeta || orbs1_eq_orbs2;
    ExEnv::out0() << "part1_equiv_part2 = " << ((part1_equiv_part2) ? "true" : "false") << endl;
    // Need to antisymmetrize 1 and 2
    const bool antisymmetrize = spincase2 != AlphaBeta;

    // some transforms can be skipped if gg1/gg2 is a subset of GG1/GG2
    // for now it's always true since can only use pq products to generate geminals
    const bool gg12_in_GG12 = true;

    Ref<TwoParticleContraction> tpcontract;
    const R12Technology::ABSMethod absmethod = r12world()->r12tech()->abs_method();
    tpcontract = new CABS_OBS_Contraction(orbs1->rank());

    const bool cabs_method = (absmethod ==  R12Technology::ABS_CABS ||
                    absmethod == R12Technology::ABS_CABSPlus);

    Ref<OrbitalSpace> rispace1, rispace2;
    if (cabs_method) {
      rispace1 = r12world()->ribs_space();
      rispace2 = r12world()->ribs_space();
    }
    else {
      throw InputError("R12IntEval::contrib_to_VXB_c_ansatz1 -- Ansatz 1 only implemented with CABS type methods.",__FILE__,__LINE__);
    }

    /// computing intermediate V
    { /// OBS term
      Ref<TwoParticleContraction> tpcontract_orbs1_orbs2 = new Direct_Contraction(orbs1->rank(),orbs2->rank(),-1.0);
      std::vector<std::string> tforms_f12;
      {
        R12TwoBodyIntKeyCreator tformkey_creator(
                    moints_runtime4(),
                    GG1_space,orbs1,GG2_space,orbs2,
                    corrfactor(),true
                    );
        fill_container(tformkey_creator,tforms_f12);
      }
      std::vector<std::string> tforms;
      if(!gg12_in_GG12) {
        const std::string tform_key = ParsedTwoBodyFourCenterIntKey::key(gg1_space->id(),gg2_space->id(),
                                                                         orbs1->id(),orbs2->id(),
                                                                         std::string("ERI"),
                                                                         std::string(TwoBodyIntLayout::b1b2_k1k2));
        tforms.push_back(tform_key);
      }
      else {
        tforms.push_back(tforms_f12[0]);
      }
      contract_tbint_tensor<ManyBodyTensors::I_to_T,ManyBodyTensors::I_to_T,ManyBodyTensors::I_to_T,true,false,false>
          (V_[s], corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_eri(),
       GG1_space, GG2_space, orbs1, orbs2,
       gg1_space, gg2_space, orbs1, orbs2,
       tpcontract_orbs1_orbs2, spincase2!=AlphaBeta, tforms_f12, tforms);
    }
    if (debug_ >= DefaultPrintThresholds::O4) {
        V_[s].print(prepend_spincase(static_cast<SpinCase2>(s),"V(diag+OBS) contribution").c_str());
    }
    { /// CABS term
      Ref<TwoParticleContraction> tpcontract_orbs1_cabs2 = new Direct_Contraction(orbs1->rank(),
                                                                                  cabs2->rank(),
                                                                                  part1_equiv_part2 ? -2.0 : -1.0);
      std::vector<std::string> tforms_f12;
      {
        R12TwoBodyIntKeyCreator tformkey_creator(
                    moints_runtime4(),
                    GG1_space,orbs1,GG2_space,cabs2,
                    corrfactor(),true
                    );
        fill_container(tformkey_creator,tforms_f12);
      }
      std::vector<std::string> tforms;
      if(!gg12_in_GG12) {
        const std::string tform_key = ParsedTwoBodyFourCenterIntKey::key(gg1_space->id(),gg2_space->id(),
                                                                         orbs1->id(),cabs2->id(),
                                                                         std::string("ERI"),
                                                                         std::string(TwoBodyIntLayout::b1b2_k1k2));
        tforms.push_back(tform_key);
      }
      else {
        tforms.push_back(tforms_f12[0]);
      }
      contract_tbint_tensor<ManyBodyTensors::I_to_T,ManyBodyTensors::I_to_T,ManyBodyTensors::I_to_T,true,false,false>
          (V_[s], corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_eri(),
       GG1_space, GG2_space, orbs1, cabs2,
       gg1_space, gg2_space, orbs1, cabs2,
       tpcontract_orbs1_cabs2, spincase2!=AlphaBeta, tforms_f12, tforms);
    }
    if(!part1_equiv_part2) { /// to be added to the second term
      Ref<TwoParticleContraction> tpcontract_cabs1_orbs2 = new Direct_Contraction(cabs1->rank(),orbs2->rank(),-1.0);
      std::vector<std::string> tforms_f12;
      {
        R12TwoBodyIntKeyCreator tformkey_creator(
                    moints_runtime4(),
                    GG1_space,cabs1,GG2_space,orbs2,
                    corrfactor(),true
                    );
        fill_container(tformkey_creator,tforms_f12);
      }
      std::vector<std::string> tforms;
      if(!gg12_in_GG12) {
        const std::string tform_key = ParsedTwoBodyFourCenterIntKey::key(gg1_space->id(),gg2_space->id(),
                                                                         cabs1->id(),orbs2->id(),
                                                                         std::string("ERI"),
                                                                         std::string(TwoBodyIntLayout::b1b2_k1k2));
        tforms.push_back(tform_key);
      }
      else {
        tforms.push_back(tforms_f12[0]);
      }
      contract_tbint_tensor<ManyBodyTensors::I_to_T,ManyBodyTensors::I_to_T,ManyBodyTensors::I_to_T,true,false,false>
          (V_[s], corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_eri(),
       GG1_space, GG2_space, cabs1, orbs2,
       gg1_space, gg2_space, cabs1, orbs2,
       tpcontract_cabs1_orbs2, spincase2!=AlphaBeta, tforms_f12, tforms);
    }  // endif !part1_equiv_part2
    if (!antisymmetrize && part1_equiv_part2) {
        symmetrize<false>(V_[s],V_[s],GG1_space,gg1_space);
    }
    if (debug_ >= DefaultPrintThresholds::O4) {
        ExEnv::out0() << indent << __FILE__ << ": "<<__LINE__<<"\n";
        V_[s].print(prepend_spincase(static_cast<SpinCase2>(s),"V(diag+OBS+ABS) contribution").c_str());
    }

    /// computing intermediate X
    { /// OBS term
      Ref<TwoParticleContraction> tpcontract_orbs1_orbs2 = new Direct_Contraction(orbs1->rank(),orbs2->rank(),-1.0);
      std::vector<std::string> tforms_f12;
      {
        R12TwoBodyIntKeyCreator tformkey_creator(
                    moints_runtime4(),
                    GG1_space,orbs1,GG2_space,orbs2,
                    corrfactor(),true
                    );
        fill_container(tformkey_creator,tforms_f12);
      }
      contract_tbint_tensor<ManyBodyTensors::I_to_T,ManyBodyTensors::I_to_T,ManyBodyTensors::I_to_T,true,true,false>
          (X_[s], corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_f12(),
       GG1_space, GG2_space, orbs1, orbs2,
       GG1_space, GG2_space, orbs1, orbs2,
       tpcontract_orbs1_orbs2, spincase2!=AlphaBeta, tforms_f12, tforms_f12);
    }
    if (debug_ >= DefaultPrintThresholds::O4) {
        X_[s].print(prepend_spincase(static_cast<SpinCase2>(s),"X(diag+OBS) contribution").c_str());
    }
    { /// CABS term
      Ref<TwoParticleContraction> tpcontract_orbs1_cabs2 = new Direct_Contraction(orbs1->rank(),
                                                                                  cabs2->rank(),
                                                                                  part1_equiv_part2 ? -2.0 : -1.0);
      std::vector<std::string> tforms_f12;
      {
        R12TwoBodyIntKeyCreator tformkey_creator(
                    moints_runtime4(),
                    GG1_space,orbs1,GG2_space,cabs2,
                    corrfactor(),true
                    );
        fill_container(tformkey_creator,tforms_f12);
      }
      contract_tbint_tensor<ManyBodyTensors::I_to_T,ManyBodyTensors::I_to_T,ManyBodyTensors::I_to_T,true,true,false>
          (X_[s], corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_f12(),
       GG1_space, GG2_space, orbs1, cabs2,
       GG1_space, GG2_space, orbs1, cabs2,
       tpcontract_orbs1_cabs2, spincase2!=AlphaBeta, tforms_f12, tforms_f12);
    }
    if(!part1_equiv_part2) { /// to be added to the second term
      Ref<TwoParticleContraction> tpcontract_cabs1_orbs2 = new Direct_Contraction(cabs1->rank(),orbs2->rank(),-1.0);
      std::vector<std::string> tforms_f12;
      {
        R12TwoBodyIntKeyCreator tformkey_creator(
                    moints_runtime4(),
                    GG1_space,cabs1,GG2_space,orbs2,
                    corrfactor(),true
                    );
        fill_container(tformkey_creator,tforms_f12);
      }
      contract_tbint_tensor<ManyBodyTensors::I_to_T,ManyBodyTensors::I_to_T,ManyBodyTensors::I_to_T,true,true,false>
          (X_[s], corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_f12(),
       GG1_space, GG2_space, cabs1, orbs2,
       GG1_space, GG2_space, cabs1, orbs2,
       tpcontract_cabs1_orbs2, spincase2!=AlphaBeta, tforms_f12, tforms_f12);
    } // endif !part1_equiv_part2
    if (!antisymmetrize && part1_equiv_part2) {
        symmetrize<false>(X_[s],X_[s],GG1_space,GG1_space);
    }
    if (debug_ >= DefaultPrintThresholds::O4) {
        X_[s].print(prepend_spincase(static_cast<SpinCase2>(s),"X(diag+OBS+ABS) contribution").c_str());
    }

  }

  timer.exit();
}
