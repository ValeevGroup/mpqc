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

#include <cassert>
#include <math/scmat/local.h>
#include <chemistry/qc/mbptr12/r12int_eval.h>
#include <chemistry/qc/mbptr12/creator.h>
#include <chemistry/qc/mbptr12/container.h>
#include <chemistry/qc/mbptr12/compute_tbint_tensor.h>
#include <chemistry/qc/mbptr12/contract_tbint_tensor.h>
#include <math/scmat/util.h>

using namespace std;
using namespace sc;

namespace {
  void _print(SpinCase2 spin,
              const Ref<DistArray4>& mat,
              const char* label);
}

RefSCMatrix
R12IntEval::V(SpinCase2 spincase2,
              const Ref<OrbitalSpace>& p1,
              const Ref<OrbitalSpace>& p2)
{
  return this->V_cabs(spincase2, p1, p2);
}



RefSCMatrix
R12IntEval::V_cabs(SpinCase2 spincase2,
                   const Ref<OrbitalSpace>& p1,
                   const Ref<OrbitalSpace>& p2)
{
  Ref<LocalSCMatrixKit> local_matrix_kit = new LocalSCMatrixKit();

  const bool obs_eq_vbs = r12world()->obs_eq_vbs();
  const bool obs_eq_ribs = r12world()->obs_eq_ribs();

  const bool p1_eq_p2 = (p1 == p2);
  // are particles 1 and 2 equivalent?
  const bool part1_equiv_part2 =  (spincase2 != AlphaBeta || p1_eq_p2);
  // Need to antisymmetrize 1 and 2
  const bool antisymmetrize = (spincase2 != AlphaBeta);
  //if (!p1_eq_p2 && antisymmetrize)
  //  throw FeatureNotImplemented("R12IntEval::V() -- p1 == p2 must be true if AA or BB spin case",__FILE__,__LINE__);

  const SpinCase1 spin1 = case1(spincase2);
  const SpinCase1 spin2 = case2(spincase2);

  RefSCMatrix V;
  const unsigned int np1 = p1->rank();
  const unsigned int np2 = p2->rank();
  if (!spin_polarized() && (spincase2 == AlphaAlpha || spincase2 == BetaBeta)) {
    V = local_matrix_kit->matrix(dim_f12(spincase2),
                 new SCDimension(np1*(np1-1)/2));
    RefSCMatrix Vab = this->V(AlphaBeta,p1,p1);
    sc::antisymmetrize(V,Vab,GGspace(Alpha),p1);
    return V;
  }

  Timer tim("R12 intermeds (tensor contract): Vpqxy");

  const Ref<OrbitalSpace>& x1 = GGspace(spin1);
  const Ref<OrbitalSpace>& x2 = GGspace(spin2);
  const Ref<OrbitalSpace>& orbs1 = orbs(spin1);
  const Ref<OrbitalSpace>& orbs2 = orbs(spin2);

  // some transforms can be skipped if p1/p2 equals x1/x2
  const bool p1p2_eq_x1x2 = (p1 == x1) && (p2 == x2);

  // allocate V
  const unsigned int np12 = (p1_eq_p2 && spincase2 != AlphaBeta) ? (np1*(np1-1)/2) : (np1*np2);
  RefSCDimension dim_p12 = new SCDimension(np12);
  V = local_matrix_kit->matrix(dim_f12(spincase2), dim_p12);
  V.assign(0.0);

  // The diagonal contribution
  Ref<R12Technology::G12CorrelationFactor> g12ptr; g12ptr << corrfactor();
  Ref<R12Technology::G12NCCorrelationFactor> g12ncptr; g12ncptr << corrfactor();
  Ref<R12Technology::R12CorrelationFactor> r12ptr; r12ptr << corrfactor();
  if (r12ptr) {
    RefSCMatrix I = compute_I_(x1,x2,p1,p2);
    if (!antisymmetrize)
      V.accumulate(I);
    else
      sc::antisymmetrize<true>(V,I,x1,x2,p1,p2);
  }
  else if (g12ptr || g12ncptr) {
    std::vector<std::string> tforms_f12_xmyn;
    {
      R12TwoBodyIntKeyCreator tformkey_creator(moints_runtime4(),
                    x1,p1,
                    x2,p2,
                    corrfactor(),true
                    );
      fill_container(tformkey_creator,tforms_f12_xmyn);
    }
    compute_tbint_tensor<ManyBodyTensors::I_to_T,true,false>(
      V, corrfactor()->tbint_type_f12eri(),
      x1, p1,
      x2, p2,
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
                      x1,
                      orbs1,
                      x2,
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

  RefSCMatrix Vobs = V.clone();
  Vobs.assign(0.0);
  contract_tbint_tensor<true,false>
    (
     Vobs, corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_eri(),
     -1.0,
     x1, x2,
     orbs1, orbs2,
     p1, p2,
     orbs1, orbs2,
     spincase2!=AlphaBeta, tforms_f12, tforms
     );
  V.accumulate(Vobs);

  if (debug_ >= DefaultPrintThresholds::O4) {
    V.print(prepend_spincase(spincase2,"Vpqxy: diag+OBS contribution").c_str());
  }

  if (debug_ >= DefaultPrintThresholds::O4) {
    std::vector< Ref<DistArray4> > vobs_da4;
    contract_tbint_tensor<true,false>
      (
       vobs_da4, corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_eri(),
       -1.0,
       x1, x2,
       orbs1, orbs2,
       p1, p2,
       orbs1, orbs2,
       spincase2!=AlphaBeta, tforms_f12, tforms
       );

    RefSCMatrix Vobs_da4 = Vobs.clone();
    Vobs_da4.assign(0.0);
    Vobs_da4 << vobs_da4[0];
    print_scmat_norms(Vobs_da4-Vobs, prepend_spincase(spincase2,"Vpqxy: OBS contribution distarray4-incore (should be 0)").c_str());
  }

  // These terms only contribute if Projector=2
  if (!obs_eq_ribs && ansatz()->projector() == R12Technology::Projector_2) {

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
                    x1,
                    occ1,
                    x2,
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

    RefSCMatrix Vcabs1 = V.clone(); Vcabs1.assign(0.0);
    contract_tbint_tensor<true,false>
      (
       Vcabs1, corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_eri(),
       perm_factor,
       x1, x2,
       occ1, rispace2,
       p1, p2,
       occ1, rispace2,
       antisymmetrize, tforms_f12_xmyP, tforms_imjP
       );
    V.accumulate(Vcabs1);

    if (debug_ >= DefaultPrintThresholds::O4) {
      std::vector< Ref<DistArray4> > vcabs1_da4;
      contract_tbint_tensor<true,false>
      (
       vcabs1_da4, corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_eri(),
       perm_factor,
       x1, x2,
       occ1, rispace2,
       p1, p2,
       occ1, rispace2,
       antisymmetrize, tforms_f12_xmyP, tforms_imjP
       );

      RefSCMatrix Vcabs1_da4 = Vcabs1.clone();
      Vcabs1_da4.assign(0.0);
      Vcabs1_da4 << vcabs1_da4[0];
      print_scmat_norms(Vcabs1_da4 - Vcabs1, prepend_spincase(spincase2,"Vpqxy: CABS1 contribution distarray4-incore (should be 0)").c_str());
    }

    // If particles 1 and 2 are not equivalent, also need another set of terms
    if (!part1_equiv_part2) {

      std::vector<std::string> tforms_iPjm;
      std::vector<std::string> tforms_f12_xPym;
      {
    R12TwoBodyIntKeyCreator tformkey_creator(
                      moints_runtime4(),
                      x1,
                      rispace1,
                      x2,
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

      RefSCMatrix Vcabs2 = V.clone(); Vcabs2.assign(0.0);
      contract_tbint_tensor<true,false>
    (
     Vcabs2, corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_eri(),
     -1.0,
     x1, x2,
     rispace1, occ2,
     p1, p2,
     rispace1, occ2,
     antisymmetrize, tforms_f12_xPym, tforms_iPjm
     );
      V.accumulate(Vcabs2);

      if (debug_ >= DefaultPrintThresholds::O4) {
        std::vector< Ref<DistArray4> > vcabs2_da4;
        contract_tbint_tensor<true,false>
        (
         vcabs2_da4, corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_eri(),
         -1.0,
         x1, x2,
         rispace1, occ2,
         p1, p2,
         rispace1, occ2,
         antisymmetrize, tforms_f12_xPym, tforms_iPjm
         );

        RefSCMatrix Vcabs2_da4 = Vcabs2.clone();
        Vcabs2_da4.assign(0.0);
        Vcabs2_da4 << vcabs2_da4[0];
        print_scmat_norms(Vcabs2_da4 - Vcabs2, prepend_spincase(spincase2,"Vpqxy: CABS2 contribution distarray4-incore (should be 0)").c_str());
      }
    } // if part1_equiv_part2
  } // ABS != OBS

  if (!antisymmetrize && part1_equiv_part2) {
    symmetrize<false>(V,V,x1,p1);
  }

  if (debug_ >= DefaultPrintThresholds::O4) {
    V.print(prepend_spincase(spincase2,"Vpqxy: diag+OBS+ABS contribution").c_str());
  }

  return V;
}


RefSCMatrix
R12IntEval::V_genref_spinfree(const Ref<OrbitalSpace>& p1,
                   const Ref<OrbitalSpace>& p2)
{
  MPQC_ASSERT(r12world()->spinadapted());
  ExEnv::out0() << std::endl << std::endl << indent << "Entered V_genref_spinfree\n\n";

  const bool debugprint = false;
  Ref<LocalSCMatrixKit> local_matrix_kit = new LocalSCMatrixKit();

  const bool p1_eq_p2 = (p1 == p2);
  const unsigned int np1 = p1->rank();
  const unsigned int np2 = p2->rank();
  const Ref<OrbitalSpace>& GG1 = GGspace(Alpha);
  const Ref<OrbitalSpace>& GG2 = GGspace(Alpha);
  const Ref<OrbitalSpace>& obs = orbs(Alpha);
  const Ref<OrbitalSpace>& occspace = occ(Alpha);
  const Ref<OrbitalSpace>& cabs = r12world()->cabs_space(Alpha);
  const Ref<OrbitalSpace>& g_m_m_av = gamma_m_m_av(); //'p' menas obs; 'A': cabs; 'P': cbs; 'm': occ
  const unsigned int np12 = np1*np2;
  RefSCDimension dim_p12 = new SCDimension(np12);

  RefSCMatrix V = local_matrix_kit->matrix(dim_f12(AlphaBeta), dim_p12);
  V.assign(0.0);
  Timer tim("R12 intermeds (tensor contract): genref spinfree Vpqxy");
  {
    std::vector<std::string> tforms_f12_xmyn;
    {
      R12TwoBodyIntKeyCreator tformkey_creator(moints_runtime4(),
                    GG1,p1,
                    GG2,p2, corrfactor(),true);
      fill_container(tformkey_creator,tforms_f12_xmyn);
    }
    compute_tbint_tensor<ManyBodyTensors::I_to_T,true,false>(V, corrfactor()->tbint_type_f12eri(),
      GG1, p1,
      GG2, p2, false, tforms_f12_xmyn);

    if (debug_ >= DefaultPrintThresholds::O4 || debugprint)
      V.print(std::string("Vpqxy spin free genref: diag contribution").c_str());
  }
  {
    std::vector<std::string> tforms;
    std::vector<std::string> tforms_f12;
    {
      R12TwoBodyIntKeyCreator tformkey_creator(moints_runtime4(),
                        GG1,obs,
                        GG2,obs,corrfactor(),true);
      fill_container(tformkey_creator,tforms_f12);
    }
    const std::string tform_key = ParsedTwoBodyFourCenterIntKey::key(p1->id(),p2->id(), obs->id(),obs->id(),
                                                                     std::string("ERI"), std::string(TwoBodyIntLayout::b1b2_k1k2));
    tforms.push_back(tform_key);

    RefSCMatrix Vobs = V.clone();
    Vobs.assign(0.0);
    contract_tbint_tensor<true,false>
      (Vobs, corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_eri(),
       -1.0,
       GG1,GG2,obs,obs,
       p1, p2, obs,obs,
       false, tforms_f12, tforms);
    V.accumulate(Vobs);

    if (debug_ >= DefaultPrintThresholds::O4 || debugprint) {
      V.print(prepend_spincase(AlphaBeta,"Vpqxy spinfree genref: diag+OBS contribution").c_str());
    }
  }

  {
    std::vector<std::string> tforms_f12;
    {
      R12TwoBodyIntKeyCreator tform_creator(moints_runtime4(),
                          GG1,g_m_m_av,GG2,cabs, corrfactor(),true);
      fill_container(tform_creator,tforms_f12);
    }
    std::vector<std::string> tforms;
    {
      const std::string tform_key = ParsedTwoBodyFourCenterIntKey::key(p1->id(),p2->id(),
                                                                       occspace->id(),cabs->id(),
                                                                       std::string("ERI"),
                                                                       std::string(TwoBodyIntLayout::b1b2_k1k2));
      tforms.push_back(tform_key);
    }
    // 2.0 due to using "average" 1-RDM (hence half of spin-free 1-RDM)
    // obviously this form is not symmetric w.r.t p1 <-> p2
    // will symmetry below just for cleanliness
    contract_tbint_tensor<true,false>(V, corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_eri(),
        -2.0,
        GG1,GG2,g_m_m_av,cabs,
        p1, p2,occspace,cabs,
        false, tforms_f12, tforms);
  }

  if(p1_eq_p2) symmetrize<false>(V,V,GG1, p1);
  tim.exit();
  if (debug_ >= DefaultPrintThresholds::O4 || debugprint) V.print(prepend_spincase(AlphaBeta,"Vpqxy genref spinfree: diag+OBS+CABS contribution").c_str());
  ExEnv::out0() << "\n\n" << indent << "Exited V_genref_spinfree\n\n";
  return V;
}


std::vector<Ref<DistArray4> > R12IntEval::V_distarray4(
                                                       SpinCase2 spincase2,
                                                       const Ref<OrbitalSpace>& p1,
                                                       const Ref<OrbitalSpace>& p2) {

  const bool include_diag =  true; // set to false to skip the "diagonal" (f12/r12) constribution
  const bool include_obs  =  true; // set to false to skip the P1 P2 part of the projector
  const bool include_cabs =  true; // set to false to skip the A1 O2 + O1 A2 part of the projector

  const bool obs_eq_ribs = r12world()->obs_eq_ribs();
  const bool obs_eq_vbs = r12world()->obs_eq_vbs();

  const bool p1_eq_p2 = (p1 == p2);
  // are particles 1 and 2 equivalent?
  const bool part1_equiv_part2 = spincase2 != AlphaBeta || p1_eq_p2;
  // Need to antisymmetrize 1 and 2
  const bool antisymmetrize = spincase2 != AlphaBeta;

  const SpinCase1 spin1 = case1(spincase2);
  const SpinCase1 spin2 = case2(spincase2);

  std::vector<Ref<DistArray4> > V;
  if (!spin_polarized() && (spincase2 == AlphaAlpha || spincase2 == BetaBeta)) {
    const unsigned int nx1 = GGspace(spin1)->rank();
    const unsigned int np1 = p1->rank();
    V = this->V_distarray4(AlphaBeta, p1, p1);
    for (int s = 0; s < V.size(); ++s)
      sc::antisymmetrize(V[s]);
    return V;
  }

  Timer tim("R12 intermeds (tensor contract): Vpqxy");

  const Ref<OrbitalSpace>& x1 = GGspace(spin1);
  const Ref<OrbitalSpace>& x2 = GGspace(spin2);
  const Ref<OrbitalSpace>& orbs1 = orbs(spin1);
  const Ref<OrbitalSpace>& orbs2 = orbs(spin2);

  // some transforms can be skipped if p1/p2 equals x1/x2
  const bool p1p2_eq_x1x2 = (p1 == x1) && (p2 == x2);

  // The diagonal contribution
  if (include_diag) {
  Ref<R12Technology::G12CorrelationFactor> g12ptr;
  g12ptr << corrfactor();
  Ref<R12Technology::G12NCCorrelationFactor> g12ncptr;
  g12ncptr << corrfactor();
  Ref<R12Technology::R12CorrelationFactor> r12ptr;
  r12ptr << corrfactor();
  if (r12ptr) {
    MPQC_ASSERT(false);
#if 0
    RefSCMatrix I = compute_I_(x1,x2,p1,p2);
    V.accumulate(I);
#endif
  } else if (g12ptr || g12ncptr) {
    std::vector<std::string> tforms_f12_xmyn;
    {
      R12TwoBodyIntKeyCreator tformkey_creator(moints_runtime4(), x1, p1,
                                               x2, p2, corrfactor(), true);
      fill_container(tformkey_creator, tforms_f12_xmyn);
    }

    for (int t = 0; t < tforms_f12_xmyn.size(); ++t) {
      Ref<TwoBodyMOIntsTransform> tform =
          moints_runtime4()->get(tforms_f12_xmyn[t]);
      tform->compute();
      // at this point V(diag) is only spatial integrals
      V.push_back(
                  extract(
                          tform->ints_distarray4(),
                          tform->intdescr()->intset(
                                                    corrfactor()->tbint_type_f12eri())));

      // may need to antisymmetrize
      if (antisymmetrize) {
          sc::antisymmetrize(V.back());
      }
    }
  }
  if (debug_ >= DefaultPrintThresholds::allO2N2) {
    for (int s = 0; s < V.size(); ++s)
      _print(spincase2, V[s],
             prepend_spincase(spincase2, "Vpqxy: diag contribution").c_str());
  }
  }

  if (include_obs) {
    const bool SplitPQContributions = false;  // set to true to evaluate PQ as AB + AN + MB + MN

    if (not SplitPQContributions) {
      std::vector<std::string> tforms;
      std::vector<std::string> tforms_f12;
      {
        R12TwoBodyIntKeyCreator tformkey_creator(moints_runtime4(), x1, orbs1,
                                                 x2, orbs2, corrfactor(), true);
        fill_container(tformkey_creator, tforms_f12);
      }
      if (!p1p2_eq_x1x2) {
        const std::string tform_key = ParsedTwoBodyFourCenterIntKey::key(
            p1->id(), p2->id(), orbs1->id(), orbs2->id(), std::string("ERI"),
            std::string(TwoBodyIntLayout::b1b2_k1k2));
        tforms.push_back(tform_key);
      } else
        tforms.push_back(tforms_f12[0]);

      contract_tbint_tensor<true, false>(V, corrfactor()->tbint_type_f12(),
                                         corrfactor()->tbint_type_eri(), -1.0,
                                         x1, x2, orbs1, orbs2, p1, p2, orbs1,
                                         orbs2, antisymmetrize, tforms_f12,
                                         tforms);
    } // PQ
    else { // split PQ = AB + AN + MB + MN
      const bool include_AB = true;
      const bool include_AN = true;
      const bool include_MB = true;
      const bool include_MN = true;

      const Ref<OrbitalSpace>& vir1 = vir(spin1);
      const Ref<OrbitalSpace>& vir2 = vir(spin2);
      const Ref<OrbitalSpace>& occ1 = occ(spin1);
      const Ref<OrbitalSpace>& occ2 = occ(spin2);

      if (include_AB) {
        std::vector<std::string> tforms;
        std::vector<std::string> tforms_f12;
        const Ref<OrbitalSpace>& i1 = vir1;
        const Ref<OrbitalSpace>& i2 = vir2;
        {
          R12TwoBodyIntKeyCreator tformkey_creator(moints_runtime4(), x1, i1,
                                                   x2, i2, corrfactor(), true);
          fill_container(tformkey_creator, tforms_f12);
        }
        {
          const std::string tform_key = ParsedTwoBodyFourCenterIntKey::key(
              p1->id(), p2->id(), i1->id(), i2->id(), std::string("ERI"),
              std::string(TwoBodyIntLayout::b1b2_k1k2));
          tforms.push_back(tform_key);
        }
        contract_tbint_tensor<true, false>(V, corrfactor()->tbint_type_f12(),
                                           corrfactor()->tbint_type_eri(), -1.0,
                                           x1, x2, i1, i2,
                                           p1, p2, i1, i2,
                                           antisymmetrize, tforms_f12,
                                           tforms);
      }
      if (include_MN) {
        std::vector<std::string> tforms;
        std::vector<std::string> tforms_f12;
        const Ref<OrbitalSpace>& i1 = occ1;
        const Ref<OrbitalSpace>& i2 = occ2;
        {
          R12TwoBodyIntKeyCreator tformkey_creator(moints_runtime4(), x1, i1,
                                                   x2, i2, corrfactor(), true);
          fill_container(tformkey_creator, tforms_f12);
        }
        {
          const std::string tform_key = ParsedTwoBodyFourCenterIntKey::key(
              p1->id(), p2->id(), i1->id(), i2->id(), std::string("ERI"),
              std::string(TwoBodyIntLayout::b1b2_k1k2));
          tforms.push_back(tform_key);
        }
        contract_tbint_tensor<true, false>(V, corrfactor()->tbint_type_f12(),
                                           corrfactor()->tbint_type_eri(), -1.0,
                                           x1, x2, i1, i2,
                                           p1, p2, i1, i2,
                                           antisymmetrize, tforms_f12,
                                           tforms);
      }
      if (include_MB) {
        std::vector<std::string> tforms;
        std::vector<std::string> tforms_f12;
        const Ref<OrbitalSpace>& i1 = occ1;
        const Ref<OrbitalSpace>& i2 = vir2;
        {
          R12TwoBodyIntKeyCreator tformkey_creator(moints_runtime4(), x1, i1,
                                                   x2, i2, corrfactor(), true);
          fill_container(tformkey_creator, tforms_f12);
        }
        {
          const std::string tform_key = ParsedTwoBodyFourCenterIntKey::key(
              p1->id(), p2->id(), i1->id(), i2->id(), std::string("ERI"),
              std::string(TwoBodyIntLayout::b1b2_k1k2));
          tforms.push_back(tform_key);
        }
        contract_tbint_tensor<true, false>(V, corrfactor()->tbint_type_f12(),
                                           corrfactor()->tbint_type_eri(), -1.0,
                                           x1, x2, i1, i2,
                                           p1, p2, i1, i2,
                                           antisymmetrize, tforms_f12,
                                           tforms);
      }
      if (include_AN) {
        std::vector<std::string> tforms;
        std::vector<std::string> tforms_f12;
        const Ref<OrbitalSpace>& i1 = vir1;
        const Ref<OrbitalSpace>& i2 = occ2;
        {
          R12TwoBodyIntKeyCreator tformkey_creator(moints_runtime4(), x1, i1,
                                                   x2, i2, corrfactor(), true);
          fill_container(tformkey_creator, tforms_f12);
        }
        {
          const std::string tform_key = ParsedTwoBodyFourCenterIntKey::key(
              p1->id(), p2->id(), i1->id(), i2->id(), std::string("ERI"),
              std::string(TwoBodyIntLayout::b1b2_k1k2));
          tforms.push_back(tform_key);
        }
        contract_tbint_tensor<true, false>(V, corrfactor()->tbint_type_f12(),
                                           corrfactor()->tbint_type_eri(), -1.0,
                                           x1, x2, i1, i2,
                                           p1, p2, i1, i2,
                                           antisymmetrize, tforms_f12,
                                           tforms);
      }

    }

    if (debug_ >= DefaultPrintThresholds::allO2N2) {
      for (int s = 0; s < V.size(); ++s)
        _print(
            spincase2,
            V[s],
            prepend_spincase(spincase2, "Vpqxy: diag+OBS contribution").c_str());
    }
  }

  if (include_cabs) {
  // These terms only contribute if Projector=2
  if (!obs_eq_ribs && ansatz()->projector() == R12Technology::Projector_2) {

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
      R12TwoBodyIntKeyCreator tformkey_creator(moints_runtime4(), x1,
                                               occ1, x2, rispace2,
                                               corrfactor(), true);
      fill_container(tformkey_creator, tforms_f12_xmyP);
    }
    if (!p1p2_eq_x1x2) {
      const std::string
          tform_key =
              ParsedTwoBodyFourCenterIntKey::key(
                                                 p1->id(),
                                                 p2->id(),
                                                 occ1->id(),
                                                 rispace2->id(),
                                                 std::string("ERI"),
                                                 std::string(
                                                             TwoBodyIntLayout::b1b2_k1k2));
      tforms_imjP.push_back(tform_key);
    } else
      tforms_imjP.push_back(tforms_f12_xmyP[0]);

    contract_tbint_tensor<true, false> (V, corrfactor()->tbint_type_f12(),
                                        corrfactor()->tbint_type_eri(),
                                        perm_factor, x1, x2, occ1,
                                        rispace2, p1, p2, occ1, rispace2,
                                        antisymmetrize, tforms_f12_xmyP,
                                        tforms_imjP);

    // If particles 1 and 2 are not equivalent, also need another set of terms
    if (!part1_equiv_part2) {

      std::vector<std::string> tforms_iPjm;
      std::vector<std::string> tforms_f12_xPym;
      {
        R12TwoBodyIntKeyCreator tformkey_creator(moints_runtime4(), x1,
                                                 rispace1, x2, occ2,
                                                 corrfactor(), true);
        fill_container(tformkey_creator, tforms_f12_xPym);
      }
      if (!p1p2_eq_x1x2) {
        const std::string
            tform_key =
                ParsedTwoBodyFourCenterIntKey::key(
                                                   p1->id(),
                                                   p2->id(),
                                                   rispace1->id(),
                                                   occ2->id(),
                                                   std::string("ERI"),
                                                   std::string(
                                                               TwoBodyIntLayout::b1b2_k1k2));
        tforms_iPjm.push_back(tform_key);
      } else
        tforms_iPjm.push_back(tforms_f12_xPym[0]);

      contract_tbint_tensor<true, false> (V, corrfactor()->tbint_type_f12(),
                                          corrfactor()->tbint_type_eri(), -1.0,
                                          x1, x2, rispace1, occ2, p1,
                                          p2, rispace1, occ2, antisymmetrize,
                                          tforms_f12_xPym, tforms_iPjm);

    } // if part1_equiv_part2
  } // ABS != OBS

  if (!antisymmetrize && part1_equiv_part2) {
    for(int s=0; s<V.size(); ++s)
      symmetrize(V[s]);
  }

  if (debug_ >= DefaultPrintThresholds::allO2N2) {
    for (int s = 0; s < V.size(); ++s)
      _print(
             spincase2,
             V[s],
             prepend_spincase(spincase2, "Vpqxy: diag+OBS+ABS contribution").c_str());
  }
  }

  return V;
}

std::vector<Ref<DistArray4> >
R12IntEval::U_distarray4(
    SpinCase2 spincase2,
    const Ref<OrbitalSpace>& p1,
    const Ref<OrbitalSpace>& p2) {

  const bool obs_eq_ribs = r12world()->obs_eq_ribs();
  const bool obs_eq_vbs = r12world()->obs_eq_vbs();

  const bool p1_eq_p2 = (p1 == p2);
  // are particles 1 and 2 equivalent?
  const bool part1_equiv_part2 = spincase2 != AlphaBeta || p1_eq_p2;
  // Need to antisymmetrize 1 and 2
  const bool antisymmetrize = spincase2 != AlphaBeta;

  const SpinCase1 spin1 = case1(spincase2);
  const SpinCase1 spin2 = case2(spincase2);

  std::vector<Ref<DistArray4> > U;
  if (!spin_polarized() && (spincase2 == AlphaAlpha || spincase2 == BetaBeta)) {
    const unsigned int nx1 = GGspace(spin1)->rank();
    const unsigned int np1 = p1->rank();
    U = this->U_distarray4(AlphaBeta, p1, p1);
    for (int s = 0; s < U.size(); ++s)
      sc::antisymmetrize(U[s]);
    return U;
  }

  Timer tim("R12 intermeds (tensor contract): Upqxy");

  const Ref<OrbitalSpace>& x1 = GGspace(spin1);
  const Ref<OrbitalSpace>& x2 = GGspace(spin2);
  const Ref<OrbitalSpace>& orbs1 = orbs(spin1);
  const Ref<OrbitalSpace>& orbs2 = orbs(spin2);
  const Ref<OrbitalSpace>& occ1_act = occ_act(spin1);
  const Ref<OrbitalSpace>& occ2_act = occ_act(spin2);
  const Ref<OrbitalSpace>& cabs1 = r12world()->cabs_space(spin1);
  const Ref<OrbitalSpace>& cabs2 = r12world()->cabs_space(spin2);

  // some transforms can be skipped if p1/p2 equals x1/x2
  const bool p1p2_eq_x1x2 = (p1 == x1) && (p2 == x2);

  std::vector<std::string> tforms_x1i1_p1A1;
  std::vector<std::string> tforms_x1i2_p1A2;
  std::vector<std::string> tforms_x2i1_p2A1;
  std::vector<std::string> tforms_x2i2_p2A2;
  std::vector<std::string> tforms_x1p1_i1A1;
  std::vector<std::string> tforms_x2p2_i2A2;

  {
    R12TwoBodyIntKeyCreator
        tformkey_creator(moints_runtime4(), x1, occ1_act, p1, cabs1,
                         corrfactor(), true, false,
                         TwoBodyIntLayout::b1k1_b2k2);
    fill_container(tformkey_creator, tforms_x1i1_p1A1);
  }
  {
    R12TwoBodyIntKeyCreator
        tformkey_creator(moints_runtime4(), x1, occ2_act, p1, cabs2,
                         corrfactor(), true, false,
                         TwoBodyIntLayout::b1k1_b2k2);
    fill_container(tformkey_creator, tforms_x1i2_p1A2);
  }
  {
    R12TwoBodyIntKeyCreator
        tformkey_creator(moints_runtime4(), x2, occ1_act, p2, cabs1,
                         corrfactor(), true, false,
                         TwoBodyIntLayout::b1k1_b2k2);
    fill_container(tformkey_creator, tforms_x2i1_p2A1);
  }
  {
    R12TwoBodyIntKeyCreator
        tformkey_creator(moints_runtime4(), x2, occ2_act, p2, cabs2,
                         corrfactor(), true, false,
                         TwoBodyIntLayout::b1k1_b2k2);
    fill_container(tformkey_creator, tforms_x2i2_p2A2);
  }
  {
    R12TwoBodyIntKeyCreator
        tformkey_creator(moints_runtime4(), x1, p1, occ1_act, cabs1,
                         corrfactor(), true, false,
                         TwoBodyIntLayout::b1b2_k1k2);
    fill_container(tformkey_creator, tforms_x1p1_i1A1);
  }
  {
    R12TwoBodyIntKeyCreator
        tformkey_creator(moints_runtime4(), x2, p2, occ2_act, cabs2,
                         corrfactor(), true, false,
                         TwoBodyIntLayout::b1b2_k1k2);
    fill_container(tformkey_creator, tforms_x2p2_i2A2);
  }

  std::vector<Ref<DistArray4> > U_b1k1_b2k2;
  contract_tbint_tensor<true, false> (U_b1k1_b2k2, corrfactor()->tbint_type_f12(),
                                      corrfactor()->tbint_type_eri(), +1.0,
                                      x1, p1, occ1_act, cabs1,
                                      x2, p2, occ1_act, cabs1, false,
                                      tforms_x1i1_p1A1, tforms_x2i1_p2A1);
  contract_tbint_tensor<true, false> (U_b1k1_b2k2, corrfactor()->tbint_type_f12(),
                                      corrfactor()->tbint_type_eri(), -1.0,
                                      x1, p1, occ1_act, cabs1,
                                      x2, p2, occ1_act, cabs1, false,
                                      tforms_x1p1_i1A1, tforms_x2i1_p2A1);
  contract_tbint_tensor<true, false> (U_b1k1_b2k2, corrfactor()->tbint_type_f12(),
                                      corrfactor()->tbint_type_eri(), +1.0,
                                      x1, p1, occ2_act, cabs2,
                                      x2, p2, occ2_act, cabs2, false,
                                      tforms_x1i2_p1A2, tforms_x2i2_p2A2);
  contract_tbint_tensor<true, false> (U_b1k1_b2k2, corrfactor()->tbint_type_f12(),
                                      corrfactor()->tbint_type_eri(), -1.0,
                                      x1, p1, occ2_act, cabs2,
                                      x2, p2, occ2_act, cabs2, false,
                                      tforms_x1i2_p1A2, tforms_x2p2_i2A2);

  if (part1_equiv_part2 == false) {
    std::vector<Ref<DistArray4> > U_b2k2_b1k1;
    contract_tbint_tensor<true, false> (U_b2k2_b1k1, corrfactor()->tbint_type_f12(),
                                        corrfactor()->tbint_type_eri(), +1.0,
                                        x2, p2, occ1_act, cabs1,
                                        x1, p1, occ1_act, cabs1, false,
                                        tforms_x2i1_p2A1, tforms_x1i1_p1A1);
    contract_tbint_tensor<true, false> (U_b2k2_b1k1, corrfactor()->tbint_type_f12(),
                                        corrfactor()->tbint_type_eri(), -1.0,
                                        x2, p2, occ1_act, cabs1,
                                        x1, p1, occ1_act, cabs1, false,
                                        tforms_x2i1_p2A1, tforms_x1p1_i1A1);
    contract_tbint_tensor<true, false> (U_b2k2_b1k1, corrfactor()->tbint_type_f12(),
                                        corrfactor()->tbint_type_eri(), +1.0,
                                        x2, p2, occ2_act, cabs2,
                                        x1, p1, occ2_act, cabs2, false,
                                        tforms_x2i2_p2A2, tforms_x1i2_p1A2);
    contract_tbint_tensor<true, false> (U_b2k2_b1k1, corrfactor()->tbint_type_f12(),
                                        corrfactor()->tbint_type_eri(), -1.0,
                                        x2, p2, occ2_act, cabs2,
                                        x1, p1, occ2_act, cabs2, false,
                                        tforms_x2p2_i2A2, tforms_x1i2_p1A2);
    // add U_b2k2_b1k1 to U
    MPQC_ASSERT(false);
  }
  else {
    // spin-restricted closed-shell case -- scale U_b1k1_b2k2 by 2
    MPQC_ASSERT(false);
  }
  // add U_b1k1_b2k2 to U
  MPQC_ASSERT(false);

  std::vector<std::string> tforms_x2p1_i2A1;
  std::vector<std::string> tforms_p2x1_i2A1;
  {
    R12TwoBodyIntKeyCreator
        tformkey_creator(moints_runtime4(), p2, occ2_act, x1, cabs1,
                         corrfactor(), true, false,
                         TwoBodyIntLayout::b1b2_k1k2);
    fill_container(tformkey_creator, tforms_p2x1_i2A1);
  }
  {
    R12TwoBodyIntKeyCreator
        tformkey_creator(moints_runtime4(), x2, occ2_act, p1, cabs1,
                         corrfactor(), true, false,
                         TwoBodyIntLayout::b1b2_k1k2);
    fill_container(tformkey_creator, tforms_x2p1_i2A1);
  }

  std::vector<Ref<DistArray4> > U_k2b1_b2k1;
  contract_tbint_tensor<true, false> (U_k2b1_b2k1, corrfactor()->tbint_type_f12(),
                                      corrfactor()->tbint_type_eri(), -1.0,
                                      p2, x1, occ2_act, cabs1,
                                      x2, p1, occ2_act, cabs1, false,
                                      tforms_p2x1_i2A1, tforms_x2p1_i2A1);

  if (part1_equiv_part2 == false) {
    std::vector<std::string> tforms_p1x2_i1A2;
    std::vector<std::string> tforms_x1p2_i1A2;
    {
      R12TwoBodyIntKeyCreator
          tformkey_creator(moints_runtime4(), x1, occ1_act, p2, cabs2,
                           corrfactor(), true, false,
                           TwoBodyIntLayout::b1b2_k1k2);
      fill_container(tformkey_creator, tforms_x1p2_i1A2);
    }
    {
      R12TwoBodyIntKeyCreator
          tformkey_creator(moints_runtime4(), p1, occ1_act, x2, cabs2,
                           corrfactor(), true, false,
                           TwoBodyIntLayout::b1b2_k1k2);
      fill_container(tformkey_creator, tforms_p1x2_i1A2);
    }

    std::vector<Ref<DistArray4> > U_k1b2_b1k2;
    contract_tbint_tensor<true, false> (U_k1b2_b1k2, corrfactor()->tbint_type_f12(),
                                        corrfactor()->tbint_type_eri(), -1.0,
                                        p1, x2, occ1_act, cabs2,
                                        x1, p2, occ1_act, cabs2, false,
                                        tforms_p1x2_i1A2, tforms_x1p2_i1A2);
    // add U_k1b2_b1k2 to U
    MPQC_ASSERT(false);
  }
  else {
    // scale U_k2b1_b2k1 by 2
    MPQC_ASSERT(false);
  }

  // add U_k2b1_b2k1 to U
  MPQC_ASSERT(false);

  if (!antisymmetrize && part1_equiv_part2) {
    for(int s=0; s<U.size(); ++s)
      symmetrize(U[s]);
  }

  if (debug_ >= DefaultPrintThresholds::O4) {
    for (int s = 0; s < U.size(); ++s)
      _print(
             spincase2,
             U[s],
             prepend_spincase(spincase2, "Upqxy").c_str());
  }

  return U;
}

RefSymmSCMatrix
R12IntEval::P(SpinCase2 spincase2)
{
  

  Ref<LocalSCMatrixKit> local_matrix_kit = new LocalSCMatrixKit();

  const bool obs_eq_vbs = r12world()->obs_eq_vbs();
  const bool obs_eq_ribs = r12world()->obs_eq_ribs();

  const SpinCase1 spin1 = case1(spincase2);
  const SpinCase1 spin2 = case2(spincase2);
  const Ref<OrbitalSpace>& x1 = GGspace(spin1);
  const Ref<OrbitalSpace>& x2 = GGspace(spin2);
  const Ref<OrbitalSpace>& orbs1 = orbs(spin1);
  const Ref<OrbitalSpace>& orbs2 = orbs(spin2);
  const Ref<OrbitalSpace>& occ1 = occ(spin1);
  const Ref<OrbitalSpace>& occ2 = occ(spin2);
  const Ref<OrbitalSpace>& cabs1 = r12world()->cabs_space(spin1);
  const Ref<OrbitalSpace>& cabs2 = r12world()->cabs_space(spin2);
  const RefSCDimension dimf12 = dim_f12(spincase2);

  // are particles 1 and 2 equivalent?
  const bool part1_equiv_part2 =  spincase2 != AlphaBeta || (x1 == x2);
  // Need to antisymmetrize 1 and 2
  const bool antisymmetrize = spincase2 != AlphaBeta;
  // For RHF compute alpha-alpha and beta-beta P from alpha-beta P
  if (!spin_polarized() && (spincase2 == AlphaAlpha || spincase2 == BetaBeta)) {
    RefSymmSCMatrix Paa = local_matrix_kit->symmmatrix(dimf12);
    RefSymmSCMatrix Pab = this->P(AlphaBeta);
    sc::antisymmetrize<false>(Paa,Pab,x1);
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
  Ref<R12Technology::G12CorrelationFactor> g12ptr; g12ptr << corrfactor();
  Ref<R12Technology::G12NCCorrelationFactor> g12ncptr; g12ncptr << corrfactor();
  Ref<R12Technology::R12CorrelationFactor> r12ptr; r12ptr << corrfactor();

  //
  // Diagonal contribution: P_{xy}^{wz} = (RGR)_{xy}^{wz}
  //
  if (r12ptr) {
    std::vector<std::string> tforms_f12_xoyw;
    {
      R12TwoBodyIntKeyCreator tformkey_creator(moints_runtime4(),
                    x1,orbs1,
                    x2,orbs2,
                    corrfactor(),
                    true
                    );
      fill_container(tformkey_creator,tforms_f12_xoyw);
    }
    compute_tbint_tensor<ManyBodyTensors::I_to_T,true,false>(
      P, corrfactor()->tbint_type_f12(),
      x1, x1,
      x2, x2,
      antisymmetrize,
      tforms_f12_xoyw);
  }
  else if (g12ptr || g12ncptr) {
    std::vector<std::string> tforms_f12f12_xoyw;
    {
      R12TwoBodyIntKeyCreator tformkey_creator(
        moints_runtime4(),
        x1,
        x1,
        x2,
        x2,
        corrfactor(),true,true
        );
      fill_container(tformkey_creator,tforms_f12f12_xoyw);
    }
    compute_tbint_tensor<ManyBodyTensors::I_to_T,true,true>(
      P, corrfactor()->tbint_type_f12eri(),
      x1, x1,
      x2, x2,
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
  if (r12ptr) {
    RefSCMatrix I = compute_I_(x1,x2,orbs1,orbs2);
    if (!antisymmetrize)
      V_pp.accumulate(I);
    else
      sc::antisymmetrize<true>(V_pp,I,x1,x2,orbs1,orbs2);
  }
  else if (g12ptr || g12ncptr) {
    std::vector<std::string> tforms_f12_xpyq;
    {
      R12TwoBodyIntKeyCreator tformkey_creator(
        moints_runtime4(),
        x1,
        orbs1,
        x2,
        orbs2,
        corrfactor(),true,false
        );
      fill_container(tformkey_creator,tforms_f12_xpyq);
    }
    compute_tbint_tensor<ManyBodyTensors::I_to_T,true,false>(
      V_pp, corrfactor()->tbint_type_f12eri(),
      x1, orbs1,
      x2, orbs2,
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
        x1,
        orbs1,
        x2,
        orbs2,
        corrfactor(),true,false
        );
      fill_container(tformkey_creator,tforms_f12_xpyq);
    }
    compute_tbint_tensor<ManyBodyTensors::I_to_T,true,false>(
      R_pp, corrfactor()->tbint_type_f12(),
      x1, orbs1,
      x2, orbs2,
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
  if (r12ptr) {
    RefSCMatrix I = compute_I_(x1,x2,occ1,cabs2);
    if (!antisymmetrize)
      V_iA.accumulate(I);
    else
      sc::antisymmetrize<true>(V_iA,I,x1,x2,occ1,cabs2);
  }
  else if (g12ptr || g12ncptr) {
    std::vector<std::string> tforms_f12_xiyA;
    {
      R12TwoBodyIntKeyCreator tformkey_creator(
        moints_runtime4(),
        x1,
        occ1,
        x2,
        cabs2,
        corrfactor(),true,false
        );
      fill_container(tformkey_creator,tforms_f12_xiyA);
    }
    compute_tbint_tensor<ManyBodyTensors::I_to_T,true,false>(
      V_iA, corrfactor()->tbint_type_f12eri(),
      x1, occ1,
      x2, cabs2,
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
        x1,
        occ1,
        x2,
        cabs2,
        corrfactor(),true,false
        );
      fill_container(tformkey_creator,tforms_f12_xiyA);
    }
    compute_tbint_tensor<ManyBodyTensors::I_to_T,true,false>(
      R_iA, corrfactor()->tbint_type_f12(),
      x1, occ1,
      x2, cabs2,
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
    if (r12ptr) {
      RefSCMatrix I = compute_I_(x1,x2,cabs1,occ2);
      if (!antisymmetrize)
        V_Ai.accumulate(I);
      else
        sc::antisymmetrize<true>(V_Ai,I,x1,x2,cabs1,occ2);
    }
    else if (g12ptr || g12ncptr) {
      std::vector<std::string> tforms_f12_xAyi;
      {
        R12TwoBodyIntKeyCreator tformkey_creator(
          moints_runtime4(),
          x1,
          cabs1,
          x2,
          occ2,
          corrfactor(),true,false
          );
        fill_container(tformkey_creator,tforms_f12_xAyi);
      }
      compute_tbint_tensor<ManyBodyTensors::I_to_T,true,false>(
        V_Ai, corrfactor()->tbint_type_f12eri(),
        x1, cabs1,
        x2, occ2,
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
          x1,
          cabs1,
          x2,
          occ2,
          corrfactor(),true,false
          );
        fill_container(tformkey_creator,tforms_f12_xAyi);
      }
      compute_tbint_tensor<ManyBodyTensors::I_to_T,true,false>(
        R_Ai, corrfactor()->tbint_type_f12(),
        x1, cabs1,
        x2, occ2,
        antisymmetrize,
        tforms_f12_xAyi);
    }
    V_Ai.scale(-1.0);
    P.accumulate(V_Ai * R_Ai.t());
    V_Ai = 0;  R_Ai = 0;
  }

  // Because we skipped Ai term if particles 1 and 2 are equivalent, must make particles equivalent explicitly
  if (part1_equiv_part2)
    symmetrize<false>(P,P,x1,x1);

  if (debug_ >= DefaultPrintThresholds::O4) {
    P.print(prepend_spincase(spincase2,"Pxyow: diag+OBS+ABS contribution").c_str());
  }

  return to_lower_triangle(P);
}

std::vector< Ref<DistArray4> >
sc::A_distarray4(SpinCase2 spincase2, const Ref<R12IntEval>& r12eval) {

  const SpinCase1 spin1 = case1(spincase2);
  const SpinCase1 spin2 = case2(spincase2);

  Ref<OrbitalSpace> v1 = r12eval->vir_act(spin1);
  Ref<OrbitalSpace> v2 = r12eval->vir_act(spin2);
  Ref<OrbitalSpace> fv1 = r12eval->F_a_A(spin1);
  Ref<OrbitalSpace> fv2 = r12eval->F_a_A(spin2);
  Ref<OrbitalSpace> x1 = r12eval->GGspace(spin1);
  Ref<OrbitalSpace> x2 = r12eval->GGspace(spin2);

  // are particles 1 and 2 equivalent?
  const bool part1_equiv_part2 = (v1==v2 && x1==x2 && fv1==fv2);

  std::vector<std::string> tform4f_keys; // get 1 3 |F12| 2 4_f
  {
    R12TwoBodyIntKeyCreator tform_creator(r12eval->moints_runtime4(),
      x1,  v1,
      x2, fv2,
      r12eval->corrfactor(),
      true);
    fill_container(tform_creator,tform4f_keys);
  }

  std::vector< Ref<DistArray4> > A;

  const double pfac = part1_equiv_part2 ? 2.0 : 1.0;
  for (int t = 0; t < tform4f_keys.size(); ++t) {
    Ref<TwoBodyMOIntsTransform> tform =
        r12eval->moints_runtime4()->get(tform4f_keys[t]);
    tform->compute();
    A.push_back(extract(tform->ints_distarray4(),
                        tform->intdescr()->intset(r12eval->corrfactor()->tbint_type_f12()),
                        pfac));
    // transform computes spatial integrals -- antisymmetrize if necessary
    if (spincase2 != AlphaBeta) sc::antisymmetrize(A.back());
  }

  if (part1_equiv_part2) {
    if (spincase2 == AlphaBeta) // symmetrization is implicit when computing antisymmetric tensors
      for(int s=0; s<A.size(); ++s)
        symmetrize(A[s]);
  }
  else {
    std::vector<std::string> tform2f_keys; // get 1 3 |F12| 2 4_f
    {
      R12TwoBodyIntKeyCreator tform_creator(r12eval->moints_runtime4(),
        x1, fv1,
        x2,  v2,
        r12eval->corrfactor(),
        true);
      fill_container(tform_creator,tform2f_keys);
    }
    for (int t = 0; t < tform4f_keys.size(); ++t) {
      Ref<TwoBodyMOIntsTransform> tform =
          r12eval->moints_runtime4()->get(tform4f_keys[t]);
      tform->compute();
      Ref<DistArray4> A_2 = extract(tform->ints_distarray4(),
                                    tform->intdescr()->intset(r12eval->corrfactor()->tbint_type_f12()));
      axpy(A_2, 1.0, A[t]);
    }
  }

  return A;
}

////////////////////////

namespace {
  void _print(SpinCase2 spin,
             const Ref<DistArray4>& mat,
             const char* label) {
    if (mat->msg()->me() == 0) {
      const size_t nij = (spin != AlphaBeta && mat->ni() == mat->nj()) ? mat->ni() * (mat->ni()-1) / 2 : mat->ni() * mat->nj();
      const size_t nxy = (spin != AlphaBeta && mat->nx() == mat->ny()) ? mat->nx() * (mat->nx()-1) / 2 : mat->nx() * mat->ny();
      RefSCMatrix scmat = SCMatrixKit::default_matrixkit()->matrix(new SCDimension(nij), new SCDimension(nxy));
      scmat << mat;
      scmat.print(label);
    }
  }

} // anonymous namespace
