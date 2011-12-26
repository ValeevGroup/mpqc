//
// dkh.cc
//
// Copyright (C) 2008 Edward Valeev
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
#include <chemistry/qc/basis/petite.h>

// tests against analytic results from Mathematica suggest that M contribution comes out too large by a factor of 2
// TODO: figure out where we made an apparent error for a factor of 2
#define P_INCLUDES_Mover2 1

using namespace std;
using namespace sc;

namespace sc {
namespace detail {

  // produces kinetic energy-transformed intspace
  Ref<OrbitalSpace> t_x_X(const Ref<OrbitalSpace>& extspace,
                          const Ref<OrbitalSpace>& intspace) {

    if (!extspace->integral()->equiv(intspace->integral()))
      throw ProgrammingError("two OrbitalSpaces use incompatible Integral factories");
    const Ref<GaussianBasisSet> bs1 = extspace->basis();
    const Ref<GaussianBasisSet> bs2 = intspace->basis();
    const bool bs1_eq_bs2 = (bs1 == bs2);
    int nshell1 = bs1->nshell();
    int nshell2 = bs2->nshell();

    RefSCMatrix vec1t = extspace->coefs().t();
    RefSCMatrix vec2 = intspace->coefs();

    Ref<Integral> localints = extspace->integral()->clone();
    localints->set_basis(bs1,bs2);
#if 1
    Ref<SCElementOp> op = new OneBodyIntOp(localints->kinetic());
#else
    Ref<SCElementOp> op = new OneBodyIntOp(localints->overlap());
#endif
    RefSCDimension aodim1 = vec1t.coldim();
    RefSCDimension aodim2 = vec2.rowdim();
    Ref<SCMatrixKit> aokit = bs1->so_matrixkit();
    RefSCMatrix t12_ao(aodim1, aodim2, aokit);
    t12_ao.assign(0.0);
    t12_ao.element_op(op);
    op = 0;

    // finally, transform
    RefSCMatrix t12 = vec1t * t12_ao * vec2;
    t12_ao = 0;

    RefSCMatrix t12_coefs = intspace->coefs() * t12.t();
    std::string id = extspace->id();  id += "_T(";  id += intspace->id();  id += ")";
    ExEnv::out0() << "id = " << id << endl;
    std::string name = "(T)-weighted space";
    Ref<OrbitalSpace> result = new OrbitalSpace(id, name, extspace, t12_coefs, bs2);

    return result;
  }

}}

void R12IntEval::compute_B_DKH_() {
  if (evaluated_)
    return;

  // verify that this contribution should be included
  const R12Technology::H0_dk_approx_pauli H0_dk_approx_pauli = r12world()->r12tech()->H0_dk_approx();
  if (H0_dk_approx_pauli == R12Technology::H0_dk_approx_pauli_false)
    return;

  const bool obs_eq_vbs = r12world()->obs_eq_vbs();
  const bool obs_eq_ribs = r12world()->obs_eq_ribs();
  const unsigned int maxnabs = r12world()->r12tech()->maxnabs();

  // Check if the requested calculation is implemented
  if (!obs_eq_vbs && maxnabs < 1)
    throw FeatureNotImplemented("OBS!=VBS & maxnabs == 0 is not supported yet in relativistic calculations",__FILE__,__LINE__);

  Timer tim_B_DKH2("Analytic double-commutator B(DKH2) intermediate");
  ExEnv::out0() << endl << indent
      << "Entered analytic double-commutator B(DKH2) intermediate evaluator" << endl;
  ExEnv::out0() << incindent;

  const double c = 137.0359895;
  const double minus_one_over_8c2 = -1.0 / (8.0 * c * c);

  //
  // Compute kinetic energy integrals and obtain geminal-generator spaces transformed with them
  //
  Ref<OrbitalSpace> t_x_P[NSpinCases1];
  Ref<OrbitalSpace> rispace = (maxnabs < 1) ? r12world()->refwfn()->orbs() : r12world()->ribs_space();
  for(int s=0; s<NSpinCases1; ++s) {
    
    const SpinCase1 spin = static_cast<SpinCase1>(s);
    Ref<OrbitalSpace> x = GGspace(spin);

    t_x_P[s] = sc::detail::t_x_X(x, rispace);
    this->r12world()->world()->tfactory()->orbital_registry()->add(make_keyspace_pair(t_x_P[s]));
  }

  // Loop over every 2-e spincase
  for (int s=0; s<nspincases2(); s++) {
    
    const SpinCase2 spincase2 = static_cast<SpinCase2>(s);
    const SpinCase1 spin1 = case1(spincase2);
    const SpinCase1 spin2 = case2(spincase2);
    Ref<RefWavefunction> ref = r12world()->refwfn();

    const Ref<OrbitalSpace>& GGspace1 = GGspace(spin1);
    const Ref<OrbitalSpace>& GGspace2 = GGspace(spin2);
    const bool x1_eq_x2 = (GGspace1 == GGspace2);

    // are particles 1 and 2 equivalent?
    const bool part1_equiv_part2 = (spincase2 != AlphaBeta) || x1_eq_x2;
    // Need to antisymmetrize 1 and 2
    const bool antisymmetrize = (spincase2 != AlphaBeta);

    RefSCMatrix B_DKH = B_[s].clone();
    B_DKH.assign(0.0);

    Ref<OrbitalSpace> t_x1 = t_x_P[spin1];
    Ref<OrbitalSpace> t_x2 = t_x_P[spin2];

    // <xy|z T> tforms
    std::vector<std::string> tforms_xyzT;
    {
      R12TwoBodyIntKeyCreator tformkey_creator(moints_runtime4(), GGspace1, GGspace1, GGspace2,
          t_x2, corrfactor(), true, true);
      fill_container(tformkey_creator, tforms_xyzT);
    }
    compute_tbint_tensor<ManyBodyTensors::I_to_T, true, true>(
        B_DKH,
        corrfactor()->tbint_type_f12t1f12(),
        GGspace1, GGspace1,
        GGspace2, t_x2,
        antisymmetrize,
        tforms_xyzT);
    if (!part1_equiv_part2) {
      // <xy|T z> tforms
      std::vector<std::string> tforms_xyTz;
      {
        R12TwoBodyIntKeyCreator tformkey_creator(moints_runtime4(), GGspace1, t_x1, GGspace2,
            GGspace2, corrfactor(), true, true);
        fill_container(tformkey_creator, tforms_xyTz);
      }
      compute_tbint_tensor<ManyBodyTensors::I_to_T, true, true>(
          B_DKH,
          corrfactor()->tbint_type_f12t1f12(),
          GGspace1, t_x1,
          GGspace2, GGspace2,
          antisymmetrize,
          tforms_xyTz);
    }
    else {
      B_DKH.scale(2.0);
      if (spincase2 == AlphaBeta) {
        symmetrize<false>(B_DKH,B_DKH,GGspace1,GGspace1);
      }
    }

    // symmetrize bra and ket
    B_DKH.scale(0.5);
    RefSCMatrix B_DKH_t = B_DKH.t();
    B_DKH.accumulate(B_DKH_t);  B_DKH_t = 0;

    // and scale by the prefactor
    // M1 = 2 ( f12(T1+T2)f12 (T1 + T2) + (T1 + T2) f12(T1+T2)f12 ) = 4 * ( f12T1f12 (T1 + T2) + (T1 + T2) f12T1f12 )
    // what I have computed so far is 1/2 * (f12T1f12 (T1 + T2) + (T1 + T2) f12T1f12)
    // hence multiply by 8
    B_DKH.scale(8.0 * minus_one_over_8c2);

#if P_INCLUDES_Mover2
    B_DKH.scale(0.5);
#endif

    if (debug_ >= DefaultPrintThresholds::O4) {
      B_DKH.print(prepend_spincase(spincase2,"B(DKH2) contribution (M1)").c_str());
    }
    B_[s].accumulate(B_DKH);
    B_DKH.assign(0.0);

    // Transforms for this type of integrals does not yet exists
    // will not use RangeCreator here because its input is hardwired to corrfactor()->tbintdescr()
    // will create a set of descriptors so that compute_tbint_tensor can construct transforms
    // only 1 correlation factor of G12 type is accepted at the moment, hence completely manual work here
    std::vector<std::string> tforms_g12dkh;
    if (corrfactor()->nfunctions() > 1)
      throw FeatureNotImplemented("B(DKH2) evaluator can only work with one correlation factor",__FILE__,__LINE__);
    Ref<R12Technology::G12CorrelationFactor> g12corrfact;   g12corrfact << corrfactor();
    Ref<R12Technology::G12NCCorrelationFactor> g12nccorrfact; g12nccorrfact << corrfactor();
    if (g12nccorrfact.null() && g12corrfact.null())
      throw FeatureNotImplemented("B(DKH2) evaluator can only work with Gaussian (or Gaussian-expanded) correlation factors",__FILE__,__LINE__);
    Ref<IntParamsG12> params = g12corrfact.nonnull() ? new IntParamsG12(g12corrfact->function(0),
                                                                      g12corrfact->function(0)) :
                                                     new IntParamsG12(g12nccorrfact->function(0),
                                                                      g12nccorrfact->function(0));
    Ref<TwoBodyIntDescr> descr_g12dkh = new TwoBodyIntDescrG12DKH(r12world()->integral(), params);
    const std::string descr_key = moints_runtime4()->descr_key(descr_g12dkh);
    const std::string tform_key = ParsedTwoBodyFourCenterIntKey::key(GGspace1->id(),GGspace2->id(),
                                                           GGspace1->id(),GGspace2->id(),
                                                           descr_key,
                                                           std::string(TwoBodyIntLayout::b1b2_k1k2));
    tforms_g12dkh.push_back(tform_key);

    // M2H + M3H
    compute_tbint_tensor<ManyBodyTensors::I_to_T, true, true>(
                                                              B_DKH,
                                                              TwoBodyOper::g12p4g12_m_g12t1g12t1,
                                                              GGspace1, GGspace1,
                                                              GGspace2, GGspace2,
                                                              antisymmetrize,
                                                              tforms_g12dkh);
    B_DKH.scale(minus_one_over_8c2);
#if P_INCLUDES_Mover2
    B_DKH.scale(0.5);
#endif
    if (debug_ >= DefaultPrintThresholds::O4) {
      B_DKH.print(prepend_spincase(spincase2,"B(DKH2) contribution (M2+M3)").c_str());
    }
    B_[s].accumulate(B_DKH);
    B_DKH.assign(0.0);

  } // end of spincase2 loop

  ExEnv::out0() << decindent;
  ExEnv::out0() << indent << "Exited analytic double-commutator B(DKH2) intermediate evaluator"
      << endl;

  tim_B_DKH2.exit();
  checkpoint_();

  return;
}

void R12IntEval::contrib_to_B_DKH_a_() {
  if (evaluated_)
    return;

  const bool obs_eq_vbs = r12world()->obs_eq_vbs();
  const bool obs_eq_ribs = r12world()->obs_eq_ribs();

  // verify that this contribution should be included
  const R12Technology::H0_dk_approx_pauli H0_dk_approx_pauli = r12world()->r12tech()->H0_dk_approx();
  if (H0_dk_approx_pauli == R12Technology::H0_dk_approx_pauli_false)
    return;

  // Check if the requested calculation is implemented
  if (obs_eq_ribs)
    throw FeatureNotImplemented("OBS==RIBS is not supported yet in relativistic calculations",__FILE__,__LINE__);

  Timer tim_B_DKH2("Analytic single-commutator B(DKH2) intermediate");
  ExEnv::out0() << endl << indent
      << "Entered analytic single-commutator B(DKH2) intermediate evaluator" << endl;
  ExEnv::out0() << incindent;

  const double c = 137.0359895;
  const double minus_one_over_8c2 = -1.0 / (8.0 * c * c);

  //
  // Compute kinetic energy integrals and x_t(A) spaces
  //
  Ref<OrbitalSpace> t_x_P[NSpinCases1];
  Ref<OrbitalSpace> t_p_P[NSpinCases1];
  Ref<OrbitalSpace> t_m_P[NSpinCases1];
  Ref<OrbitalSpace> t_A_P[NSpinCases1];
  for(int s=0; s<nspincases1(); ++s) {
    
    const SpinCase1 spin = static_cast<SpinCase1>(s);
    Ref<OrbitalSpace> x = GGspace(spin);
    Ref<OrbitalSpace> m = occ(spin);
    Ref<OrbitalSpace> obs = r12world()->refwfn()->orbs(spin);
    Ref<OrbitalSpace> ribs = r12world()->ribs_space();
    Ref<OrbitalSpace> cabs = r12world()->cabs_space(spin);

    Ref<OrbitalSpaceRegistry> oreg = this->r12world()->world()->tfactory()->orbital_registry();

    t_x_P[s] = sc::detail::t_x_X(x, ribs);
    oreg->add(make_keyspace_pair(t_x_P[s]));

    t_m_P[s] = sc::detail::t_x_X(m, ribs);
    oreg->add(make_keyspace_pair(t_m_P[s]));

    t_p_P[s] = sc::detail::t_x_X(obs, ribs);
    oreg->add(make_keyspace_pair(t_p_P[s]));

    t_A_P[s] = sc::detail::t_x_X(cabs, ribs);
    oreg->add(make_keyspace_pair(t_A_P[s]));
  }
  if (t_x_P[Beta].null()) t_x_P[Beta] = t_x_P[Alpha];
  if (t_m_P[Beta].null()) t_m_P[Beta] = t_m_P[Alpha];
  if (t_p_P[Beta].null()) t_p_P[Beta] = t_p_P[Alpha];
  if (t_A_P[Beta].null()) t_A_P[Beta] = t_A_P[Alpha];

  // Loop over every 2-e spincase
  for (int s=0; s<nspincases2(); s++) {
    using namespace sc::mbptr12;
    const SpinCase2 spincase2 = static_cast<SpinCase2>(s);
    const SpinCase1 spin1 = case1(spincase2);
    const SpinCase1 spin2 = case2(spincase2);
    Ref<RefWavefunction> ref = r12world()->refwfn();

    const Ref<OrbitalSpace>& GGspace1 = GGspace(spin1);
    const Ref<OrbitalSpace>& GGspace2 = GGspace(spin2);
    const Ref<OrbitalSpace>& x_tP1 = t_x_P[spin1];
    const Ref<OrbitalSpace>& x_tP2 = t_x_P[spin2];
    const bool x1_eq_x2 = (GGspace1 == GGspace2);

    // are particles 1 and 2 equivalent?
    const bool part1_equiv_part2 = (spincase2 != AlphaBeta) || x1_eq_x2;
    // Need to antisymmetrize 1 and 2
    const bool antisymmetrize = (spincase2 != AlphaBeta);

    RefSCMatrix B_DKH = B_[s].clone();
    B_DKH.assign(0.0);

    { // pq contribution
      const Ref<OrbitalSpace>& p_tP1 = t_p_P[spin1];
      const Ref<OrbitalSpace>& p_tP2 = t_p_P[spin2];
      const Ref<OrbitalSpace>& orbs1 = ref->orbs(spin1);
      const Ref<OrbitalSpace>& orbs2 = ref->orbs(spin2);

    Ref<TwoParticleContraction> tpcontract = new CABS_OBS_Contraction(ref->orbs(spin1)->rank());
    // <x y|p q> tforms
    std::vector<std::string> tforms_xy_pq;
    {
      R12TwoBodyIntKeyCreator tformkey_creator(
        moints_runtime4(),
        GGspace1,
        orbs1,
        GGspace2,
        orbs2,
        corrfactor(),true
        );
      fill_container(tformkey_creator, tforms_xy_pq);
    }
    // <x y|p_T q> tforms
    std::vector<std::string> tforms_xy_pTq;
    {
      R12TwoBodyIntKeyCreator tformkey_creator(
        moints_runtime4(),
        GGspace1,
        p_tP1,
        GGspace2,
        orbs2,
        corrfactor(),true
        );
      fill_container(tformkey_creator, tforms_xy_pTq);
    }
    // <x y|p q_T> tforms
    std::vector<std::string> tforms_xy_pqT;
    {
      R12TwoBodyIntKeyCreator tformkey_creator(
        moints_runtime4(),
        GGspace1,
        orbs1,
        GGspace2,
        p_tP2,
        corrfactor(),true
        );
      fill_container(tformkey_creator, tforms_xy_pqT);
    }
    // <x_T y|p q> tforms
    std::vector<std::string> tforms_xTy_pq;
    {
      R12TwoBodyIntKeyCreator tformkey_creator(
        moints_runtime4(),
        x_tP1,
        orbs1,
        GGspace2,
        orbs2,
        corrfactor(),true
        );
      fill_container(tformkey_creator, tforms_xTy_pq);
    }
    // <x y_T|p q> tforms
    std::vector<std::string> tforms_xyT_pq;
    {
      R12TwoBodyIntKeyCreator tformkey_creator(
        moints_runtime4(),
        GGspace1,
        orbs1,
        x_tP2,
        orbs2,
        corrfactor(),true
        );
      fill_container(tformkey_creator, tforms_xyT_pq);
    }

    contract_tbint_tensor<ManyBodyTensors::I_to_T,
        ManyBodyTensors::I_to_T,
        ManyBodyTensors::I_to_T,
        true,true,false>
        (
        B_DKH, corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_t1f12(),
        GGspace1, GGspace2,
        orbs1, orbs2,
        GGspace1, GGspace2,
        p_tP1, orbs2,
        tpcontract,
        spincase2!=AlphaBeta, tforms_xy_pq, tforms_xy_pTq
        );
    contract_tbint_tensor<ManyBodyTensors::I_to_T,
        ManyBodyTensors::I_to_T,
        ManyBodyTensors::I_to_T,
        true,true,false>
        (
        B_DKH, corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_t2f12(),
        GGspace1, GGspace2,
        orbs1, orbs2,
        GGspace1, GGspace2,
        orbs1, p_tP2,
        tpcontract,
        spincase2!=AlphaBeta, tforms_xy_pq, tforms_xy_pqT
        );

    contract_tbint_tensor<ManyBodyTensors::I_to_T,
        ManyBodyTensors::I_to_T,
        ManyBodyTensors::I_to_T,
        true,true,false>
        (
        B_DKH, corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_t1f12(),
        GGspace1, GGspace2,
        orbs1, orbs2,
        x_tP1, GGspace2,
        orbs1, orbs2,
        tpcontract,
        spincase2!=AlphaBeta, tforms_xy_pq, tforms_xTy_pq
        );
    contract_tbint_tensor<ManyBodyTensors::I_to_T,
        ManyBodyTensors::I_to_T,
        ManyBodyTensors::I_to_T,
        true,true,false>
        (
        B_DKH, corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_t2f12(),
        GGspace1, GGspace2,
        orbs1, orbs2,
        GGspace1, x_tP2,
        orbs1, orbs2,
        tpcontract,
        spincase2!=AlphaBeta, tforms_xy_pq, tforms_xyT_pq
        );

    // symmetrize bra and ket
    B_DKH.scale(0.5);
    RefSCMatrix B_DKH_t = B_DKH.t();
    B_DKH.accumulate(B_DKH_t);  B_DKH_t = 0;

    // 0.5 due to double counting contributions from electrons 1 and 2
    // 4 from T^2 -> p^4
    B_DKH.scale(0.5 * 4.0 * minus_one_over_8c2);

    if (debug_ >= DefaultPrintThresholds::O4) {
      globally_sum_intermeds_();
      B_DKH.print(prepend_spincase(static_cast<SpinCase2>(s),"B(DKH pq) contribution").c_str());
    }

    B_[s].accumulate(B_DKH);
    B_DKH.assign(0.0);
    } // end of pq contribution

    { // mA contribution
      const Ref<OrbitalSpace>& occ1 = occ(spin1);
      const Ref<OrbitalSpace>& occ2 = occ(spin2);
      const Ref<OrbitalSpace>& cabs1 = r12world()->cabs_space(spin1);
      const Ref<OrbitalSpace>& cabs2 = r12world()->cabs_space(spin2);
      const Ref<OrbitalSpace>& m_tP1 = t_m_P[spin1];
      const Ref<OrbitalSpace>& m_tP2 = t_m_P[spin2];
      const Ref<OrbitalSpace>& A_tP1 = t_A_P[spin1];
      const Ref<OrbitalSpace>& A_tP2 = t_A_P[spin2];

      // TODO take advantage of the permutational symmetry
      Ref<TwoParticleContraction> tpcontract_mA =
        new Direct_Contraction(occ1->rank(), cabs2->rank(), -1.0);
      Ref<TwoParticleContraction> tpcontract_Am =
        new Direct_Contraction(cabs1->rank(), occ2->rank(), -1.0);

      // <x y|m A> tforms
      std::vector<std::string> tforms_xy_mA;
      {
        R12TwoBodyIntKeyCreator tformkey_creator(
          moints_runtime4(),
          GGspace1,
          occ1,
          GGspace2,
          cabs2,
          corrfactor(),true
          );
        fill_container(tformkey_creator, tforms_xy_mA);
      }
      // <x y|m_T A> tforms
      std::vector<std::string> tforms_xy_mTA;
      {
        R12TwoBodyIntKeyCreator tformkey_creator(
          moints_runtime4(),
          GGspace1,
          m_tP1,
          GGspace2,
          cabs2,
          corrfactor(),true
          );
        fill_container(tformkey_creator, tforms_xy_mTA);
      }
      // <x y|m A_T> tforms
      std::vector<std::string> tforms_xy_mAT;
      {
        R12TwoBodyIntKeyCreator tformkey_creator(
          moints_runtime4(),
          GGspace1,
          occ1,
          GGspace2,
          A_tP2,
          corrfactor(),true
          );
        fill_container(tformkey_creator, tforms_xy_mAT);
      }
      // <x_T y|m A> tforms
      std::vector<std::string> tforms_xTy_mA;
      {
        R12TwoBodyIntKeyCreator tformkey_creator(
          moints_runtime4(),
          x_tP1,
          occ1,
          GGspace2,
          cabs2,
          corrfactor(),true
          );
        fill_container(tformkey_creator, tforms_xTy_mA);
      }
      // <x y_T|m A> tforms
      std::vector<std::string> tforms_xyT_mA;
      {
        R12TwoBodyIntKeyCreator tformkey_creator(
          moints_runtime4(),
          GGspace1,
          occ1,
          x_tP2,
          cabs2,
          corrfactor(),true
          );
        fill_container(tformkey_creator, tforms_xyT_mA);
      }

      // <x y|A m> tforms
      std::vector<std::string> tforms_xy_Am;
      {
        R12TwoBodyIntKeyCreator tformkey_creator(
          moints_runtime4(),
          GGspace1,
          cabs1,
          GGspace2,
          occ2,
          corrfactor(),true
          );
        fill_container(tformkey_creator, tforms_xy_Am);
      }
      // <x y|A m_T> tforms
      std::vector<std::string> tforms_xy_AmT;
      {
        R12TwoBodyIntKeyCreator tformkey_creator(
          moints_runtime4(),
          GGspace1,
          cabs1,
          GGspace2,
          m_tP2,
          corrfactor(),true
          );
        fill_container(tformkey_creator, tforms_xy_AmT);
      }
      // <x y|A_T m> tforms
      std::vector<std::string> tforms_xy_ATm;
      {
        R12TwoBodyIntKeyCreator tformkey_creator(
          moints_runtime4(),
          GGspace1,
          A_tP1,
          GGspace2,
          occ2,
          corrfactor(),true
          );
        fill_container(tformkey_creator, tforms_xy_ATm);
      }
      // <x_T y|A m> tforms
      std::vector<std::string> tforms_xTy_Am;
      {
        R12TwoBodyIntKeyCreator tformkey_creator(
          moints_runtime4(),
          x_tP1,
          cabs1,
          GGspace2,
          occ2,
          corrfactor(),true
          );
        fill_container(tformkey_creator, tforms_xTy_Am);
      }
      // <x y_T|A m> tforms
      std::vector<std::string> tforms_xyT_Am;
      {
        R12TwoBodyIntKeyCreator tformkey_creator(
          moints_runtime4(),
          GGspace1,
          cabs1,
          x_tP2,
          occ2,
          corrfactor(),true
          );
        fill_container(tformkey_creator, tforms_xyT_Am);
      }

      contract_tbint_tensor<ManyBodyTensors::I_to_T,
          ManyBodyTensors::I_to_T,
          ManyBodyTensors::I_to_T,
          true,true,false>
          (
          B_DKH, corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_t1f12(),
          GGspace1, GGspace2,
          occ1, cabs2,
          GGspace1, GGspace2,
          m_tP1, cabs2,
          tpcontract_mA,
          spincase2!=AlphaBeta, tforms_xy_mA, tforms_xy_mTA
          );
      contract_tbint_tensor<ManyBodyTensors::I_to_T,
          ManyBodyTensors::I_to_T,
          ManyBodyTensors::I_to_T,
          true,true,false>
          (
          B_DKH, corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_t2f12(),
          GGspace1, GGspace2,
          occ1, cabs2,
          GGspace1, GGspace2,
          occ1, A_tP2,
          tpcontract_mA,
          spincase2!=AlphaBeta, tforms_xy_mA, tforms_xy_mAT
          );
      contract_tbint_tensor<ManyBodyTensors::I_to_T,
          ManyBodyTensors::I_to_T,
          ManyBodyTensors::I_to_T,
          true,true,false>
          (
          B_DKH, corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_t1f12(),
          GGspace1, GGspace2,
          occ1, cabs2,
          x_tP1, GGspace2,
          occ1, cabs2,
          tpcontract_mA,
          spincase2!=AlphaBeta, tforms_xy_mA, tforms_xTy_mA
          );
      contract_tbint_tensor<ManyBodyTensors::I_to_T,
          ManyBodyTensors::I_to_T,
          ManyBodyTensors::I_to_T,
          true,true,false>
          (
          B_DKH, corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_t2f12(),
          GGspace1, GGspace2,
          occ1, cabs2,
          GGspace1, x_tP2,
          occ1, cabs2,
          tpcontract_mA,
          spincase2!=AlphaBeta, tforms_xy_mA, tforms_xyT_mA
          );

          contract_tbint_tensor<ManyBodyTensors::I_to_T,
              ManyBodyTensors::I_to_T,
              ManyBodyTensors::I_to_T,
              true,true,false>
              (
              B_DKH, corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_t2f12(),
              GGspace1, GGspace2,
              cabs1, occ2,
              GGspace1, GGspace2,
              cabs1, m_tP2,
              tpcontract_Am,
              spincase2!=AlphaBeta, tforms_xy_Am, tforms_xy_AmT
              );
          contract_tbint_tensor<ManyBodyTensors::I_to_T,
              ManyBodyTensors::I_to_T,
              ManyBodyTensors::I_to_T,
              true,true,false>
              (
              B_DKH, corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_t1f12(),
              GGspace1, GGspace2,
              cabs1, occ2,
              GGspace1, GGspace2,
              A_tP1, occ2,
              tpcontract_Am,
              spincase2!=AlphaBeta, tforms_xy_Am, tforms_xy_ATm
              );
          contract_tbint_tensor<ManyBodyTensors::I_to_T,
              ManyBodyTensors::I_to_T,
              ManyBodyTensors::I_to_T,
              true,true,false>
              (
              B_DKH, corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_t1f12(),
              GGspace1, GGspace2,
              cabs1, occ2,
              x_tP1, GGspace2,
              cabs1, occ2,
              tpcontract_Am,
              spincase2!=AlphaBeta, tforms_xy_Am, tforms_xTy_Am
              );
          contract_tbint_tensor<ManyBodyTensors::I_to_T,
              ManyBodyTensors::I_to_T,
              ManyBodyTensors::I_to_T,
              true,true,false>
              (
              B_DKH, corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_t2f12(),
              GGspace1, GGspace2,
              cabs1, occ2,
              GGspace1, x_tP2,
              cabs1, occ2,
              tpcontract_Am,
              spincase2!=AlphaBeta, tforms_xy_Am, tforms_xyT_Am
              );

    // symmetrize bra and ket
    B_DKH.scale(0.5);
    RefSCMatrix B_DKH_t = B_DKH.t();
    B_DKH.accumulate(B_DKH_t);  B_DKH_t = 0;

    // 0.5 due to double counting contributions from electrons 1 and 2
    // 4 from T^2 -> p^4
    B_DKH.scale(0.5 * 4.0 * minus_one_over_8c2);

    if (debug_ >= DefaultPrintThresholds::O4) {
      globally_sum_intermeds_();
      B_DKH.print(prepend_spincase(static_cast<SpinCase2>(s),"B(DKH mA+Am) contribution").c_str());
    }

    B_[s].accumulate(B_DKH);
    B_DKH.assign(0.0);
    } // end of mA contribution

    B_[s].print(prepend_spincase(static_cast<SpinCase2>(s),"B(diag+OBS+ABS+MV) contribution").c_str());
  }

  ExEnv::out0() << decindent;
  ExEnv::out0() << indent << "Exited analytic single-commutator B(DKH2) intermediate evaluator"
      << endl;

  tim_B_DKH2.exit();
  checkpoint_();

  return;
}

////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
