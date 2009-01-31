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

void R12IntEval::compute_B_DKH_() {
  if (evaluated_)
    return;

  const bool obs_eq_vbs = r12info_->basis_vir()->equiv(r12info_->basis());
  const bool obs_eq_ribs = r12info()->basis_ri()->equiv(r12info()->basis());

  // Check if the requsted calculation is implemented
  if (!obs_eq_vbs)
    throw FeatureNotImplemented("OBS!=VBS is not supported yet in relativistic calculations",__FILE__,__LINE__);

  Timer tim_B_DKH2("Analytic B(DKH2) intermediate");
  ExEnv::out0() << endl << indent
      << "Entered analytic B(DKH2) intermediate evaluator" << endl;
  ExEnv::out0() << incindent;

  const double c = 137.0359895;
  const double minus_one_over_8c2 = -1.0 / (8.0 * c * c);
  //
  // Compute kinetic energy integrals and obtain geminal-generator spaces transformed with them
  //
  Ref<OrbitalSpace> t_x_P[NSpinCases1];
  // first need kinetic energy matrix between OBS and RIBS
  Ref<GaussianBasisSet> obs = r12info()->basis();
  const unsigned int maxnabs = r12info()->maxnabs();
  Ref<OrbitalSpace> rispace = (maxnabs < 1) ? r12info()->refinfo()->orbs() : r12info()->ribs_space();
  Ref<GaussianBasisSet> ribs = rispace->basis();
  // Symmetry blocked
  RefSCMatrix t_obs_ribs;
  {
    Ref<Integral> ref_integral = r12info()->refinfo()->ref()->integral();
    ref_integral->set_basis(obs, ribs);
    Ref<SCElementOp> op = new OneBodyIntOp(ref_integral->kinetic());
    //Ref<SCElementOp> op = new OneBodyIntOp(ref_integral->overlap());
    RefSCMatrix t_obs_ribs_ao(obs->basisdim(), ribs->basisdim(),
                              ribs->matrixkit());
    t_obs_ribs_ao.assign(0.0);
    t_obs_ribs_ao.element_op(op);
    op = 0;
    //t_obs_ribs_ao.print("T(OBS/RIBS) in AO basis");

    // now must the OBS dimension into the SO basis
    ref_integral->set_basis(ribs);
    Ref<PetiteList> pl_ribs = ref_integral->petite_list();
    RefSCDimension ribs_aodim_sb = pl_ribs->AO_basisdim();
    ref_integral->set_basis(obs);
    Ref<PetiteList> pl_obs = ref_integral->petite_list();
    RefSCDimension obs_aodim_sb = pl_obs->AO_basisdim();
    // first convert to a blocked matrix
    RefSCMatrix blocked_t_obs_ribs_ao(obs_aodim_sb, ribs_aodim_sb,
                                      ribs->so_matrixkit());
    blocked_t_obs_ribs_ao->convert(t_obs_ribs_ao);  t_obs_ribs_ao = 0;
    // now must transform the RIBS dimension with S^-1 = U . Ut, where U = rispace->coefs()
    RefSCMatrix U = rispace->coefs();
    t_obs_ribs = blocked_t_obs_ribs_ao * U * U.t();
    //t_obs_ribs.print("T(OBS/RIBS) * S(RIBS/RIBS)^-1 in AO basis");
  }
  for(int s=0; s<NSpinCases1; ++s) {
    using namespace sc::LinearR12;
    const SpinCase1 spin = static_cast<SpinCase1>(s);
    // t_x_P = xspace.coef() * t_obs_ribs
    Ref<OrbitalSpace> x = xspace(spin);
    if (x->basis() != r12info()->basis())
      throw FeatureNotImplemented("Cannot handle geminal generating spaces supported by non-orbital basis sets",__FILE__,__LINE__);
    RefSCMatrix t_x_P_coefs = t_obs_ribs.t() * x->coefs();
    //x->coefs().print("xspace coefficients");
    //t_x_P_coefs.print("T-weighted xspace coefficients");
    std::string id = x->id();  id += "_T(";  id += rispace->id();  id += ")";
    ExEnv::out0() << "id = " << id << endl;
    std::string name = "(T)-weighted space";
    name = prepend_spincase(spin,name);
    t_x_P[s] = new OrbitalSpace(id, name, x, t_x_P_coefs, ribs);
    OrbitalSpaceRegistry::instance()->add(make_keyspace_pair(t_x_P[s]));
  }

  //
  // Compute transformed integrals
  //

  // Loop over every 2-e spincase
  for (int s=0; s<nspincases2(); s++) {
    using namespace sc::LinearR12;
    const SpinCase2 spincase2 = static_cast<SpinCase2>(s);
    const SpinCase1 spin1 = case1(spincase2);
    const SpinCase1 spin2 = case2(spincase2);
    Ref<SingleRefInfo> refinfo = r12info()->refinfo();

    const Ref<OrbitalSpace>& xspace1 = xspace(spin1);
    const Ref<OrbitalSpace>& xspace2 = xspace(spin2);
    const bool x1_eq_x2 = (xspace1 == xspace2);

    // are particles 1 and 2 equivalent?
    const bool part1_equiv_part2 = (spincase2 != AlphaBeta) || x1_eq_x2;
    // Need to antisymmetrize 1 and 2
    const bool antisymmetrize = (spincase2 != AlphaBeta);

    RefSCMatrix B_DKH = B_[s].clone();
    B_DKH.assign(0.0);

    Ref<OrbitalSpace> t_x1 = t_x_P[spin1];
    Ref<OrbitalSpace> t_x2 = t_x_P[spin2];

    // <xy|T z> tforms
    std::vector<std::string> tforms_xyTz;
    {
      R12TwoBodyIntKeyCreator tformkey_creator(r12info()->moints_runtime(), xspace1, t_x1, xspace2,
          xspace2, r12info()->corrfactor(), true, true);
      fill_container(tformkey_creator, tforms_xyTz);
    }
    compute_tbint_tensor<ManyBodyTensors::I_to_T, true, true>(
        B_DKH,
        corrfactor()->tbint_type_f12t1f12(),
        xspace1, t_x1,
        xspace2, xspace2,
        antisymmetrize,
        tforms_xyTz);
    if (!part1_equiv_part2) {
      // <xy|z T> tforms
      std::vector<std::string> tforms_xyzT;
      {
        R12TwoBodyIntKeyCreator tformkey_creator(r12info()->moints_runtime(), xspace1, xspace1, xspace2,
            t_x2, r12info()->corrfactor(), true, true);
        fill_container(tformkey_creator, tforms_xyzT);
      }
      compute_tbint_tensor<ManyBodyTensors::I_to_T, true, true>(
          B_DKH,
          corrfactor()->tbint_type_f12t1f12(),
          xspace1,
          xspace1,
          xspace2, t_x2,
          antisymmetrize,
          tforms_xyzT);
    }
    else {
      B_DKH.scale(2.0);
      symmetrize<false>(B_DKH,B_DKH,xspace1,xspace1);
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
    {
      Ref<G12NCCorrelationFactor> g12corrfact;
      g12corrfact << corrfactor();
      if (g12corrfact.null())
        throw FeatureNotImplemented("B(DKH2) evaluator can only work in standard approximation C",__FILE__,__LINE__);
      if (g12corrfact->nfunctions() > 1)
        throw FeatureNotImplemented("B(DKH2) evaluator can only work with one correlation factor",__FILE__,__LINE__);
      Ref<TwoBodyIntDescr> descr_g12dkh = new TwoBodyIntDescrG12DKH(r12info()->integral(),
                                                                    new IntParamsG12(g12corrfact->function(0),
                                                                                     g12corrfact->function(0))
      );
      const std::string descr_key = r12info()->moints_runtime()->descr_key(descr_g12dkh);
      const std::string tform_key = ParsedTwoBodyIntKey::key(xspace1->id(),xspace2->id(),
                                                             xspace1->id(),xspace2->id(),
                                                             descr_key,
                                                             std::string(MOIntsRuntime::Layout_b1b2_k1k2));
      tforms_g12dkh.push_back(tform_key);
    }

    // M2H + M3H
    compute_tbint_tensor<ManyBodyTensors::I_to_T, true, true>(
                                                              B_DKH,
                                                              TwoBodyInt::g12p4g12_m_g12t1g12t1,
                                                              xspace1, xspace1,
                                                              xspace2, xspace2,
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
  ExEnv::out0() << indent << "Exited analytic B(DKH2) intermediate evaluator"
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
