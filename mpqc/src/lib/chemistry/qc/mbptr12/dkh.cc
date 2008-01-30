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

using namespace std;
using namespace sc;

void R12IntEval::compute_B_DKH_() {
  if (evaluated_)
    return;
  
  // get smart ptr to this, but how? Just unmanage it for now
  Ref<R12IntEval> thisref(this);
  const bool obs_eq_vbs = r12info_->basis_vir()->equiv(r12info_->basis());
  const bool obs_eq_ribs = r12info()->basis_ri()->equiv(r12info()->basis());
  
  // Check if the requsted calculation is implemented
  if (!obs_eq_vbs)
    throw FeatureNotImplemented("OBS!=VBS is not supported yet in relativistic calculations",__FILE__,__LINE__);
  
  tim_enter("Analytic B(DKH2) intermediate");
  ExEnv::out0() << endl << indent
      << "Entered analytic B(DKH2) intermediate evaluator" << endl;
  ExEnv::out0() << incindent;
  
  // Test new tensor compute function
  for (int s=0; s<nspincases2(); s++) {
    using namespace sc::LinearR12;
    const SpinCase2 spincase2 = static_cast<SpinCase2>(s);
    const SpinCase1 spin1 = case1(spincase2);
    const SpinCase1 spin2 = case2(spincase2);
    Ref<SingleRefInfo> refinfo = r12info()->refinfo();
    
    const Ref<MOIndexSpace>& xspace1 = xspace(spin1);
    const Ref<MOIndexSpace>& xspace2 = xspace(spin2);
    const bool x1_eq_x2 = (xspace1 == xspace2);
    
    // are particles 1 and 2 equivalent?
    const bool part1_equiv_part2 = (spincase2 != AlphaBeta) || x1_eq_x2;
    // Need to antisymmetrize 1 and 2
    const bool antisymmetrize = (spincase2 != AlphaBeta);

    RefSCMatrix B_DKH = B_[s].clone();
    B_DKH.assign(0.0);

    const unsigned int maxnabs = r12info()->maxnabs();
    Ref<MOIndexSpace> t_x1, t_x2;
    if (maxnabs < 1) {
      //t_x1 = t_x_p(spin1);
      //t_x2 = t_x_p(spin2);
      t_x1 = xspace1;
      t_x2 = xspace2;
    }
    else {
      //t_x1 = t_x_P(spin1);
      //t_x2 = t_x_P(spin2);
      t_x1 = xspace1;
      t_x2 = xspace2;
    }

    // <xy|T z> tforms
    std::vector< Ref<TwoBodyMOIntsTransform> > tforms_xyTz;
    {
      NewTransformCreator tform_creator(thisref, xspace1, t_x1, xspace2,
          xspace2, true, true);
      fill_container(tform_creator, tforms_xyTz);
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
      std::vector< Ref<TwoBodyMOIntsTransform> > tforms_xyzT;
      {
        NewTransformCreator tform_creator(thisref, xspace1, xspace1, xspace2,
            t_x2, true, true);
        fill_container(tform_creator, tforms_xyzT);
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
    // M1 = 3/2 ( f12(T1+T2)f12 (T1 + T2) + (T1 + T2) f12(T1+T2)f12 ) = 3 * ( f12T1f12 (T1 + T2) + (T1 + T2) f12T1f12 )
    // what I have computed so far is 1/2 * (f12T1f12 (T1 + T2) + (T1 + T2) f12T1f12)
    // hence multiply by 6
    B_DKH.scale(6.0);
    
    if (debug_ >= DefaultPrintThresholds::O4) {
      B_DKH.print(prepend_spincase(spincase2,"B(DKH2) contribution (M1)").c_str());
    }
    //B_[s].accumulate(B_DKH);
    B_DKH.assign(0.0);
    
    // Transforms for this type of integrals does not yet exists
    // will not use RangeCreator here because its input is hardwired to corrfactor()->tbintdescr()
    // will create a set of descriptors so that compute_tbint_tensor can construct transforms
    // only 1 correlation factor of G12 type is accepted at the moment, hence completely manual work here
    std::vector< Ref<TwoBodyMOIntsTransform> > tforms_g12dkh;
    std::vector< Ref<TwoBodyIntDescr> > descrs_g12dkh;
    {
      Ref<G12NCCorrelationFactor> g12corrfact;
      g12corrfact << corrfactor();
      if (g12corrfact.null())
        throw FeatureNotImplemented("B(DKH2) evaluator can only work in standard approximation C",__FILE__,__LINE__);
      descrs_g12dkh.push_back(new TwoBodyIntDescrG12DKH(r12info()->integral(),new IntParamsG12(g12corrfact->function(0),
              g12corrfact->function(0)
          )
      ));
    }
    
    // M2H + M3H
    compute_tbint_tensor<ManyBodyTensors::I_to_T, true, true>(
                                                              B_DKH,
                                                              TwoBodyInt::g12p4g12_m_g12t1g12t1,
                                                              xspace1, xspace1,
                                                              xspace2, xspace2,
                                                              antisymmetrize,
                                                              tforms_g12dkh,
                                                              descrs_g12dkh);
    if (debug_ >= DefaultPrintThresholds::O4) {
      B_DKH.print(prepend_spincase(spincase2,"B(DKH2) contribution (M2+M3)").c_str());
    }
    //B_[s].accumulate(B_DKH);
    B_DKH.assign(0.0);
        
  } // end of spincase2 loop

  ExEnv::out0() << decindent;
  ExEnv::out0() << indent << "Exited analytic B(DKH2) intermediate evaluator"
      << endl;
  
  tim_exit("Analytic B(DKH2) intermediate");
  checkpoint_();
  
  return;
}

////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
