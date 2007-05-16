//
// intermeds_g12.cc
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

#include <chemistry/qc/mbptr12/r12int_eval.h>
#include <chemistry/qc/mbptr12/creator.h>
#include <chemistry/qc/mbptr12/container.h>
#include <chemistry/qc/mbptr12/compute_tbint_tensor.h>
#include <chemistry/qc/mbptr12/contract_tbint_tensor.h>

using namespace std;
using namespace sc;

void
R12IntEval::init_intermeds_g12_()
{
  if (evaluated_)
    return;

  // get smart ptr to this, but how? Just unmanage it for now
  Ref<R12IntEval> thisref(this);
  const bool obs_eq_vbs = r12info_->basis_vir()->equiv(r12info_->basis());
  const bool obs_eq_ribs = r12info()->basis_ri()->equiv(r12info()->basis());

  tim_enter("\"diagonal\" part of G12 intermediates");
  ExEnv::out0() << endl << indent
		<< "Entered G12 diagonal intermediates evaluator" << endl;
  ExEnv::out0() << incindent;

  // Test new tensor compute function
  for(int s=0; s<nspincases2(); s++) {
      using namespace sc::LinearR12;
      const SpinCase2 spincase2 = static_cast<SpinCase2>(s);
      const SpinCase1 spin1 = case1(spincase2);
      const SpinCase1 spin2 = case2(spincase2);
      Ref<SingleRefInfo> refinfo = r12info()->refinfo();

      const Ref<MOIndexSpace>& occ1 = occ(spin1);
      const Ref<MOIndexSpace>& occ2 = occ(spin2);
      const Ref<MOIndexSpace>& occ1_act = occ_act(spin1);
      const Ref<MOIndexSpace>& occ2_act = occ_act(spin2);
      const Ref<MOIndexSpace>& xspace1 = xspace(spin1);
      const Ref<MOIndexSpace>& xspace2 = xspace(spin2);

      // for now geminal-generating products must have same equivalence as the occupied orbitals
      const bool occ1_eq_occ2 = (occ1 == occ2);
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

      std::vector<  Ref<TwoBodyMOIntsTransform> > tforms_f12f12_xmyn;
      {
	  NewTransformCreator tform_creator(
	      thisref,
	      xspace1,
	      xspace1,
	      xspace2,
	      xspace2,true,true
	      );
	  fill_container(tform_creator,tforms_f12f12_xmyn);
      }
      std::vector<  Ref<TwoBodyMOIntsTransform> > tforms_f12_xmyn;
      {
	  // use xmyn, not xiyj, cause if OBS != VBS the former is needed
	  NewTransformCreator tform_creator(
	      thisref,
	      xspace1,
	      occ1,
	      xspace2,
	      occ2,true
	      );
	  fill_container(tform_creator,tforms_f12_xmyn);
      }

      compute_tbint_tensor<ManyBodyTensors::I_to_T,true,false>(
	  V_[s], corrfactor()->tbint_type_f12eri(),
	  xspace1, occ1_act,
	  xspace2, occ2_act,
	  antisymmetrize,
	  tforms_f12_xmyn
	  );
      // g12*g12' operator
      compute_tbint_tensor<ManyBodyTensors::I_to_T,true,true>(
	  X_[s], corrfactor()->tbint_type_f12f12(),
	  xspace1, xspace1,
	  xspace2, xspace2,
	  antisymmetrize,
	  tforms_f12f12_xmyn
	  );
      // 0.5 ( g12*[T,g12'] + [g12,T]*g12' ) = [g12,[t1,g12']] + 0.5 ( g12*[T,g12'] - g12'*[T,g12] )
      //   = [g12,[t1,g12']] - 0.5 ( g12*[g12',T] - g12'*[g12,T] ) = [g12,[t1,g12']] - 0.5 ( (beta-alpha)/(beta+alpha) * [g12*g12',T] ),
      // where the last step valid is for 2 primitive gaussians only (must be contracted otherwise)
      // The first term is symmetric with respect to permutation of g12 and g12' and computed directly,
      compute_tbint_tensor<ManyBodyTensors::I_to_T,true,true>(
	  B_[s], corrfactor()->tbint_type_f12t1f12(),
	  xspace1, xspace1,
	  xspace2, xspace2,
	  antisymmetrize,
	  tforms_f12f12_xmyn
	  );
      // the second is antisymmetric wrt such permutation and is only needed when number of geminals > 1
      if (r12info()->corrfactor()->nfunctions() > 1) {
	  RefSCMatrix Banti = B_[s].clone(); Banti.assign(0.0);
	  // the handling of the second term differs between standard approximations {A,A',B} and {C}
	  if (stdapprox() != LinearR12::StdApprox_C) {
	      // 1) in standard approximations A and B the commutators are explicitly evaluated:
	      //    the second term is exactly what is computed in G12Libint2 under the name t1f12 and t2f12 when g12!=g12'
	      compute_tbint_tensor<ManyBodyTensors::I_to_T,true,true>(
		  Banti, corrfactor()->tbint_type_t1f12(),
		  xspace1, xspace1,
		  xspace2, xspace2,
		  antisymmetrize,
		  tforms_f12f12_xmyn
		  );
	      compute_tbint_tensor<ManyBodyTensors::I_to_T,true,true>(
		  Banti, corrfactor()->tbint_type_t2f12(),
		  xspace1, xspace1,
		  xspace2, xspace2,
		  antisymmetrize,
		  tforms_f12f12_xmyn
		  );
	  }
	  else {
	      // 2) in standard approximation C the commutators are evaluated via RI:
	      // Firstly, instead of T we can use h+J (it will be used later anyway)
	      // Second let's designate (beta-alpha)/(beta+alpha) * g12*g12 as A12. A12 is computed as anti_f12f12 by G12NCLibint2
	      // [A12,T] = [A12,h+J] = A12 (hJ_1 + hJ_2) -  (hJ_1 + hJ_2) A12
	      Ref<MOIndexSpace> hj_x1 = hj_x_P(spin1);
	      Ref<MOIndexSpace> hj_x2 = hj_x_P(spin2);

	      // <xy|hJ z> tforms
	      std::vector<  Ref<TwoBodyMOIntsTransform> > tforms_xyHz;
	      {
		  NewTransformCreator tform_creator(thisref,xspace1,hj_x1,xspace2,xspace2,true);
		  fill_container(tform_creator,tforms_xyHz);
	      }
	      // <hJ z|xy> tforms
	      std::vector<  Ref<TwoBodyMOIntsTransform> > tforms_Hzxy;
	      {
		  NewTransformCreator tform_creator(thisref,hj_x1,xspace1,xspace2,xspace2,true);
		  fill_container(tform_creator,tforms_Hzxy);
	      }

	      compute_tbint_tensor<ManyBodyTensors::I_to_T,true,true>(
		  Banti, corrfactor()->tbint_type_f12f12_anti(),
		  xspace1, hj_x1,
		  xspace2, xspace2,
		  antisymmetrize,
		  tforms_xyHz
		  );
	      compute_tbint_tensor<ManyBodyTensors::I_to_mT,true,true>(
		  Banti, corrfactor()->tbint_type_f12f12_anti(),
		  hj_x1, xspace1,
		  xspace2, xspace2,
		  antisymmetrize,
		  tforms_Hzxy
		  );

	      if (!x1_eq_x2) {
		  // <xy|z hJ> tforms
		  std::vector<  Ref<TwoBodyMOIntsTransform> > tforms_xyzH;
		  {
		      NewTransformCreator tform_creator(thisref,xspace1,xspace1,xspace2,hj_x2,true);
		      fill_container(tform_creator,tforms_xyzH);
		  }
		  // <z hJ|xy> tforms
		  std::vector<  Ref<TwoBodyMOIntsTransform> > tforms_zHxy;
		  {
		      NewTransformCreator tform_creator(thisref,xspace1,xspace1,hj_x2,xspace2,true);
		      fill_container(tform_creator,tforms_zHxy);
		  }
		  
		  compute_tbint_tensor<ManyBodyTensors::I_to_T,true,true>(
		      Banti, corrfactor()->tbint_type_f12f12_anti(),
		      xspace1, xspace2,
		      xspace2, hj_x2,
		      antisymmetrize,
		      tforms_xyzH
		      );
		  compute_tbint_tensor<ManyBodyTensors::I_to_mT,true,true>(
		      Banti, corrfactor()->tbint_type_f12f12_anti(),
		      xspace1, xspace1,
		      hj_x2, xspace2,
		      antisymmetrize,
		      tforms_zHxy
		      );
	      }
	      else {
		  Banti.scale(2.0);
	      }
	  }
	  Banti.scale(-0.5);
	  B_[s].accumulate(Banti);
      }
      // Finally, copy B to BC, since their "diagonal" parts are the same
      if (stdapprox() == LinearR12::StdApprox_C)
	  BC_[s].assign(B_[s]);
  }

  ExEnv::out0() << decindent;
  ExEnv::out0() << indent << "Exited G12 diagonal intermediates evaluator" << endl;

  tim_exit("\"diagonal\" part of G12 intermediates");
  checkpoint_();
  
  return;
}

////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
