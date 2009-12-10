//
// intermeds_g12.cc
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

void
R12IntEval::init_intermeds_g12_()
{
  if (evaluated_)
    return;

  // get smart ptr to this, but how? Just unmanage it for now
  Ref<R12IntEval> thisref(this);
  const bool obs_eq_vbs = r12world()->obs_eq_vbs();
  const bool obs_eq_ribs = r12world()->obs_eq_ribs();

  Timer tim_diagonal("\"diagonal\" part of G12 intermediates");
  ExEnv::out0() << endl << indent
		<< "Entered G12 diagonal intermediates evaluator" << endl;
  ExEnv::out0() << incindent;

  // Test new tensor compute function
  for(int s=0; s<nspincases2(); s++) {
      using namespace sc::LinearR12;
      const SpinCase2 spincase2 = static_cast<SpinCase2>(s);
      const SpinCase1 spin1 = case1(spincase2);
      const SpinCase1 spin2 = case2(spincase2);

      const Ref<OrbitalSpace>& occ1 = occ(spin1);
      const Ref<OrbitalSpace>& occ2 = occ(spin2);
      const Ref<OrbitalSpace>& occ1_act = occ_act(spin1);
      const Ref<OrbitalSpace>& occ2_act = occ_act(spin2);
      const Ref<OrbitalSpace>& xspace1 = xspace(spin1);
      const Ref<OrbitalSpace>& xspace2 = xspace(spin2);
      const Ref<OrbitalSpace>& gg1space = ggspace(spin1);
      const Ref<OrbitalSpace>& gg2space = ggspace(spin2);
      const Ref<OrbitalSpace>& GG1space = GGspace(spin1);
      const Ref<OrbitalSpace>& GG2space = GGspace(spin2);
      
      // for now geminal-generating products must have same equivalence as the occupied orbitals
      const bool occ1_eq_occ2 = (occ1 == occ2);
      const bool x1_eq_x2 = (xspace1 == xspace2);
      const bool gg1_eq_gg2 = (gg1space == gg2space);
      const bool GG1_eq_GG2 = (GG1space == GG2space);
      if(gg1_eq_gg2 ^ GG1_eq_GG2) {
        throw ProgrammingError("R12IntEval::init_intermeds_g12_ -- this orbital_product cannot be handled yet",__FILE__,__LINE__);
      }

      // are particles 1 and 2 equivalent?
      const bool part1_equiv_part2 =  spincase2 != AlphaBeta || gg1_eq_gg2;
      // Need to antisymmetrize 1 and 2
      const bool antisymmetrize = spincase2 != AlphaBeta;

      // some transforms can be skipped if occ1/occ2 is a subset of x1/x2
      // for now it's always true since can only use ij and pq products to generate geminals
      const bool occ12_in_x12 = true;
      const bool gg1_in_gg2 = true;

      std::vector<std::string> tforms_f12_xmyn_keys;
        R12TwoBodyIntKeyCreator tformkey_creator(
          moints_runtime4(),
          GG1space,
          gg1space,
          GG2space,
          gg2space,
          corrfactor(),
          true
          );
        fill_container(tformkey_creator,tforms_f12_xmyn_keys);

      compute_tbint_tensor<ManyBodyTensors::I_to_T,true,false>(
	  V_[s], corrfactor()->tbint_type_f12eri(),
      GG1space, gg1space,
      GG2space, gg2space,
	  antisymmetrize,
	  tforms_f12_xmyn_keys
	  );

      // g12*g12' operator
      std::vector<std::string> tforms_f12f12_xzyw_keys;
      {
      R12TwoBodyIntKeyCreator tformkey_creator(
          moints_runtime4(),
          GG1space,
          GG1space,
          GG2space,
          GG2space,
          corrfactor(),
          true,true
          );
      fill_container(tformkey_creator,tforms_f12f12_xzyw_keys);
      }
      compute_tbint_tensor<ManyBodyTensors::I_to_T,true,true>(
      X_[s], corrfactor()->tbint_type_f12f12(),
      GG1space, GG1space,
      GG2space, GG2space,
      antisymmetrize,
      tforms_f12f12_xzyw_keys
      );
      // 0.5 ( g12*[T,g12'] + [g12,T]*g12' ) = [g12,[t1,g12']] + 0.5 ( g12*[T,g12'] - g12'*[T,g12] )
      //   = [g12,[t1,g12']] - 0.5 ( g12*[g12',T] - g12'*[g12,T] ) = [g12,[t1,g12']] - 0.5 ( (beta-alpha)/(beta+alpha) * [g12*g12',T] ),
      // where the last step valid is for 2 primitive gaussians only (must be contracted otherwise)
      // The first term is symmetric with respect to permutation of g12 and g12' and computed directly,
      compute_tbint_tensor<ManyBodyTensors::I_to_T,true,true>(
      B_[s], corrfactor()->tbint_type_f12t1f12(),
      GG1space, GG1space,
      GG2space, GG2space,
      antisymmetrize,
	  tforms_f12f12_xzyw_keys
      );
      // the second is antisymmetric wrt such permutation and is only needed when number of geminals > 1
      if (corrfactor()->nfunctions() > 1) {

          if (debug_ >= DefaultPrintThresholds::mostO4)
            B_[s].print(prepend_spincase(spincase2,"B(diag;symm)").c_str());
	  RefSCMatrix Banti = B_[s].clone(); Banti.assign(0.0);
	  // the handling of the second term differs between standard approximations {A,A',B} and {C,C'}
	  if (stdapprox() == LinearR12::StdApprox_Ap ||
	      stdapprox() == LinearR12::StdApprox_App ||
	      stdapprox() == LinearR12::StdApprox_B) {
	      // 1) in standard approximations A and B the commutators are explicitly evaluated:
	      //    the second term is exactly what is computed in G12Libint2 under the name t1f12 and t2f12 when g12!=g12'
	       compute_tbint_tensor<ManyBodyTensors::I_to_T,true,true>(
	       Banti, corrfactor()->tbint_type_t1f12(),
	       GG1space, GG1space,
	       GG2space, GG2space,
	       antisymmetrize,
           tforms_f12f12_xzyw_keys
	       );
	       compute_tbint_tensor<ManyBodyTensors::I_to_T,true,true>(
           Banti, corrfactor()->tbint_type_t2f12(),
           GG1space, GG1space,
           GG2space, GG2space,
           antisymmetrize,
           tforms_f12f12_xzyw_keys
           );
	  }
	  else {
	      // 2) in standard approximation C the commutators are evaluated via RI:
	      // Firstly, instead of T we can use h+J (it will be used later anyway)
	      // Second let's designate (beta-alpha)/(beta+alpha) * g12*g12 as A12. A12 is computed as anti_f12f12 by G12NCLibint2
	      // [A12,T] = [A12,h+J] = A12 (hJ_1 + hJ_2) -  (hJ_1 + hJ_2) A12
	      Ref<OrbitalSpace> hj_x1 = hj_x_P(spin1);
	      Ref<OrbitalSpace> hj_x2 = hj_x_P(spin2);

	      // <xy|hJ z> tforms
	      std::vector<std::string> tforms_xyHz_keys;
	      {
		  R12TwoBodyIntKeyCreator tformkey_creator(moints_runtime4(),
		                                           GG1space,hj_x1,GG2space,GG2space,
		                                           corrfactor(),
		                                           true,true);
		  fill_container(tformkey_creator,tforms_xyHz_keys);
	      }
	      // <hJ z|xy> tforms
	      std::vector<std::string> tforms_Hzxy_keys;
	      {
		  R12TwoBodyIntKeyCreator tformkey_creator(moints_runtime4(),
		                                           hj_x1,GG1space,GG2space,GG2space,
		                                           corrfactor(),
		                                           true,true);
		  fill_container(tformkey_creator,tforms_Hzxy_keys);
	      }

	      compute_tbint_tensor<ManyBodyTensors::I_to_T,true,true>(
		  Banti, corrfactor()->tbint_type_f12f12_anti(),
          GG1space, hj_x1,
          GG2space, GG2space,
		  antisymmetrize,
		  tforms_xyHz_keys
		  );
	      compute_tbint_tensor<ManyBodyTensors::I_to_mT,true,true>(
		  Banti, corrfactor()->tbint_type_f12f12_anti(),
          hj_x1, GG1space,
          GG2space, GG2space,
		  antisymmetrize,
		  tforms_Hzxy_keys
		  );

	      if (!x1_eq_x2) {
		  // <xy|z hJ> tforms
		  std::vector<std::string> tforms_xyzH_keys;
		  {
		      R12TwoBodyIntKeyCreator tformkey_creator(moints_runtime4(),
		                                               GG1space,GG1space,GG2space,hj_x2,
		                                               corrfactor(),
		                                               true,true);
		      fill_container(tformkey_creator,tforms_xyzH_keys);
		  }
		  // <z hJ|xy> tforms
		  std::vector<std::string> tforms_zHxy_keys;
		  {
		      R12TwoBodyIntKeyCreator tformkey_creator(moints_runtime4(),
		                                               GG1space,GG1space,hj_x2,GG2space,
		                                               corrfactor(),
		                                               true,true);
		      fill_container(tformkey_creator,tforms_zHxy_keys);
		  }

		  compute_tbint_tensor<ManyBodyTensors::I_to_T,true,true>(
		      Banti, corrfactor()->tbint_type_f12f12_anti(),
		      GG1space, GG1space,
              GG2space, hj_x2,
		      antisymmetrize,
		      tforms_xyzH_keys
		      );
		  compute_tbint_tensor<ManyBodyTensors::I_to_mT,true,true>(
		      Banti, corrfactor()->tbint_type_f12f12_anti(),
              GG1space, GG1space,
              hj_x2, GG2space,
		      antisymmetrize,
		      tforms_zHxy_keys
		      );
	      }
	      else {
                  // contribution from particle 2 is the same as from particle 1
		    Banti.scale(2.0);
            symmetrize<false>(Banti,Banti,GG1space,GG1space);
	      }
	  }
	  Banti.scale(-0.5);
          if (debug_ >= DefaultPrintThresholds::mostO4)
            Banti.print(prepend_spincase(spincase2,"B(diag;antisymm)").c_str());
	  B_[s].accumulate(Banti);
      }
  }

  ExEnv::out0() << decindent;
  ExEnv::out0() << indent << "Exited G12 diagonal intermediates evaluator" << endl;

  tim_diagonal.exit();
  checkpoint_();

  return;
}

////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
