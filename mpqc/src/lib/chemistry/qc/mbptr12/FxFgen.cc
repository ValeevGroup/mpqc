//
// FxFgen.cc
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

#include <stdexcept>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#include <scconfig.h>
#include <util/misc/formio.h>
#include <util/misc/timer.h>
#include <util/class/class.h>
#include <util/state/state.h>
#include <util/state/state_text.h>
#include <util/state/state_bin.h>
#include <math/scmat/matrix.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/mbptr12/blas.h>
#include <chemistry/qc/mbptr12/r12ia.h>
#include <chemistry/qc/mbptr12/vxb_eval_info.h>
#include <chemistry/qc/mbptr12/pairiter.h>
#include <chemistry/qc/mbptr12/r12int_eval.h>
#include <chemistry/qc/mbptr12/creator.h>
#include <chemistry/qc/mbptr12/container.h>
#include <chemistry/qc/mbptr12/compute_tbint_tensor.h>
#include <chemistry/qc/mbptr12/contract_tbint_tensor.h>
#include <chemistry/qc/mbptr12/twoparticlecontraction.h>
#include <chemistry/qc/mbptr12/utils.h>
#include <chemistry/qc/mbptr12/utils.impl.h>

using namespace std;
using namespace sc;

void
R12IntEval::compute_FxF_(RefSCMatrix& FxF,
                         SpinCase2 spincase2,
                         const Ref<MOIndexSpace>& bra1,
                         const Ref<MOIndexSpace>& bra2,
                         const Ref<MOIndexSpace>& ket1,
                         const Ref<MOIndexSpace>& ket2,
                         const Ref<MOIndexSpace>& int1,
                         const Ref<MOIndexSpace>& int2,
                         const Ref<MOIndexSpace>& intk1,
                         const Ref<MOIndexSpace>& intk2,
                         const Ref<MOIndexSpace>& intkx1,
                         const Ref<MOIndexSpace>& intkx2
                        )
{
  const bool abs_eq_obs = r12info()->basis()->equiv(r12info()->basis_ri());
  const bool part1_equiv_part2 = (bra1 == bra2 && ket1 == ket2);
  
  // Check semantics
  bool correct_semantics = (spincase2 != AlphaBeta && part1_equiv_part2) ||
                           (spincase2 == AlphaBeta);
  correct_semantics = correct_semantics &&
                      intk1->rank() == intkx1->rank() &&
                      intk2->rank() == intkx2->rank();
  if (!correct_semantics)
    throw ProgrammingError("R12IntEval::compute_FxF_() -- incorrect call semantics",__FILE__,__LINE__);
  
  // check if summation spaces consistent with part1_equiv_part2
  const bool int_1_equiv_2 = (int1 == int2 && intk1 == intk2 && intkx1 == intkx2);
  if (part1_equiv_part2 ^ int_1_equiv_2)
    throw ProgrammingError("R12IntEval::compute_FxF_() -- contraction spaces must have same permutational symmetry as outer spaces",__FILE__,__LINE__);
  
  tim_enter("generic FxF intermediate");
  ExEnv::out0() << indent << "Entered generic FxF intermediate evaluator" << endl;
  ExEnv::out0() << incindent;
  
  const unsigned int nf12 = corrfactor()->nfunctions();
  SpinMOPairIter braiter(bra1,bra2,spincase2);
  SpinMOPairIter ketiter(ket1,ket2,spincase2);
  const unsigned int nbra = nf12 * braiter.nij();
  const unsigned int nket = nf12 * ketiter.nij();

  if (FxF.null()) {
    // use the same matrix kit as the intermediates
    FxF = B_[AlphaBeta].kit()->matrix(new SCDimension(nbra),
                                      new SCDimension(nket));
    FxF.assign(0.0);
  }
  else {
    if (FxF.rowdim().n() != nbra)
      throw ProgrammingError("R12IntEval::compute_FxF_() -- row dimension of the given FxF doesn't match given bra dimensions",__FILE__,__LINE__);
    if (FxF.coldim().n() != nket)
      throw ProgrammingError("R12IntEval::compute_FxF_() -- column dimension of the given FxF doesn't match given ket dimensions",__FILE__,__LINE__);
  }
  
  const SpinCase1 spin1 = case1(spincase2);
  const SpinCase1 spin2 = case2(spincase2);
  
#if 0
  Ref<SingleRefInfo> refinfo = r12info()->refinfo();
  Ref<MOIndexSpace> occ1 = refinfo->occ(spin1);
  Ref<MOIndexSpace> occ2 = refinfo->occ(spin2);
  Ref<MOIndexSpace> vir1 = vir(spin1);
  Ref<MOIndexSpace> vir2 = vir(spin2);
  Ref<MOIndexSpace> orbs1 = refinfo->orbs(spin1);
  Ref<MOIndexSpace> orbs2 = refinfo->orbs(spin2);
  // if orbs1 and orbs2 have different rank -- something is TERRIBLY wrong
  if (orbs1->rank() != orbs2->rank())
    throw ProgrammingError("R12IntEval::compute_FxF_() -- orbs1 and orbs2 have different ranks",__FILE__,__LINE__);
  const unsigned int nobs = orbs1->rank();
#endif
  
  Ref<R12IntEval> thisref(this);
  
  using LinearR12::TwoParticleContraction;
  using LinearR12::Direct_Contraction;
  // If particles are not equivalent, will add the 21 contribution too
  const double perm_pfac = (part1_equiv_part2 ? 2.0 : 1.0);
  Ref<TwoParticleContraction> dircontract_k1 =
    new Direct_Contraction(intk1->rank(),int2->rank(),perm_pfac);
  // (bra1 intkx1 |bra2 int2) tforms
  std::vector<  Ref<TwoBodyMOIntsTransform> > tforms_Fbra_k1;
  {
    NewTransformCreator tform_creator(thisref,bra1,intkx1,bra2,int2,true);
    fill_container(tform_creator,tforms_Fbra_k1);
  }
  // (ket1 intk1 |ket2 int2) tforms
  std::vector<  Ref<TwoBodyMOIntsTransform> > tforms_Fket_k1;
  {
    NewTransformCreator tform_creator(thisref,ket1,intk1,ket2,int2,true);
    fill_container(tform_creator,tforms_Fket_k1);
  }
  // contract
  contract_tbint_tensor<
    ManyBodyTensors::I_to_T,
    ManyBodyTensors::I_to_T,
    ManyBodyTensors::I_to_T,
    true,true,false>
    (
      FxF, corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_f12(),
      bra1, bra2,
      intkx1, int2,
      ket1, ket2,
      intk1, int2,
      dircontract_k1,
      spincase2!=AlphaBeta, tforms_Fbra_k1, tforms_Fket_k1
    );
  
  // Do the same for particle 2
  if (!part1_equiv_part2) {
    
    Ref<TwoParticleContraction> dircontract_k2 =
      new Direct_Contraction(int1->rank(),intk2->rank(),1.0);
    // (bra1 intb1 |bra2 intkx2) tforms
    std::vector<  Ref<TwoBodyMOIntsTransform> > tforms_Fbra_k2;
    {
      NewTransformCreator tform_creator(thisref,bra1,int1,bra2,intkx2,true);
      fill_container(tform_creator,tforms_Fbra_k2);
    }
    // (ket1 int1 |ket2 intk2) tforms
    std::vector<  Ref<TwoBodyMOIntsTransform> > tforms_Fket_k2;
    {
      NewTransformCreator tform_creator(thisref,ket1,int1,ket2,intk2,true);
      fill_container(tform_creator,tforms_Fket_k2);
    }
    // contract
    contract_tbint_tensor<
      ManyBodyTensors::I_to_T,
      ManyBodyTensors::I_to_T,
      ManyBodyTensors::I_to_T,
      true,true,false>
      (
        FxF, corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_f12(),
        bra1, bra2,
        int1, intkx2,
        ket1, ket2,
        int1, intk2,
        dircontract_k2,
        spincase2!=AlphaBeta, tforms_Fbra_k2, tforms_Fket_k2
      );
  }
  // particles equivalent -- just symmetrize
  else {
    symmetrize<false>(FxF,FxF,bra1,ket1);
  }
  
  if (debug_ > 1) {
    std::string label = prepend_spincase(spincase2,"generic FxF");
    FxF.print(label.c_str());
  }
  
  // Bra-Ket symmetrize
  FxF.scale(0.5);
  RefSCMatrix FxF_t = FxF.t();
  FxF.accumulate(FxF_t);
  
  globally_sum_scmatrix_(FxF);
  
  ExEnv::out0() << decindent;
  ExEnv::out0() << indent << "Exited generic FxF intermediate evaluator" << endl;

  tim_exit("generic FxF intermediate");
}

////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
