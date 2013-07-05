//
// FxFgen.cc
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

#include <stdexcept>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#include <mpqc_config.h>
#include <util/misc/formio.h>
#include <util/misc/regtime.h>
#include <util/class/class.h>
#include <util/state/state.h>
#include <util/state/state_text.h>
#include <util/state/state_bin.h>
#include <math/scmat/matrix.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/basis/integral.h>
#include <math/distarray4/distarray4.h>
#include <chemistry/qc/mbptr12/r12wfnworld.h>
#include <math/mmisc/pairiter.h>
#include <chemistry/qc/mbptr12/r12int_eval.h>
#include <chemistry/qc/mbptr12/creator.h>
#include <chemistry/qc/mbptr12/container.h>
#include <chemistry/qc/mbptr12/compute_tbint_tensor.h>
#include <chemistry/qc/mbptr12/contract_tbint_tensor.h>
#include <chemistry/qc/mbptr12/twoparticlecontraction.h>
#include <chemistry/qc/lcao/utils.h>
#include <chemistry/qc/lcao/utils.impl.h>
#include <util/misc/print.h>

using namespace std;
using namespace sc;

void
R12IntEval::compute_FxF_(RefSCMatrix& FxF,
                         SpinCase2 spincase2,
                         const Ref<OrbitalSpace>& bra1,
                         const Ref<OrbitalSpace>& bra2,
                         const Ref<OrbitalSpace>& ket1,
                         const Ref<OrbitalSpace>& ket2,
                         const Ref<OrbitalSpace>& int1,
                         const Ref<OrbitalSpace>& int2,
                         const Ref<OrbitalSpace>& intk1,
                         const Ref<OrbitalSpace>& intk2,
                         const Ref<OrbitalSpace>& intkx1,
                         const Ref<OrbitalSpace>& intkx2
                        )
{
  const bool abs_eq_obs = r12world()->basis()->equiv(r12world()->basis_ri());
  const bool part1_equiv_part2 = (bra1 == bra2 && ket1 == ket2);

  // Check semantics
  bool correct_semantics = (spincase2 != AlphaBeta && part1_equiv_part2) ||
                           (spincase2 == AlphaBeta);
  correct_semantics = correct_semantics &&
                      intk1->rank() == intkx1->rank() &&
                      intk2->rank() == intkx2->rank();
  if (!correct_semantics)
    throw ProgrammingError("R12IntEval::compute_FxF_() -- incorrect call semantics",__FILE__,__LINE__);

#if 0
  // check if summation spaces consistent with part1_equiv_part2
  const bool int_1_equiv_2 = (int1 == int2 && intk1 == intk2 && intkx1 == intkx2);
  if (part1_equiv_part2 ^ int_1_equiv_2)
    throw ProgrammingError("R12IntEval::compute_FxF_() -- contraction spaces must have same permutational symmetry as outer spaces",__FILE__,__LINE__);
#endif

  // heuristic rules to make supporting half-transformed integral types less numerous
  // 1) if p1 and p2 are not equivalent, both types of integrals (x1_2 and 1_x2) will be needed -- nothing can be done here
  // 2) if p1 and p2 are equivalent, choose x1_2 or 1_x2 path so that AO basis rank of p2 is greater than of p1
  const bool reorder = intkx1->basis()->nbasis() > int2->basis()->nbasis();
  const bool do_x1_2 = part1_equiv_part2 ? (!reorder) : true;
  const bool do_1_x2 = part1_equiv_part2 ? (reorder) : true;

  Timer tim("generic FxF intermediate");
  ExEnv::out0() << indent << "Entered generic FxF intermediate evaluator" << endl;
  ExEnv::out0() << incindent;

  const unsigned int nf12 = corrfactor()->nfunctions();
  SpinMOPairIter braiter(bra1->rank(),bra2->rank(),spincase2);
  SpinMOPairIter ketiter(ket1->rank(),ket2->rank(),spincase2);
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

  using mbptr12::TwoParticleContraction;
  using mbptr12::Direct_Contraction;

  // If only x1_2 or 1_x2 is computed (this happens only if p1 is equiv to p2), need to double its weight; symmetrization later will take care of things
  const double perm_pfac = (part1_equiv_part2 ? 2.0 : 1.0);

  if (do_x1_2) {
  // (bra1 intkx1 |bra2 int2) tforms
  std::vector<std::string> tforms_Fbra_k1;
  {
    R12TwoBodyIntKeyCreator tformkey_creator(moints_runtime4(),bra1,intkx1,bra2,int2,
                                             corrfactor(),true);
    fill_container(tformkey_creator,tforms_Fbra_k1);
  }
  // (ket1 intk1 |ket2 int2) tforms
  std::vector<std::string> tforms_Fket_k1;
  {
    R12TwoBodyIntKeyCreator tformkey_creator(moints_runtime4(),ket1,intk1,ket2,int2,
                                             corrfactor(),true);
    fill_container(tformkey_creator,tforms_Fket_k1);
  }
  // contract
  contract_tbint_tensor<true,true>
    (
      FxF, corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_f12(),
      perm_pfac,
      bra1, bra2,
      intkx1, int2,
      ket1, ket2,
      intk1, int2,
      spincase2!=AlphaBeta, tforms_Fbra_k1, tforms_Fket_k1
    );
  }

  if (do_1_x2) {

    // (bra1 intb1 |bra2 intkx2) tforms
    std::vector<std::string> tforms_Fbra_k2;
    {
      R12TwoBodyIntKeyCreator tformkey_creator(moints_runtime4(),bra1,int1,bra2,intkx2,
                                               corrfactor(),true);
      fill_container(tformkey_creator,tforms_Fbra_k2);
    }
    // (ket1 int1 |ket2 intk2) tforms
    std::vector<std::string> tforms_Fket_k2;
    {
      R12TwoBodyIntKeyCreator tformkey_creator(moints_runtime4(),ket1,int1,ket2,intk2,
                                               corrfactor(),true);
      fill_container(tformkey_creator,tforms_Fket_k2);
    }
    // contract
    contract_tbint_tensor<true,true>
      (
        FxF, corrfactor()->tbint_type_f12(), corrfactor()->tbint_type_f12(),
        perm_pfac,
        bra1, bra2,
        int1, intkx2,
        ket1, ket2,
        int1, intk2,
        spincase2!=AlphaBeta, tforms_Fbra_k2, tforms_Fket_k2
      );
  }

  // if particles are equivalent and spincase alpha-beta -- symmetrize just in case
  if (part1_equiv_part2 && spincase2 == AlphaBeta)
    symmetrize<false>(FxF,FxF,bra1,ket1);

  if (debug_ >= DefaultPrintThresholds::allO4) {
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
}

////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
