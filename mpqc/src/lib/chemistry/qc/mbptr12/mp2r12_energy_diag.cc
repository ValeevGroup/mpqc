//
// mp2r12_energy_diag.cc
//
// Copyright (C) 2010 Jinmei Zhang
//
// Author: Jinmei Zhang <jmzhang@vt.edu>
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

#if 0
#ifdef __GNUG__
#pragma implementation
#endif

#include <chemistry/qc/mbptr12/mp2r12_energy.h>
#include <chemistry/qc/mbptr12/print.h>

using namespace std;
using namespace sc;
#endif

/**********************************
 * class MP2R12Energy_Diag *
 **********************************/

void MP2R12Energy_Diag::compute_ef12() {
  if (evaluated_)
    return;

  Ref<R12WavefunctionWorld> r12world = r12eval_->r12world();
  Ref<MessageGrp> msg = r12world->world()->msg();
  int me = msg->me();
  int ntasks = msg->n();

  const bool obs_eq_ribs = r12world->obs_eq_ribs();
  const bool cabs_empty = obs_eq_ribs;
  // Only diagonal ansatz is supported
  const bool diag = r12world->r12tech()->ansatz()->diag();
  if (diag == false)
    throw ProgrammingError("only diagonal ansatz supported",__FILE__,__LINE__);

  // obtain some preliminaries
  Ref<TwoBodyFourCenterMOIntsRuntime> moints4_rtime = r12world->world()->moints_runtime4();
  Ref<TwoBodyIntDescr> descr_f12 = r12world->r12tech()->corrfactor()->tbintdescr(r12world->integral(),0);
  Ref<TwoBodyIntDescr> descr_f12f12 = r12world->r12tech()->corrfactor()->tbintdescr(r12world->integral(),0,0);
  const std::string descr_f12_key = moints4_rtime->descr_key(descr_f12);
  const std::string descr_f12f12_key = moints4_rtime->descr_key(descr_f12f12);

  //
  // Evaluate pair energies:
  // distribute workload among nodes by pair index
  //
  const int num_unique_spincases2 = (r12eval()->spin_polarized() ? 3 : 2);
  for (int spin = 0; spin < num_unique_spincases2; spin++) {

    const SpinCase2 spincase = static_cast<SpinCase2> (spin);
    if (r12eval()->dim_oo(spincase).n() == 0)
      continue;
    const SpinCase1 spin1 = case1(spincase);
    const SpinCase1 spin2 = case2(spincase);

    const Ref<OrbitalSpace>& occ1_act = r12eval()->occ_act(spin1);
    const Ref<OrbitalSpace>& vir1_act = r12eval()->vir_act(spin1);
    const Ref<OrbitalSpace>& xspace1 = r12eval()->xspace(spin1);
    const Ref<OrbitalSpace>& occ2_act = r12eval()->occ_act(spin2);
    const Ref<OrbitalSpace>& vir2_act = r12eval()->vir_act(spin2);
    const Ref<OrbitalSpace>& xspace2 = r12eval()->xspace(spin2);
    int nocc1_act = occ1_act->rank();
    int nvir1_act = vir1_act->rank();
    int nx1 = xspace1->rank();
    int nocc2_act = occ2_act->rank();
    int nvir2_act = vir2_act->rank();
    int nx2 = xspace2->rank();

    // compute intermediates V, X, B
    //
    // V = (f12/r12) -
    const std::string iiii_key = ParsedTwoBodyFourCenterIntKey::key(occ1_act->id(), occ2_act->id(),
                                                                    occ1_act->id(), occ2_act->id(),
                                                                    descr_f12_key,
                                                                    TwoBodyIntLayout::b1b2_k1k2);
    Ref<TwoBodyMOIntsTransform> iiii_tform = moints4_rtime->get(iiii_key);
    iiii_tform->compute();
    Ref<DistArray4> iiii_ints = iiii_tform->ints_acc();
    iiii_ints->activate();
#if 0
    for(int p=0; p<n; ++p) {
      for(int q=0; q<n; ++q) {
        const double* rs_blk = iiii_ints->retrieve_pair_block(p, q, r12world->r12tech()->corrfactor()->tbint_type_f12eri());

        iiii_ints->release_pair_block(p, q, TwoBodyOper::eri);
      }
    }
#endif

    // for each ij pair compute its contribution to the Hylleraas functional for second-order R12 energy

  } // end of spincase loop

  // Set beta-beta energies to alpha-alpha for closed-shell
  if (!r12world->ref()->spin_polarized()) {
    C_[BetaBeta] = C_[AlphaAlpha];
    emp2f12_[BetaBeta] = emp2f12_[AlphaAlpha];
    ef12_[BetaBeta] = ef12_[AlphaAlpha];
  }

  evaluated_ = true;

  return;
}

