//
// mp2_pair_energies.cc
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

#include <sstream>

#include <util/misc/timer.h>
#include <chemistry/qc/mbptr12/r12int_eval.h>
#include <chemistry/qc/mbptr12/pairiter.h>

using namespace std;
using namespace sc;

void
R12IntEval::compute_mp2_pair_energies_(SpinCase2 S)
{
  const SpinCase1 spin1 = case1(S);
  const SpinCase1 spin2 = case2(S);

  const Ref<SingleRefInfo> refinfo = r12info()->refinfo();
  Ref<MOIndexSpace> occ1_act = refinfo->occ_act(spin1);
  Ref<MOIndexSpace> occ2_act = refinfo->occ_act(spin2);
  Ref<MOIndexSpace> vir1_act = r12info()->vir_act(spin1);
  Ref<MOIndexSpace> vir2_act = r12info()->vir_act(spin2);
  const int nocc1_act = occ1_act->rank();
  const int nocc2_act = occ2_act->rank();
  const int nvir1_act = vir1_act->rank();
  const int nvir2_act = vir2_act->rank();
  RefDiagSCMatrix evals_occ1_act = occ1_act->evals();
  RefDiagSCMatrix evals_occ2_act = occ2_act->evals();
  RefDiagSCMatrix evals_vir1_act = vir1_act->evals();
  RefDiagSCMatrix evals_vir2_act = vir2_act->evals();

  // Figure out which transform to use -- if VBS != OBS then use (ia|jb), else (ip|jq)
  const bool vbs_neq_obs = (refinfo->uocc(Alpha) != r12info()->vir(Alpha));
  // If closed-shell then should compute all pair energies
  const bool compute_all_spincases = !refinfo->ref()->spin_polarized();

  Ref<MOIndexSpace> xspace, yspace;
  int x_offset, y_offset;
  if (vbs_neq_obs) {
    xspace = r12info()->vir_act(spin1);
    yspace = r12info()->vir_act(spin2);
    x_offset = 0;
    y_offset = 0;
  }
  else {
    xspace = refinfo->orbs(spin1);
    yspace = refinfo->orbs(spin2);
    x_offset = refinfo->occ(spin1)->rank();
    y_offset = refinfo->occ(spin2)->rank();
  }
  const int nx = xspace->rank();
  const std::string tform_name = transform_label(occ1_act,
                                                 xspace,
                                                 occ2_act,
                                                 yspace,
                                                 0);
  
  // get the transform object
  // NOTE needs to become standalone function
  Ref<TwoBodyMOIntsTransform> ixjy_tform;
  try {
    ixjy_tform = get_tform_(tform_name);
  }
  catch (TransformNotFound& ex) {
    Ref<MOIntsTransformFactory> tfactory = r12info()->tfactory();
    tfactory->set_ints_method((MOIntsTransformFactory::StoreMethod)r12info()->ints_method());
    tfactory->set_spaces(occ1_act,xspace,occ2_act,yspace);
    ixjy_tform = tfactory->twobody_transform_13(tform_name);
  }
  
  Ref<R12IntsAcc> ijxy_acc = ixjy_tform->ints_acc();
  if (ijxy_acc.null() || !ijxy_acc->is_committed()) {
    // only need ERIs
    Ref<TwoBodyIntDescr> tbintdescr = new TwoBodyIntDescrERI(r12info()->integral());
    ixjy_tform->compute(tbintdescr);
  }
  // Should make something like this possible:
  //ixjy_tform->compute(corrfactor()->function(0));
  if (!ijxy_acc->is_active())
    ijxy_acc->activate();

  SpinMOPairIter ij_iter(occ1_act,occ2_act,S);
  SpinMOPairIter xy_iter(xspace,yspace,S);

  vector<int> proc_with_ints;
  const int nproc_with_ints = tasks_with_ints_(ijxy_acc,proc_with_ints);
  const int me = r12info()->msg()->me();

  if (ijxy_acc->has_access(me)) {
    RefSCVector emp2pair = emp2pair_[S];
    for(ij_iter.start(); ij_iter; ij_iter.next()) {
      const int ij = ij_iter.ij();

      const int ij_proc = ij%nproc_with_ints;
      if (ij_proc != proc_with_ints[me])
        continue;

      const int i = ij_iter.i();
      const int j = ij_iter.j();

      if (debug_)
          ExEnv::outn() << indent << "task " << me << ": working on (i,j) = " << i << "," << j << " " << endl;
      tim_enter("MO ints retrieve");
      const double *ijxy_buf_eri = ijxy_acc->retrieve_pair_block(i,j,corrfactor()->tbint_type_eri());
      tim_exit("MO ints retrieve");
      if (debug_)
        ExEnv::outn() << indent << "task " << me << ": obtained ij blocks" << endl;

      double emp2 = 0.0;
      double emp2_aa = 0.0;
      double emp2_ab = 0.0;
      for(xy_iter.start(); xy_iter; xy_iter.next()) {
        const int x = xy_iter.i();
        const int y = xy_iter.j();
        const int a = x - x_offset;
        if (a<0)
          continue;
        const int b = y - y_offset;
        if (b<0)
          continue;
        const int xy = x*nx+y;
        const int yx = y*nx+x;

        const double ERI_xy = ijxy_buf_eri[xy];
        const double ERI_yx = ijxy_buf_eri[yx];
        const double denom = 1.0/(evals_occ1_act(i) + evals_occ2_act(j) - evals_vir1_act(a) - evals_vir2_act(b));
        
        if (debug_ > 2)
          ExEnv::out0() << "i = " << i << " j = " << j << " a = " << x << " b = " << y
                        << " <ij|ab> = " << ERI_xy << " <ij|ba> = " << ERI_yx
                        << " denom = " << denom << endl;

        if (compute_all_spincases) {
          const double ERI_aa = ERI_xy - ERI_yx;
          emp2_aa += 0.5*ERI_aa*ERI_aa*denom;
          emp2_ab += ERI_xy*ERI_xy*denom;
        }
        else {
          if (S == AlphaBeta)
            emp2 += ERI_xy * ERI_xy * denom;
          else
            emp2 += (ERI_xy-ERI_yx) * (ERI_xy-ERI_yx) * denom;
        }
      }
      ijxy_acc->release_pair_block(i,j,corrfactor()->tbint_type_eri());

      if (compute_all_spincases) {
        emp2pair_[AlphaBeta].set_element(ij,emp2_ab);
        if (i != j) {
          const int ij_aa = i*(i-1)/2 + j;
          emp2pair_[AlphaAlpha].set_element(ij_aa,emp2_aa);
          emp2pair_[BetaBeta].set_element(ij_aa,emp2_aa);
        }
      }
      else
        emp2pair.set_element(ij,emp2);
    }
  }
}


void
R12IntEval::compute_mp2_pair_energies_(RefSCVector& emp2pair,
                                       SpinCase2 S,
                                       const Ref<MOIndexSpace>& space1,
                                       const Ref<MOIndexSpace>& space2,
                                       const Ref<MOIndexSpace>& space3,
                                       const Ref<MOIndexSpace>& space4,
                                       const Ref<TwoBodyMOIntsTransform>& transform)
{
  // Check correct semantics of this call : if S != AlphaBeta then space1==space3 and space2==space4
  const bool correct_semantics = (S != AlphaBeta && space1==space3 && space2==space4) ||
                                 S==AlphaBeta;
  if (!correct_semantics)
    throw ProgrammingError("R12IntEval::compute_mp2_pair_energies_() -- incorrect call semantics",
                           __FILE__,__LINE__);
  emp2pair.assign(0.0);

  //
  // If transform not given, construct appropriately, otherwise notify user which transform is used
  //
  Ref<TwoBodyMOIntsTransform> tform = transform;
  if (tform.null()) {
    // Only need 1/r12 integrals, hence doesn't matter which f12 to use
    tform = get_tform_(transform_label(space1,space2,space3,space4,0));
  }
  else {
    if (debug_ > 0)
      ExEnv::out0() << "Using transform " << tform->name() << std::endl;
  }
  
  
  //
  // Initialize spaces and maps
  //
  // maps spaceX to spaceX of the transform
  std::vector<unsigned int> map1, map2, map3, map4;
  // maps space4 to space2 of transform
  std::vector<unsigned int> map24;
  // maps space2 to space4 of transform
  std::vector<unsigned int> map42;
  Ref<MOIndexSpace> tspace1 = tform->space1();
  Ref<MOIndexSpace> tspace2 = tform->space2();
  Ref<MOIndexSpace> tspace3 = tform->space3();
  Ref<MOIndexSpace> tspace4 = tform->space4();
  map1 = *tspace1<<*space1;
  map2 = *tspace2<<*space2;
  map3 = *tspace3<<*space3;
  map4 = *tspace4<<*space4;
  if (S != AlphaBeta) {
    if (tspace2 == tspace4) {
      map24 = map2;
      map42 = map4;
    }
    else {
      map24 = *tspace2<<*space4;
      map42 = *tspace4<<*space2;
    }
  }
  
  const unsigned int rank1 = space1->rank();
  const unsigned int rank2 = space2->rank();
  const unsigned int rank3 = space3->rank();
  const unsigned int rank4 = space4->rank();
  const unsigned int trank4 = tspace4->rank();
  const RefDiagSCMatrix evals1 = space1->evals();
  const RefDiagSCMatrix evals2 = space2->evals();
  const RefDiagSCMatrix evals3 = space3->evals();
  const RefDiagSCMatrix evals4 = space4->evals();
  
  // Using spinorbital iterators means I don't take into account perm symmetry
  // More efficient algorithm will require generic code
  SpinMOPairIter iter13(space1,space3,S);
  SpinMOPairIter iter24(space2,space4,S);
  
  Ref<R12IntsAcc> accum = tform->ints_acc();
  if (accum.null() || !accum->is_committed()) {
    // only need ERIs
    Ref<TwoBodyIntDescr> tbintdescr = new TwoBodyIntDescrERI(r12info()->integral());
    tform->compute(tbintdescr);
  }
  accum->activate();
  
  tim_enter("MP2 pair energies");
  std::ostringstream oss;
  oss << "<" << space1->id() << " " << space3->id() << "|"
      << space2->id() << " " << space4->id() << ">";
  const std::string label = oss.str();
  ExEnv::out0() << endl << indent
	       << "Entered MP2 pair energies (" << label << ") evaluator" << std::endl;
  ExEnv::out0() << incindent;
  
  vector<int> proc_with_ints;
  const int nproc_with_ints = tasks_with_ints_(accum,proc_with_ints);
  const int me = r12info()->msg()->me();
  
  if (accum->has_access(me)) {
    for(iter13.start(); iter13; iter13.next()) {
      const int ij = iter13.ij();
      
      const int ij_proc = ij%nproc_with_ints;
      if (ij_proc != proc_with_ints[me])
        continue;
      
      const unsigned int i = iter13.i();
      const unsigned int j = iter13.j();
      const unsigned int ii = map1[i];
      const unsigned int jj = map3[j];
      
      if (debug_)
        ExEnv::outn() << indent << "task " << me << ": working on (i,j) = " << i << "," << j << " " << endl;
      tim_enter("MO ints retrieve");
      const double *ij_buf_eri = accum->retrieve_pair_block(ii,jj,corrfactor()->tbint_type_eri());
      tim_exit("MO ints retrieve");
      if (debug_)
        ExEnv::outn() << indent << "task " << me << ": obtained ij blocks" << endl;
      
      double emp2 = 0.0;
      for(iter24.start(); iter24; iter24.next()) {
        const unsigned int a = iter24.i();
        const unsigned int b = iter24.j();
        const unsigned int aa = map2[a];
        const unsigned int bb = map4[b];
        const int ab = aa*trank4+bb;
        
        const double ERI_iajb = ij_buf_eri[ab];
        const double denom = 1.0/(evals1(i) + evals3(j) - evals2(a) - evals4(b));
        
        if (S == AlphaBeta) {
          emp2 += ERI_iajb*ERI_iajb*denom;
          if (debug_ > 2) {
            ExEnv::out0() << "i = " << i << " j = " << j << " a = " << a << " b = " << b
            << " <ij|ab> = " << ERI_iajb
            << " denom = " << denom << endl;
          }
        }
        else {
          const int aa = map42[a];
          const int bb = map24[b];
          const int ba = bb*trank4+aa;
          const double ERI_ibja = ij_buf_eri[ba];
          const double ERI_anti = ERI_iajb - ERI_ibja;
          emp2 += ERI_anti*ERI_anti*denom;
          if (debug_ > 2) {
            ExEnv::out0() << "i = " << i << " j = " << j << " a = " << a << " b = " << b
            << " <ij|ab> = " << ERI_iajb
            << " <ij|ba> = " << ERI_ibja
            << " denom = " << denom << endl;
          }
        }
        
      }
      accum->release_pair_block(ii,jj,corrfactor()->tbint_type_eri());
      emp2pair.set_element(ij,emp2);
    }
  }
  
  accum->deactivate();
  ExEnv::out0() << decindent;
  ExEnv::out0() << indent << "Exited MP2 pair energies (" << label << ") evaluator" << endl;
  tim_exit("MP2 pair energies");
}

