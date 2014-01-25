//
// r12_amps.cc
//
// Copyright (C) 2004 Edward Valeev
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

#include <math/scmat/local.h>
#include <chemistry/qc/mbptr12/r12_amps.h>
#include <chemistry/qc/mbptr12/r12int_eval.h>
#include <chemistry/qc/mbptr12/creator.h>
#include <chemistry/qc/mbptr12/container.h>

using namespace std;
using namespace sc;

R12Amplitudes::R12Amplitudes(const Ref<R12IntEval>& r12eval) :
  r12eval_(r12eval)
{
}

R12Amplitudes::~R12Amplitudes()
{
}

const RefSCMatrix&
R12Amplitudes::T2(SpinCase2 S)
{
  if (T2_[S] == 0)
    compute_(S);
  return T2_[S];
}

const RefSCMatrix&
R12Amplitudes::Fvv(SpinCase2 S)
{
  if (Fvv_[S] == 0)
    compute_(S);
  return Fvv_[S];
}

const RefSCMatrix&
R12Amplitudes::Foo(SpinCase2 S)
{
  if (Foo_[S] == 0)
    compute_(S);
  return Foo_[S];
}

const RefSCMatrix&
R12Amplitudes::Fov(SpinCase2 S)
{
  if (Fov_[S] == 0)
    compute_(S);
  return Fov_[S];
}

const RefSCMatrix&
R12Amplitudes::Fox(SpinCase2 S)
{
  if (Fox_[S] == 0)
    compute_(S);
  return Fox_[S];
}

const RefSCMatrix&
R12Amplitudes::Fvo(SpinCase2 S)
{
  if (Fvo_[S] == 0)
    compute_(S);
  return Fvo_[S];
}

const RefSCMatrix&
R12Amplitudes::Fxo(SpinCase2 S)
{
  if (Fxo_[S] == 0)
    compute_(S);
  return Fxo_[S];
}

void
R12Amplitudes::compute_(SpinCase2 spincase2)
{
  Ref<R12WavefunctionWorld> r12world = r12eval_->r12world();
  Ref<RefWavefunction> ref = r12world->refwfn();
  const bool obs_eq_vbs = r12world->obs_eq_vbs();
  const bool spin_polarized = ref->spin_polarized();

  const unsigned int s = static_cast<unsigned int>(spincase2);
  const SpinCase1 spin1 = case1(spincase2);
  const SpinCase1 spin2 = case2(spincase2);
  const bool p1_neq_p2 = spin_polarized && spincase2 == AlphaBeta;

  Ref<OrbitalSpace> occ1_act = r12eval_->occ_act(spin1);
  Ref<OrbitalSpace> occ2_act = r12eval_->occ_act(spin2);
  Ref<OrbitalSpace> occ1 = r12eval_->occ(spin1);
  Ref<OrbitalSpace> occ2 = r12eval_->occ(spin2);
  Ref<OrbitalSpace> cabs1 = r12world->cabs_space(spin1);
  Ref<OrbitalSpace> cabs2 = r12world->cabs_space(spin2);
  Ref<OrbitalSpace> vir1_act = r12eval_->vir_act(spin1);
  Ref<OrbitalSpace> vir2_act = r12eval_->vir_act(spin2);
  Ref<OrbitalSpace> GGspace1 = r12eval_->GGspace(spin1);
  Ref<OrbitalSpace> GGspace2 = r12eval_->GGspace(spin2);

  // Allocate the matrices
  RefSCDimension dim_f12 = r12eval_->dim_f12(spincase2);
  RefSCDimension dim_aa = r12eval_->dim_oo(spincase2);
  RefSCDimension dim_oo = new SCDimension(spincase2 != AlphaBeta ?
                                          occ1->rank()*(occ1->rank()-1)/2 :
                                          occ1->rank() * occ2->rank());
  RefSCDimension dim_vv = r12eval_->dim_vv(spincase2);
  RefSCDimension dim_ov = new SCDimension(occ1->rank() * vir2_act->rank());
  RefSCDimension dim_ox = new SCDimension(occ1->rank() * cabs2->rank());
  RefSCDimension dim_vo = new SCDimension(occ2->rank() * vir1_act->rank());
  RefSCDimension dim_xo = new SCDimension(occ2->rank() * cabs1->rank());
  Ref<SCMatrixKit> kit = new LocalSCMatrixKit;
  T2_[s] = kit->matrix(dim_aa,dim_vv);  T2_[s].assign(0.0);
  Fvv_[s] = kit->matrix(dim_f12,dim_vv);  Fvv_[s].assign(0.0);
  Foo_[s] = kit->matrix(dim_f12,dim_oo);  Foo_[s].assign(0.0);
  Fov_[s] = kit->matrix(dim_f12,dim_ov);  Fov_[s].assign(0.0);
  Fox_[s] = kit->matrix(dim_f12,dim_ox);  Fox_[s].assign(0.0);
  if (true) {
    Fvo_[s] = kit->matrix(dim_f12,dim_vo);  Fvo_[s].assign(0.0);
    Fxo_[s] = kit->matrix(dim_f12,dim_xo);  Fxo_[s].assign(0.0);
  }
  // If no active orbital pairs for this spin case -- leave
  if (dim_f12.n() == 0) return;

  std::string tform0_pp_key;
  std::vector<std::string> tform_pp_keys;
  std::vector<std::string> tform_mx_keys;
  std::vector<std::string> tform_xm_keys;

  // if OBS == VBS then use (ip|ip) and (imjx) integrals
  if (obs_eq_vbs) {
    {
      tform0_pp_key = ParsedTwoBodyFourCenterIntKey::key(occ1_act->id(),occ2_act->id(),
                                                         ref->orbs(spin1)->id(),ref->orbs(spin2)->id(),
                                                         std::string("ERI"),
                                                         std::string(TwoBodyIntLayout::b1b2_k1k2));
    }
    {
      R12TwoBodyIntKeyCreator tformkey_creator(r12world->world()->moints_runtime4(),
        GGspace1,
        ref->orbs(spin1),
        GGspace2,
        ref->orbs(spin2),
        r12world->r12tech()->corrfactor(),
        true);
        fill_container(tformkey_creator,tform_pp_keys);
    }
    {
      R12TwoBodyIntKeyCreator tform_creator(r12world->world()->moints_runtime4(),
        GGspace1,
        occ1,
        GGspace2,
        cabs2,
        r12world->r12tech()->corrfactor(),
        true);
        fill_container(tform_creator,tform_mx_keys);
    }
    {
      R12TwoBodyIntKeyCreator tformkey_creator(r12world->world()->moints_runtime4(),
        GGspace1,
        cabs1,
        GGspace2,
        occ2,
        r12world->r12tech()->corrfactor(),
        true);
        fill_container(tformkey_creator,tform_xm_keys);
    }
  }
  // else the needed transforms already exist

  const bool antisymm = spincase2!=AlphaBeta;
  r12eval_->compute_T2_(T2_[s],occ1_act,vir1_act,occ2_act,vir2_act,antisymm,tform0_pp_key);
  r12eval_->compute_F12_(Fvv_[s],GGspace1,vir1_act,GGspace2,vir2_act,antisymm,tform_pp_keys);
  r12eval_->compute_F12_(Foo_[s],GGspace1,occ1,GGspace2,occ2,antisymm,tform_pp_keys);

  // WARNING cannot antisymmetrize matrices if p1_neq_p2. Should throw but will do nothing for now
  if (!antisymm) {
    r12eval_->compute_F12_(Fov_[s],GGspace1,occ1,GGspace2,vir2_act,antisymm,tform_pp_keys);
    r12eval_->compute_F12_(Fox_[s],GGspace1,occ1,GGspace2,cabs2,antisymm,tform_mx_keys);
    if (true) {
      Fvo_[s] = kit->matrix(dim_f12,dim_vo);
      Fxo_[s] = kit->matrix(dim_f12,dim_xo);
      r12eval_->compute_F12_(Fvo_[s],GGspace1,vir1_act,GGspace2,occ2,antisymm,tform_pp_keys);
      r12eval_->compute_F12_(Fxo_[s],GGspace1,cabs1,GGspace2,occ2,antisymm,tform_xm_keys);
    }
  }

}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
