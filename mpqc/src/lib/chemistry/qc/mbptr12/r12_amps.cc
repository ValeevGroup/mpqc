//
// r12_amps.cc
//
// Copyright (C) 2004 Edward Valeev
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

#ifdef __GNUG__
#pragma implementation
#endif

#include <math/scmat/local.h>
#include <chemistry/qc/mbptr12/r12_amps.h>
#include <chemistry/qc/mbptr12/r12int_eval.h>
#include <chemistry/qc/mbptr12/creator.h>
#include <chemistry/qc/mbptr12/container.h>

using namespace std;
using namespace sc;

F12Amplitudes::F12Amplitudes(const Ref<R12IntEval>& r12eval) :
  r12eval_(r12eval)
{
}

F12Amplitudes::~F12Amplitudes()
{
}

const RefSCMatrix&
F12Amplitudes::T2(SpinCase2 S)
{
  if (T2_[S].null())
    compute_(S);
  return T2_[S];
}

const RefSCMatrix&
F12Amplitudes::Fvv(SpinCase2 S)
{
  if (Fvv_[S].null())
    compute_(S);
  return Fvv_[S];
}

const RefSCMatrix&
F12Amplitudes::Foo(SpinCase2 S)
{
  if (Foo_[S].null())
    compute_(S);
  return Foo_[S];
}

const RefSCMatrix&
F12Amplitudes::Fov(SpinCase2 S)
{
  if (Fov_[S].null())
    compute_(S);
  return Fov_[S];
}

const RefSCMatrix&
F12Amplitudes::Fox(SpinCase2 S)
{
  if (Fox_[S].null())
    compute_(S);
  return Fox_[S];
}

const RefSCMatrix&
F12Amplitudes::Fvo(SpinCase2 S)
{
  if (Fvo_[S].null())
    compute_(S);
  return Fvo_[S];
}

const RefSCMatrix&
F12Amplitudes::Fxo(SpinCase2 S)
{
  if (Fxo_[S].null())
    compute_(S);
  return Fxo_[S];
}

void
F12Amplitudes::compute_(SpinCase2 spincase2)
{
  Ref<R12IntEvalInfo> r12info = r12eval_->r12info();
  Ref<SingleRefInfo> refinfo = r12info->refinfo();
  const bool obs_eq_vbs = r12info->basis_vir()->equiv(r12info->basis());
  const bool spin_polarized = refinfo->ref()->spin_polarized();
  
  const unsigned int s = static_cast<unsigned int>(spincase2);
  const SpinCase1 spin1 = case1(spincase2);
  const SpinCase1 spin2 = case2(spincase2);
  const bool p1_neq_p2 = spin_polarized && spincase2 == AlphaBeta;
  
  Ref<MOIndexSpace> occ1_act = r12eval_->occ_act(spin1);
  Ref<MOIndexSpace> occ2_act = r12eval_->occ_act(spin2);
  Ref<MOIndexSpace> occ1 = r12eval_->occ(spin1);
  Ref<MOIndexSpace> occ2 = r12eval_->occ(spin2);
  Ref<MOIndexSpace> ribs1 = r12info->ribs_space(spin1);
  Ref<MOIndexSpace> ribs2 = r12info->ribs_space(spin2);
  Ref<MOIndexSpace> vir1_act = r12eval_->vir_act(spin1);
  Ref<MOIndexSpace> vir2_act = r12eval_->vir_act(spin2);
  
  // Allocate the matrices
  RefSCDimension dim_f12 = r12eval_->dim_f12(spincase2);
  RefSCDimension dim_aa = r12eval_->dim_oo(spincase2);
  if (dim_f12.n() == 0) return;
  RefSCDimension dim_oo = new SCDimension(spincase2 != AlphaBeta ?
                                          occ1->rank()*(occ1->rank()+1)/2 :
                                          occ1->rank() * occ2->rank());
  RefSCDimension dim_vv = new SCDimension(spincase2 != AlphaBeta ?
                                          vir1_act->rank()*(vir1_act->rank()+1)/2 :
                                          vir1_act->rank() * vir2_act->rank());
  RefSCDimension dim_ov = new SCDimension(occ1->rank() * vir2_act->rank());
  RefSCDimension dim_ox = new SCDimension(occ1->rank() * ribs2->rank());
  RefSCDimension dim_vo = new SCDimension(occ2->rank() * vir1_act->rank());
  RefSCDimension dim_xo = new SCDimension(occ2->rank() * ribs1->rank());
  Ref<SCMatrixKit> kit = new LocalSCMatrixKit;
  T2_[s] = kit->matrix(dim_aa,dim_vv);  T2_[s].assign(0.0);
  Fvv_[s] = kit->matrix(dim_f12,dim_vv);  Fvv_[s].assign(0.0);
  Foo_[s] = kit->matrix(dim_f12,dim_oo);  Foo_[s].assign(0.0);
  Fov_[s] = kit->matrix(dim_f12,dim_ov);  Fov_[s].assign(0.0);
  Fox_[s] = kit->matrix(dim_f12,dim_ox);  Fox_[s].assign(0.0);
  if (p1_neq_p2) {
    Fvo_[s] = kit->matrix(dim_f12,dim_vo);  Fvo_[s].assign(0.0);
    Fxo_[s] = kit->matrix(dim_f12,dim_xo);  Fxo_[s].assign(0.0);
  }
  
  Ref<TwoBodyMOIntsTransform> tform0_pp;
  std::vector<  Ref<TwoBodyMOIntsTransform> > tforms_pp;
  std::vector<  Ref<TwoBodyMOIntsTransform> > tforms_mx;
  std::vector<  Ref<TwoBodyMOIntsTransform> > tforms_xm;
  
  // if OBS == VBS then use (ip|ip) and (imjx) integrals
  if (obs_eq_vbs) {
    {
      NewTransformCreator tform_creator(r12eval_,
        occ1_act,
        refinfo->orbs(spin1),
        occ2_act,
        refinfo->orbs(spin2));
        tform0_pp = tform_creator();
    }
    {
      NewTransformCreator tform_creator(r12eval_,
        occ1_act,
        refinfo->orbs(spin1),
        occ2_act,
        refinfo->orbs(spin2),true);
        fill_container(tform_creator,tforms_pp);
    }
    {
      NewTransformCreator tform_creator(r12eval_,
        occ1_act,
        occ1,
        occ2_act,
        ribs2,true);
        fill_container(tform_creator,tforms_mx);
    }
    {
      NewTransformCreator tform_creator(r12eval_,
        occ1_act,
        ribs1,
        occ2_act,
        occ2,true);
        fill_container(tform_creator,tforms_xm);
    }
  }
  // else the needed transforms already exist
  
  const bool antisymm = spincase2!=AlphaBeta;
  r12eval_->compute_T2_(T2_[s],occ1_act,vir1_act,occ2_act,vir2_act,antisymm,tform0_pp);
  r12eval_->compute_F12_(Fvv_[s],occ1_act,vir1_act,occ2_act,vir2_act,antisymm,tforms_pp);
  r12eval_->compute_F12_(Foo_[s],occ1_act,occ1,occ2_act,occ2,antisymm,tforms_pp);
  
  // WARNING cannot antisymmetrize matrices if p1_neq_p2. Should throw but will do nothing for now
  if (!antisymm) {
    r12eval_->compute_F12_(Fov_[s],occ1_act,occ1,occ2_act,vir2_act,antisymm,tforms_pp);
    r12eval_->compute_F12_(Fox_[s],occ1_act,occ1,occ2_act,ribs2,antisymm,tforms_mx);
    if (p1_neq_p2) {
      Fvo_[s] = kit->matrix(dim_f12,dim_vo);
      Fxo_[s] = kit->matrix(dim_f12,dim_xo);
      r12eval_->compute_F12_(Fvo_[s],occ1_act,vir1_act,occ2_act,occ2,antisymm,tforms_pp);
      r12eval_->compute_F12_(Fxo_[s],occ1_act,ribs1,occ2_act,occ2,antisymm,tforms_xm);
    }
  }
  
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
