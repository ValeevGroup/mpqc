//
// vxb_eval_abs_a.cc
//
// Copyright (C) 2003 Edward Valeev
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

#include <util/misc/formio.h>
#include <util/ref/ref.h>
#include <chemistry/qc/mbptr12/vxb_eval_info.h>
#include <chemistry/qc/mbptr12/vxb_eval_abs_a.h>

using namespace std;
using namespace sc;

inline int max(int a,int b) { return (a > b) ? a : b;}

/*-----------------
  R12IntEval_abs_A
 -----------------*/
static ClassDesc R12IntEval_abs_A_cd(
  typeid(R12IntEval_abs_A),"R12IntEval_abs_A",1,"virtual public SavableState",
  0, 0, create<R12IntEval_abs_A>);

R12IntEval_abs_A::R12IntEval_abs_A(Ref<R12IntEvalInfo>& r12info) :
  r12info_(r12info), r12intsacc_(NULL)
{
  evaluated_ = false;
}

R12IntEval_abs_A::R12IntEval_abs_A(StateIn& si) : SavableState(si)
{
  r12info_ << SavableState::restore_state(si);
  r12intsacc_ << SavableState::restore_state(si);

  int evaluated; si.get(evaluated); evaluated_ = (bool) evaluated;
}

R12IntEval_abs_A::~R12IntEval_abs_A()
{
  r12info_ = 0;
  r12intsacc_ = 0;
}

void R12IntEval_abs_A::save_data_state(StateOut& so)
{
  SavableState::save_state(r12info_.pointer(),so);

  so.put((int)evaluated_);
}

void R12IntEval_abs_A::obsolete()
{
  evaluated_ = false;
  r12intsacc_ = 0;
}

Ref<R12IntEvalInfo> R12IntEval_abs_A::r12info() const { return r12info_; };

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
