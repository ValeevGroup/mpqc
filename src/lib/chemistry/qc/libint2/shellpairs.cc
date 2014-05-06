//
// shellpairs.cc
//
// Copyright (C) 2001 Edward Valeev
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

#include <util/state/statein.h>
#include <util/state/stateout.h>
#include <util/misc/formio.h>
#include <util/ref/ref.h>
#include <chemistry/qc/basis/gaussshell.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/libint2/shellpairs.h>
#include <chemistry/qc/libint2/primpairs.h>

using namespace std;
using namespace sc;

inline int max(int a,int b) { return (a > b) ? a : b;}
inline void fail()
{
  ExEnv::errn() << scprintf("failing module:\n%s",__FILE__) << endl;
  abort();
}

/*----------------
  ShellPairsLibint2
  ----------------*/
static ClassDesc ShellPairsLibint2_cd(
  typeid(ShellPairsLibint2),"ShellPairsLibint2",1,"virtual public SavableState",
  0, 0, create<ShellPairsLibint2>);

ShellPairsLibint2::ShellPairsLibint2(const Ref<GaussianBasisSet>& bs1,
                                     const Ref<GaussianBasisSet>& bs2) :
    bs1_(bs1), bs2_(bs2), prim_pairs_(new PrimPairsLibint2(bs1_,bs2_)),
    shell_pair_(new ShellPairLibint2(*prim_pairs_))
{
}

ShellPairsLibint2::ShellPairsLibint2(const ShellPairsLibint2& other) :
    bs1_(other.bs1_), bs2_(other.bs2_), prim_pairs_(other.prim_pairs_),
    shell_pair_(new ShellPairLibint2(*prim_pairs_))
{
}

ShellPairsLibint2::ShellPairsLibint2(StateIn& si) :
  SavableState(si)
{
  bs1_ << SavableState::restore_state(si);
  bs2_ << SavableState::restore_state(si);
  prim_pairs_ = new PrimPairsLibint2(bs1_,bs2_);
  shell_pair_ = new ShellPairLibint2(*prim_pairs_);
}

ShellPairsLibint2::~ShellPairsLibint2()
{
  delete shell_pair_;
}

void ShellPairsLibint2::save_data_state(StateOut& so)
{
  SavableState::save_state(bs1_.pointer(),so);
  SavableState::save_state(bs2_.pointer(),so);
}


/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
