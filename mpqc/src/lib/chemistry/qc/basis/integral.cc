//
// integral.cc --- implementation of the Integral class
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Edward Seidl <seidl@janed.com>
// Maintainer: LPS
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

#ifdef __GNUC__
#pragma implementation
#endif

#include <util/state/stateio.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/shellrot.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/basis/obint.h>

SavableState_REF_def(Integral);

#define CLASSNAME Integral
#define PARENTS public SavableState
#include <util/state/statei.h>
#include <util/class/classia.h>
void *
Integral::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

Integral::Integral(const RefGaussianBasisSet &b1,
                   const RefGaussianBasisSet &b2,
                   const RefGaussianBasisSet &b3,
                   const RefGaussianBasisSet &b4)
{
  storage_ = 0;
  storage_used_ = 0;
  grp_ = MessageGrp::get_default_messagegrp();
  set_basis(b1,b2,b3,b4);
}

Integral::Integral(StateIn& s) :
  SavableState(s)
{
  storage_used_ = 0;
  bs1_.restore_state(s);
  bs2_.restore_state(s);
  bs3_.restore_state(s);
  bs4_.restore_state(s);
  s.get(storage_);
  grp_ = MessageGrp::get_default_messagegrp();
}

Integral::Integral(const RefKeyVal&)
{
  storage_used_ = 0;
  storage_ = 0;
  grp_ = MessageGrp::get_default_messagegrp();
}

void
Integral::save_data_state(StateOut&o)
{
  bs1_.save_state(o);
  bs2_.save_state(o);
  bs3_.save_state(o);
  bs4_.save_state(o);
  o.put(storage_);
}

RefPetiteList
Integral::petite_list()
{
  return new PetiteList(bs1_, this);
}

RefPetiteList
Integral::petite_list(const RefGaussianBasisSet& gbs)
{
  return new PetiteList(gbs, this);
}

ShellRotation
Integral::shell_rotation(int am, SymmetryOperation& so, int pure)
{
  this->reference();
  ShellRotation r(am, so, this, pure);
  this->dereference();
  return r;
}

void
Integral::set_basis(const RefGaussianBasisSet &b1,
                    const RefGaussianBasisSet &b2,
                    const RefGaussianBasisSet &b3,
                    const RefGaussianBasisSet &b4)
{
  bs1_ = b1;
  bs2_ = b2;
  bs3_ = b3;
  bs4_ = b4;
  if (bs2_.null()) bs2_ = bs1_;
  if (bs3_.null()) bs3_ = bs2_;
  if (bs4_.null()) bs4_ = bs3_;
}

int
Integral::storage_unused()
{
  int tmp=storage_-storage_used_;
  return (tmp<0?0:tmp);
}
      

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
