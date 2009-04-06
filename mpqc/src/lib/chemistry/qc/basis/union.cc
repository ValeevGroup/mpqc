//
// union.cc
//
// Copyright (C) 2009 Edward Valeev
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

#ifdef __GNUG__
#pragma implementation
#endif

// includes go here
#include <chemistry/qc/basis/union.h>

using namespace sc;

ClassDesc
UnionBasisSet::class_desc_(typeid(UnionBasisSet),"UnionBasisSet",1,"public GaussianBasisSet",
  0, create<UnionBasisSet>, create<UnionBasisSet>);

UnionBasisSet::UnionBasisSet(const Ref<KeyVal>& keyval) :
  GaussianBasisSet(
    *( GaussianBasisSetSum(
         require_dynamic_cast<GaussianBasisSet*>(keyval->describedclassvalue("basis1"),"UnionBasisSet:basis1"),
         require_dynamic_cast<GaussianBasisSet*>(keyval->describedclassvalue("basis2"),"UnionBasisSet:basis2")
       ).bs12()) )
{
}

UnionBasisSet::UnionBasisSet(StateIn& si) :
  GaussianBasisSet(si)
{
}

void
UnionBasisSet::save_data_state(StateOut& so)
{
  GaussianBasisSet::save_data_state(so);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
