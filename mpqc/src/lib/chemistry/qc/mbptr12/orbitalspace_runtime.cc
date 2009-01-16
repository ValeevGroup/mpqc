//
// orbitalspace_runtime.cc
//
// Copyright (C) 2008 Edward Valeev
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
#include<chemistry/qc/mbptr12/orbitalspace_runtime.h>

using namespace sc;

static ClassDesc OrbitalSpaceRuntime_cd(
  typeid(OrbitalSpaceRuntime),"OrbitalSpaceRuntime",1,"virtual public SavableState",
  create<OrbitalSpaceRuntime>, 0, create<OrbitalSpaceRuntime>);

OrbitalSpaceRuntime::OrbitalSpaceRuntime() {}

OrbitalSpaceRuntime::OrbitalSpaceRuntime(StateIn& si) : SavableState(si)
{
  basis_to_aospace_map_ = BasisToAOSpaceMap::restore_instance(si);
}

void
OrbitalSpaceRuntime::save_data_state(StateOut& so)
{
  BasisToAOSpaceMap::save_instance(basis_to_aospace_map_,so);
}

void
OrbitalSpaceRuntime::declare_aospace(const Ref<OrbitalSpace>& aospace)
{
  basis_to_aospace_map_->add(aospace->basis(),aospace);
}

Ref<OrbitalSpace>
OrbitalSpaceRuntime::aospace(const Ref<GaussianBasisSet>& basis)
{
  return basis_to_aospace_map_->value(basis);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End: