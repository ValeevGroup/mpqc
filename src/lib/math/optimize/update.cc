//
// update.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
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

#include <math.h>

#include <util/state/stateio.h>
#include <math/optimize/update.h>
#include <util/keyval/keyval.h>

using namespace sc;

/////////////////////////////////////////////////////////////////////////
// HessianUpdate

static ClassDesc HessianUpdate_cd(
  typeid(HessianUpdate),"HessianUpdate",1,"virtual public SavableState",
  0, 0, 0);

HessianUpdate::HessianUpdate() : inverse_hessian_(0)
{
}

HessianUpdate::HessianUpdate(StateIn&s):
  SavableState(s)
{
  s.get(inverse_hessian_);
}

HessianUpdate::HessianUpdate(const Ref<KeyVal>&keyval) :
  inverse_hessian_(0)
{
}

HessianUpdate::~HessianUpdate()
{
}

void
HessianUpdate::save_data_state(StateOut&s)
{
  s.put(inverse_hessian_);
}

void
HessianUpdate::set_inverse(void)
{
  inverse_hessian_ = 1;
}

void
HessianUpdate::apply_transform(const Ref<NonlinearTransform>&)
{
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
