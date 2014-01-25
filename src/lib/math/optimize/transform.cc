//
// transform.cc
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

#include <math/optimize/transform.h>
#include <util/misc/formio.h>

using namespace sc;

////////////////////////////////////////////////////////////////////////////

NonlinearTransform::~NonlinearTransform()
{
}

void
NonlinearTransform::transform_gradient(const RefSCVector& g)
{
  if (g.null()) return;
  g.assign(linear_transform_ * g);
}

void
NonlinearTransform::transform_hessian(const RefSymmSCMatrix& h)
{
  if (h.null()) return;
  ExEnv::out0() << indent
       << "WARNING: NonlinearTransform::transform_hessian: "
       << "using linear transform\n";
  RefSymmSCMatrix newh = h->clone();
  newh.assign(0.0);
  newh->accumulate_transform(linear_transform_.pointer(), h.pointer());
  h.assign(newh);
}

void
NonlinearTransform::transform_ihessian(const RefSymmSCMatrix &ih)
{
  if (ih.null()) return;
  RefSymmSCMatrix h(ih.gi());
  transform_hessian(h);
  ih.assign(h.gi());
}

////////////////////////////////////////////////////////////////////////////

IdentityTransform::~IdentityTransform()
{
}

void
IdentityTransform::transform_coordinates(const RefSCVector& x)
{
}

void
IdentityTransform::transform_gradient(const RefSCVector& g)
{
}

void
IdentityTransform::transform_hessian(const RefSymmSCMatrix& h)
{
}

void
IdentityTransform::transform_ihessian(const RefSymmSCMatrix &ih)
{
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
