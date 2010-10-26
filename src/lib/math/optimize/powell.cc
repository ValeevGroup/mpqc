//
// powell.cc
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

#include <math.h>

#include <util/state/stateio.h>
#include <math/optimize/update.h>
#include <math/optimize/transform.h>
#include <util/keyval/keyval.h>

using namespace sc;

/////////////////////////////////////////////////////////////////////////
// PowellUpdate

static ClassDesc PowellUpdate_cd(
  typeid(PowellUpdate),"PowellUpdate",1,"public HessianUpdate",
  create<PowellUpdate>, create<PowellUpdate>, create<PowellUpdate>);

PowellUpdate::PowellUpdate()
{
}

PowellUpdate::PowellUpdate(const Ref<KeyVal>&keyval):
  HessianUpdate(keyval)
{
}

PowellUpdate::PowellUpdate(StateIn&s):
  SavableState(s),
  HessianUpdate(s)
{
  Ref<SCMatrixKit> k = SCMatrixKit::default_matrixkit();
  RefSCDimension dim;
  dim << SavableState::restore_state(s);
  xprev = k->vector(dim);
  gprev = k->vector(dim);
  xprev.restore(s);
  gprev.restore(s);
}

PowellUpdate::~PowellUpdate()
{
}

void
PowellUpdate::save_data_state(StateOut&s)
{
  HessianUpdate::save_data_state(s);
  SavableState::save_state(xprev.dim().pointer(),s);
  xprev.save(s);
  gprev.save(s);
}

void
PowellUpdate::update(const RefSymmSCMatrix&hessian,const Ref<Function>&func,
                     const RefSCVector&xn,const RefSCVector&gn)
{
  RefSCVector xnew, gnew;

  // the update for the inverse hessian differs from the update for the
  // hessian in that xdisp and gdisp are exchanged
  // test this...it may only be true for the Broyden family of updates
  if (!inverse_hessian_) {
    xnew = xn;
    gnew = gn;
  } else {
    xnew = gn;
    gnew = xn;
  }
  
  if (xprev.nonnull()) {
    RefSCVector xdisp = xnew-xprev;
    RefSCVector gdisp = gnew-gprev-hessian*xdisp;
    double xdisp_xdisp = xdisp.scalar_product(xdisp);
    double gdisp_xdisp = gdisp.scalar_product(xdisp);
    hessian.accumulate(
      xdisp.symmetric_outer_product() *
      (-gdisp_xdisp/(xdisp_xdisp*xdisp_xdisp))
    );
    RefSCMatrix temp =
      (gdisp.outer_product(xdisp) + xdisp.outer_product(gdisp)) *
      (0.5/xdisp_xdisp);
    hessian.accumulate_symmetric_sum(temp);
    xprev.assign(xnew);
    gprev.assign(gnew);
  } else {
    xprev = xnew.copy();
    gprev = gnew.copy();
  }
}

void
PowellUpdate::apply_transform(const Ref<NonlinearTransform>& trans)
{
  if (trans.null()) return;
  HessianUpdate::apply_transform(trans);
  trans->transform_coordinates(xprev);
  trans->transform_gradient(gprev);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
