//
// dfp.cc
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
// DFPUpdate

static ClassDesc DFPUpdate_cd(
  typeid(DFPUpdate),"DFPUpdate",1,"public HessianUpdate",
  create<DFPUpdate>, create<DFPUpdate>, create<DFPUpdate>);

DFPUpdate::DFPUpdate()
{
}

DFPUpdate::DFPUpdate(const Ref<KeyVal>&keyval):
  HessianUpdate(keyval)
{
  if (keyval->exists("xprev") && keyval->exists("gprev")) {
    Ref<SCMatrixKit> k = SCMatrixKit::default_matrixkit();
    RefSCDimension dim = new SCDimension(keyval->count("xprev"));
    xprev = k->vector(dim);
    gprev = k->vector(dim);
    for (int i=0; i<dim.n(); i++) {
      xprev(i) = keyval->doublevalue("xprev",i);
      gprev(i) = keyval->doublevalue("gprev",i);
    }
  }
}

DFPUpdate::DFPUpdate(StateIn&s):
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

DFPUpdate::~DFPUpdate()
{
}

void
DFPUpdate::save_data_state(StateOut&s)
{
  HessianUpdate::save_data_state(s);
  SavableState::save_state(xprev.dim().pointer(),s);
  xprev.save(s);
  gprev.save(s);
}

void
DFPUpdate::update(const RefSymmSCMatrix&ihessian,const Ref<Function>&func,
                  const RefSCVector&xn,const RefSCVector&gn)
{
  RefSCVector xnew, gnew;

  // the update for the inverse hessian differs from the update for the
  // hessian in that xdisp and gdisp are exchanged
  if (inverse_hessian_) {
    xnew = xn;
    gnew = gn;
  } else {
    xnew = gn;
    gnew = xn;
  }
  
  if (xprev) {
    RefSCVector xdisp = xnew-xprev;
    RefSCVector gdisp = gnew-gprev;
    RefSCVector ihessian_gdisp = ihessian * gdisp;
    double gdisp_ihessian_gdisp = ihessian_gdisp.scalar_product(gdisp);
    double xdisp_gdisp = xdisp.scalar_product(gdisp);
    ihessian.accumulate(
        xdisp.symmetric_outer_product()*(1.0/xdisp_gdisp)
        - ihessian_gdisp.symmetric_outer_product()*(1.0/gdisp_ihessian_gdisp)
      );
    xprev.assign(xnew);
    gprev.assign(gnew);
  } else {
    xprev = xnew.copy();
    gprev = gnew.copy();
  }
}

void
DFPUpdate::apply_transform(const Ref<NonlinearTransform>& trans)
{
  if (trans == 0) return;
  HessianUpdate::apply_transform(trans);
  trans->transform_coordinates(xprev);
  trans->transform_gradient(gprev);
}

void
DFPUpdate::set_inverse(void)
{
  HessianUpdate::set_inverse();
  RefSCVector tmp;
  tmp = xprev;
  xprev = gprev;
  gprev = tmp;
}

/////////////////////////////////////////////////////////////////////////
// BFGSUpdate

static ClassDesc BFGSUpdate_cd(
  typeid(BFGSUpdate),"BFGSUpdate",1,"public HessianUpdate",
  create<BFGSUpdate>, create<BFGSUpdate>, create<BFGSUpdate>);

BFGSUpdate::BFGSUpdate()
{
}

BFGSUpdate::BFGSUpdate(const Ref<KeyVal>&keyval):
  DFPUpdate(keyval)
{
}

BFGSUpdate::BFGSUpdate(StateIn&s):
  SavableState(s),
  DFPUpdate(s)
{
}

BFGSUpdate::~BFGSUpdate()
{
}

void
BFGSUpdate::save_data_state(StateOut&s)
{
  DFPUpdate::save_data_state(s);
}

void
BFGSUpdate::update(const RefSymmSCMatrix&ihessian,const Ref<Function>&func,
                   const RefSCVector&xn,const RefSCVector&gn)
{
  RefSCVector xnew, gnew;

  // the update for the inverse hessian differs from the update for the
  // hessian in that xdisp and gdisp are exchanged
  if (inverse_hessian_) {
    xnew = xn;
    gnew = gn;
  } else {
    xnew = gn;
    gnew = xn;
  }
  
  if (xprev) {
    RefSCVector xdisp = xnew-xprev;
    RefSCVector gdisp = gnew-gprev;
    RefSCVector ihessian_gdisp = ihessian * gdisp;
    double gdisp_ihessian_gdisp = ihessian_gdisp.scalar_product(gdisp);
    double xdisp_gdisp = xdisp.scalar_product(gdisp);
    RefSCVector u =   xdisp*(1.0/xdisp_gdisp)
                      - ihessian_gdisp*(1.0/gdisp_ihessian_gdisp);
    ihessian.accumulate(
        // DFP part
        xdisp.symmetric_outer_product()*(1.0/xdisp_gdisp)
        - ihessian_gdisp.symmetric_outer_product()*(1.0/gdisp_ihessian_gdisp)
        // BFGS part
        + u.symmetric_outer_product() * gdisp_ihessian_gdisp
      );
    xprev.assign(xnew);
    gprev.assign(gnew);
  } else {
    xprev = xnew.copy();
    gprev = gnew.copy();
  }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
