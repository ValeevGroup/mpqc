
extern "C" {
#  include <math.h>
}

#include "update.h"
#include <util/keyval/keyval.h>

#define CLASSNAME DFPUpdate
#define PARENTS public HessianUpdate
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>

#define CLASSNAME BFGSUpdate
#define PARENTS public DFPUpdate
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>

/////////////////////////////////////////////////////////////////////////
// DFPUpdate

void *
DFPUpdate::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = HessianUpdate::_castdown(cd);
  return do_castdowns(casts,cd);
}

DFPUpdate::DFPUpdate()
{
}

DFPUpdate::DFPUpdate(const RefKeyVal&keyval):
  HessianUpdate(keyval)
{
}

DFPUpdate::DFPUpdate(StateIn&s):
  HessianUpdate(s)
  maybe_SavableState(s)
{
  xprev.restore_state(s);
  gprev.restore_state(s);
}

DFPUpdate::~DFPUpdate()
{
}

void
DFPUpdate::save_data_state(StateOut&s)
{
  HessianUpdate::save_data_state(s);
  xprev.save_state(s);
  gprev.save_state(s);
}

void
DFPUpdate::update(RefSymmSCMatrix&ihessian,RefNLP2&nlp,
                  RefSCVector&xn,RefSCVector&gn)
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
  
  if (xprev.nonnull()) {
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

/////////////////////////////////////////////////////////////////////////
// BFGSUpdate

void *
BFGSUpdate::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = DFPUpdate::_castdown(cd);
  return do_castdowns(casts,cd);
}

BFGSUpdate::BFGSUpdate()
{
}

BFGSUpdate::BFGSUpdate(const RefKeyVal&keyval):
  DFPUpdate(keyval)
{
}

BFGSUpdate::BFGSUpdate(StateIn&s):
  DFPUpdate(s)
  maybe_SavableState(s)
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
BFGSUpdate::update(RefSymmSCMatrix&ihessian,RefNLP2&nlp,
                   RefSCVector&xn,RefSCVector&gn)
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
  
  if (xprev.nonnull()) {
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
