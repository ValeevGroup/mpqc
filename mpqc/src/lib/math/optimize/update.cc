
extern "C" {
#  include <math.h>
}

#include "opt.h"
#include <util/keyval/keyval.h>

SavableState_REF_def(IHessianUpdate);
#define CLASSNAME IHessianUpdate
#include <util/state/statei.h>
#include <util/class/classia.h>

#define CLASSNAME DFPUpdate
#define PARENTS public IHessianUpdate
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
// IHessianUpdate

void *
IHessianUpdate::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

IHessianUpdate::IHessianUpdate()
{
}

IHessianUpdate::IHessianUpdate(StateIn&s):
  SavableState(s,IHessianUpdate::class_desc_)
{
}

IHessianUpdate::IHessianUpdate(KeyVal&keyval)
{
}

IHessianUpdate::~IHessianUpdate()
{
}

void
IHessianUpdate::save_data_state(StateOut&s)
{
}

/////////////////////////////////////////////////////////////////////////
// DFPUpdate

void *
DFPUpdate::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = IHessianUpdate::_castdown(cd);
  return do_castdowns(casts,cd);
}

DFPUpdate::DFPUpdate()
{
}

DFPUpdate::DFPUpdate(KeyVal&keyval):
  IHessianUpdate(keyval)
{
}

DFPUpdate::DFPUpdate(StateIn&s):
  SavableState(s,DFPUpdate::class_desc_),
  IHessianUpdate(s)
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
  IHessianUpdate::save_data_state(s);
  xprev.save_state(s);
  gprev.save_state(s);
}

void
DFPUpdate::update(RefSymmSCMatrix&ihessian,RefNLP2&nlp,
                  RefSCVector&xnew,RefSCVector&gnew)
{
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
    }
  xprev = xnew;
  gprev = gnew;
}

/////////////////////////////////////////////////////////////////////////
// BFGSUpdate

void *
BFGSUpdate::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = IHessianUpdate::_castdown(cd);
  return do_castdowns(casts,cd);
}

BFGSUpdate::BFGSUpdate()
{
}

BFGSUpdate::BFGSUpdate(KeyVal&keyval):
  DFPUpdate(keyval)
{
}

BFGSUpdate::BFGSUpdate(StateIn&s):
  SavableState(s,BFGSUpdate::class_desc_),
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
BFGSUpdate::update(RefSymmSCMatrix&ihessian,RefNLP2&nlp,
                   RefSCVector&xnew,RefSCVector&gnew)
{
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
    }
  xprev = xnew;
  gprev = gnew;
}
