
extern "C" {
#  include <math.h>
}

#include <math/optimize/update.h>
#include <util/keyval/keyval.h>

#define CLASSNAME PowellUpdate
#define PARENTS public HessianUpdate
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>

/////////////////////////////////////////////////////////////////////////
// PowellUpdate

void *
PowellUpdate::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = HessianUpdate::_castdown(cd);
  return do_castdowns(casts,cd);
}

PowellUpdate::PowellUpdate()
{
}

PowellUpdate::PowellUpdate(const RefKeyVal&keyval):
  HessianUpdate(keyval)
{
}

PowellUpdate::PowellUpdate(StateIn&s):
  HessianUpdate(s)
  maybe_SavableState(s)
{
  xprev.restore_state(s);
  gprev.restore_state(s);
}

PowellUpdate::~PowellUpdate()
{
}

void
PowellUpdate::save_data_state(StateOut&s)
{
  HessianUpdate::save_data_state(s);
  xprev.save_state(s);
  gprev.save_state(s);
}

void
PowellUpdate::update(RefSymmSCMatrix&hessian,RefNLP2&nlp,
                  RefSCVector&xn,RefSCVector&gn)
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
