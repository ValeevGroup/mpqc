
#ifdef __GNUC__
#pragma implementation
#endif

#include <math.h>

#include <util/misc/formio.h>
#include <util/keyval/keyval.h>
#include <math/optimize/conv.h>

/////////////////////////////////////////////////////////////////////////
// Convergence

SavableState_REF_def(Convergence);
#define CLASSNAME Convergence
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#define PARENTS virtual_base public SavableState
#include <util/state/statei.h>
#include <util/class/classi.h>

void *
Convergence::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

Convergence::Convergence()
{
  set_defaults();
}

Convergence::Convergence(double tolerance)
{
  set_defaults();
  max_disp_ = tolerance;
  max_grad_ = tolerance;
  rms_disp_ = tolerance;
  rms_grad_ = tolerance;
  graddisp_ = tolerance;
}

Convergence::Convergence(StateIn&s):
  SavableState(s)
{
  s.get(use_max_disp_);
  s.get(use_max_grad_);
  s.get(use_rms_disp_);
  s.get(use_rms_grad_);
  s.get(use_graddisp_);
  s.get(max_disp_);
  s.get(max_grad_);
  s.get(rms_disp_);
  s.get(rms_grad_);
  s.get(graddisp_);
}

Convergence::Convergence(const RefKeyVal&keyval)
{
  use_max_disp_ = keyval->exists("max_disp");
  use_max_grad_ = keyval->exists("max_grad");
  use_rms_disp_ = keyval->exists("rms_disp");
  use_rms_grad_ = keyval->exists("rms_grad");
  use_graddisp_ = keyval->exists("graddisp");
  if (use_max_disp_) max_disp_ = keyval->doublevalue("max_disp");
  if (use_max_grad_) max_grad_ = keyval->doublevalue("max_grad");
  if (use_rms_disp_) rms_disp_ = keyval->doublevalue("rms_disp");
  if (use_rms_grad_) rms_grad_ = keyval->doublevalue("rms_grad");
  if (use_graddisp_) graddisp_ = keyval->doublevalue("graddisp");

  if (!use_max_disp_ && !use_max_grad_
      && !use_rms_disp_ && !use_rms_grad_
      && !use_graddisp_) {
      set_defaults();
    }
}

Convergence::~Convergence()
{
}

void
Convergence::save_data_state(StateOut&s)
{
  s.put(use_max_disp_);
  s.put(use_max_grad_);
  s.put(use_rms_disp_);
  s.put(use_rms_grad_);
  s.put(use_graddisp_);
  s.put(max_disp_);
  s.put(max_grad_);
  s.put(rms_disp_);
  s.put(rms_grad_);
  s.put(graddisp_);
}

void
Convergence::set_defaults()
{
  use_max_disp_ = 1;
  use_max_grad_ = 1;
  use_rms_disp_ = 0;
  use_rms_grad_ = 0;
  use_graddisp_ = 1;
  max_disp_ = 1.0e-6;
  max_grad_ = 1.0e-6;
  graddisp_ = 1.0e-6;
}

void
Convergence::get_x(const RefFunction &f)
{
  x_ = f->get_x();
}

void
Convergence::get_nextx(const RefFunction &f)
{
  nextx_ = f->get_x();
}

void
Convergence::get_grad(const RefFunction &f)
{
  grad_ = f->gradient();
}

int
Convergence::converged()
{
  int fail = 0;
  int pass = 0;

  RefSCVector disp;
  if (x_.nonnull() && nextx_.nonnull()) disp = nextx_ - x_;

  if (use_max_grad_ && grad_.nonnull()) {
      check_conv("Max Gradient     ", grad_.maxabs(), max_grad_, pass, fail);
    }
  if (use_rms_grad_ && grad_.nonnull()) {
      check_conv("RMS Gradient     ",
                 sqrt(grad_.scalar_product(grad_)), rms_grad_, pass, fail);
    }
  if (use_max_disp_ && disp.nonnull()) {
      check_conv("Max Displacement ", disp.maxabs(), max_disp_, pass, fail);
    }
  if (use_rms_disp_ && disp.nonnull()) {
      check_conv("RMS Displacement ",
                 sqrt(disp.scalar_product(disp)), rms_disp_, pass, fail);
    }
  if (use_graddisp_ && disp.nonnull() && grad_.nonnull()) {
      check_conv("Gradient*Displace", fabs(disp.scalar_product(grad_)),
                 graddisp_, pass, fail);
    }
  if (fail + pass == 0) {
      cerr << "ERROR: Convergence::converged: no applicable convergence tests"
           << endl;
      abort();
    }
  if (!fail) {
      cout << node0 << indent << "All convergence criteria have been met."
           << endl;
    }
  return !fail;
}

void
Convergence::check_conv(const char *heading,
                        double val, double bound,
                        int &pass, int &fail)
{
  int converged = val <= bound;
  cout << node0 << indent << heading << ": "
       << scprintf("%14.10f ", val)
       << scprintf("%14.10f  ", bound)
       << (converged?"yes":"no")
       << endl;
  if (converged) pass++;
  else fail++;
}

void
Convergence::reset()
{
  grad_ = 0;
  x_ = 0;
  nextx_ = 0;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")
