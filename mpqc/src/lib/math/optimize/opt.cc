
#ifdef __GNUC__
#pragma implementation
#endif

extern "C" {
#  include <math.h>
}

#include "opt.h"
#include <util/keyval/keyval.h>

/////////////////////////////////////////////////////////////////////////
// Optimize

SavableState_REF_def(Optimize);

#define CLASSNAME Optimize
#define PARENTS virtual public SavableState
#include <util/state/statei.h>
#include <util/class/classia.h>
void *
Optimize::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

Optimize::Optimize()
{
}

Optimize::Optimize(StateIn&s):
  SavableState(s,Optimize::class_desc_)
{
  s.get(max_iterations_);
  n_iterations_ = 0;
}

Optimize::Optimize(KeyVal&keyval)
{
  max_iterations_ = keyval.intvalue("max_iterations");
  if (keyval.error() != KeyVal::OK) max_iterations_ = 10;
  n_iterations_ = 0;
}

Optimize::~Optimize()
{
}

void
Optimize::save_data_state(StateOut&s)
{
  s.put(max_iterations_);
}

void
Optimize::init()
{
  n_iterations_ = 0;
}

int
Optimize::optimize()
{
  int result;
  while((!(result = update())) && (n_iterations_ < max_iterations_))
      n_iterations_++;
  return result;
}

/////////////////////////////////////////////////////////////////////////
// LineOpt

SavableState_REF_def(LineOpt);

#define CLASSNAME LineOpt
#define PARENTS public Optimize
#include <util/state/statei.h>
#include <util/class/classia.h>
void *
LineOpt::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = Optimize::_castdown(cd);
  return do_castdowns(casts,cd);
}

LineOpt::LineOpt()
{
}

LineOpt::LineOpt(StateIn&s):
  SavableState(s,LineOpt::class_desc_),
  Optimize(s)
{
  search_direction_.restore_state(s);
}

LineOpt::LineOpt(KeyVal&keyval):
  Optimize(keyval)
{
}

LineOpt::~LineOpt()
{
}

void
LineOpt::save_data_state(StateOut&s)
{
  Optimize::save_data_state(s);
  search_direction_.save_state(s);
}

void
LineOpt::set_search_direction(RefSCVector&s)
{
  search_direction_ = s.copy();
}

/////////////////////////////////////////////////////////////////////////
// IHessianUpdate

SavableState_REF_def(IHessianUpdate);

#   define CLASSNAME IHessianUpdate
#   include <util/state/statei.h>
#   include <util/class/classia.h>
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

#define CLASSNAME DFPUpdate
#define PARENTS public IHessianUpdate
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
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

#define CLASSNAME BFGSUpdate
#define PARENTS public DFPUpdate
#define HAVE_CTOR
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
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

/////////////////////////////////////////////////////////////////////////
// QNewtonOpt

#define CLASSNAME QNewtonOpt
#define PARENTS public Optimize
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
QNewtonOpt::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = Optimize::_castdown(cd);
  return do_castdowns(casts,cd);
}

QNewtonOpt::QNewtonOpt(KeyVal&keyval):
  Optimize(keyval),
  maxabs_gradient(-1.0)
{
  nlp_ = keyval.describedclassvalue("function");
  update_ = keyval.describedclassvalue("update");
  convergence_ = keyval.doublevalue("convergence");
  if (keyval.error() != KeyVal::OK) convergence_ = 1.0e-6;
  lineopt_ = keyval.describedclassvalue("lineopt");
  accuracy_ = keyval.doublevalue("accuracy");
  if (keyval.error() != KeyVal::OK) accuracy_ = 0.0001;

  RefSymmSCMatrix hessian(nlp_->dimension());
  // get a guess hessian from nlp
  nlp_->guess_hessian(hessian);
  // see if any hessian matrix elements have been given in the input
  if (keyval.exists("hessian")) {
      int n = hessian.n();
      for (int i=0; i<n; i++) {
          if (keyval.exists("hessian",i)) {
              for (int j=0; j<=i; j++) {
                  double tmp = keyval.doublevalue("hessian",i,j);
                  if (keyval.error() == KeyVal::OK) hessian(i,j) = tmp;
                }
            }
        }
    }
  ihessian_ = hessian.i();
}

QNewtonOpt::QNewtonOpt(StateIn&s):
  SavableState(s,QNewtonOpt::class_desc_),
  Optimize(s)
{
  nlp_.restore_state(s);
  ihessian_.restore_state(s);
  update_.restore_state(s);
  s.get(convergence_);
  s.get(accuracy_);
  s.get(take_newton_step_);
  s.get(maxabs_gradient);
  lineopt_.restore_state(s);
}

QNewtonOpt::~QNewtonOpt()
{
}

void
QNewtonOpt::save_data_state(StateOut&s)
{
  Optimize::save_data_state(s);
  nlp_.save_state(s);
  ihessian_.save_state(s);
  update_.save_state(s);
  s.put(convergence_);
  s.put(accuracy_);
  s.put(take_newton_step_);
  s.put(maxabs_gradient);
  lineopt_.save_state(s);
}

void
QNewtonOpt::init()
{
  Optimize::init();
  take_newton_step_ = 1;
  maxabs_gradient = -1.0;
}

int
QNewtonOpt::update()
{
  // these are good candidates to be input options
  const double maxabs_gradient_to_desired_accuracy = 0.05;
  const double maxabs_gradient_to_next_desired_accuracy = 0.001;
  const double roundoff_error_factor = 1.1;

  // the gradient convergence criterion.
  double old_maxabs_gradient = maxabs_gradient;
  RefSCVector xcurrent;
  RefSCVector gcurrent;
  double value;

  SCostream::cout.flush();
    
  // get the next gradient at the required level of accuracy.
  // usually only one pass is needed, unless we happen to find
  // that the accuracy was set too low.
  int accurate_enough;
  do {
      // compute the current point
      nlp_->set_desired_gradient_accuracy(accuracy_);

      xcurrent = nlp_->get_x().copy();
      gcurrent = nlp_->gradient().copy();
      value = nlp_->value();

      // compute the gradient convergence criterion now so i can see if
      // the accuracy needs to be tighter
      maxabs_gradient = gcurrent.maxabs();
      // compute the required accuracy
      accuracy_ = maxabs_gradient * maxabs_gradient_to_desired_accuracy;

      // The roundoff_error_factor is thrown in to allow for round off making
      // the current gcurrent.maxabs() a bit smaller than the previous,
      // which would make the current required accuracy less than the
      // gradient's actual accuracy and cause everything to be recomputed.
      accurate_enough = (
        nlp_->actual_gradient_accuracy() <= accuracy_*roundoff_error_factor);

      if (!accurate_enough) {
          printf("NOTICE: nlp_->actual_gradient_accuracy() > accuracy_:\n");
          printf("  nlp_->actual_gradient_accuracy() = %15.8f\n",
                 nlp_->actual_gradient_accuracy());
          printf("  accuracy_ = %15.8f\n", accuracy_);
          fflush(stdout);
        }
    } while(!accurate_enough);

  if (old_maxabs_gradient >= 0.0 && old_maxabs_gradient < maxabs_gradient) {
      printf("NOTICE: maxabs_gradient increased from %8.4e to %8.4e\n",
             old_maxabs_gradient, maxabs_gradient);
      fflush(stdout);
    }

  // make the next gradient computed more accurate, since it will
  // be smaller
  accuracy_ = maxabs_gradient * maxabs_gradient_to_next_desired_accuracy;
  
  if (!take_newton_step_ && lineopt_.nonnull()) {
      // see if the line min step is really needed
//       if (line min really needed) {
//           take_newton_step_ = (lineopt_->update() == 1);
//           // maybe check for convergence and return here
//         }
    }

  // update the hessian
  if (update_.nonnull()) {
      update_->update(ihessian_,nlp_,xcurrent,gcurrent);
    }

  // take the step
  RefSCVector xdisp = -1.0*(ihessian_ * gcurrent);
  RefSCVector xnext = xcurrent + xdisp;
  nlp_->set_x(xnext);

  // do a line min step next time
  if (lineopt_.nonnull()) take_newton_step_ = 0;

  // compute the convergence criteria
  double con_crit1 = fabs(xdisp.scalar_product(gcurrent));
  double con_crit2 = maxabs_gradient;
  double con_crit3 = xdisp.maxabs();

  return   (con_crit1 <= convergence_)
            && (con_crit2 <= convergence_)
            && (con_crit3 <= convergence_);
}

RefNLP0
QNewtonOpt::nlp()
{
  return nlp_;
}
