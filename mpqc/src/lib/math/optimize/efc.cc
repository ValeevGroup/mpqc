
#ifdef __GNUC__
#pragma implementation
#endif

extern "C" {
#  include <math.h>
}

#include "efc.h"
#include <util/keyval/keyval.h>
#include <math/scmat/local.h>

#define CLASSNAME EFCOpt
#define PARENTS public Optimize
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>

/////////////////////////////////////////////////////////////////////////
// EFCOpt

void *
EFCOpt::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = Optimize::_castdown(cd);
  return do_castdowns(casts,cd);
}

EFCOpt::EFCOpt(const RefKeyVal&keyval):
  Optimize(keyval),
  maxabs_gradient(-1.0)
{
  nlp_ = keyval->describedclassvalue("function");
  update_ = keyval->describedclassvalue("update");
  
  convergence_ = keyval->doublevalue("convergence");
  if (keyval->error() != KeyVal::OK) convergence_ = 1.0e-6;
  
  accuracy_ = keyval->doublevalue("accuracy");
  if (keyval->error() != KeyVal::OK) accuracy_ = 0.0001;

  tstate = keyval->booleanvalue("transition_state");
  if (keyval->error() != KeyVal::OK) tstate = 0;

  if (tstate) {
    printf("\n  performing a transition state search\n\n");
  }
  
  RefSymmSCMatrix hessian(nlp_->dimension());
  // get a guess hessian from nlp
  nlp_->guess_hessian(hessian);
  
  // see if any hessian matrix elements have been given in the input
  if (keyval->exists("hessian")) {
    int n = hessian.n();
    for (int i=0; i<n; i++) {
      if (keyval->exists("hessian",i)) {
        for (int j=0; j<=i; j++) {
          double tmp = keyval->doublevalue("hessian",i,j);
          if (keyval->error() == KeyVal::OK) hessian(i,j) = tmp;
        }
      }
    }
  }
  hessian_ = hessian;
}

EFCOpt::EFCOpt(StateIn&s):
  SavableState(s,EFCOpt::class_desc_),
  Optimize(s)
{
  s.get(tstate);
  nlp_.restore_state(s);
  hessian_.restore_state(s);
  update_.restore_state(s);
  s.get(convergence_);
  s.get(accuracy_);
  s.get(take_newton_step_);
  s.get(maxabs_gradient);
}

EFCOpt::~EFCOpt()
{
}

void
EFCOpt::save_data_state(StateOut&s)
{
  Optimize::save_data_state(s);
  s.put(tstate);
  nlp_.save_state(s);
  hessian_.save_state(s);
  update_.save_state(s);
  s.put(convergence_);
  s.put(accuracy_);
  s.put(take_newton_step_);
  s.put(maxabs_gradient);
}

void
EFCOpt::init()
{
  Optimize::init();
  take_newton_step_ = 1;
  maxabs_gradient = -1.0;
}

int
EFCOpt::update()
{
  int i,j,ii,jj;
  
  // these are good candidates to be input options
  const double maxabs_gradient_to_desired_accuracy = 0.05;
  const double maxabs_gradient_to_next_desired_accuracy = 0.001;
  const double roundoff_error_factor = 1.1;

  // the gradient convergence criterion.
  double old_maxabs_gradient = maxabs_gradient;
  RefSCVector xcurrent;
  RefSCVector gcurrent;

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

    // compute the gradient convergence criterion now so i can see if
    // the accuracy needs to be tighter
    maxabs_gradient = gcurrent.maxabs();
    // compute the required accuracy
    accuracy_ = maxabs_gradient * maxabs_gradient_to_desired_accuracy;

    // The roundoff_error_factor is thrown in to allow for round off making
    // the current gcurrent.maxabs() a bit smaller than the previous,
    // which would make the current required accuracy less than the
    // gradient's actual accuracy and cause everything to be recomputed.
    accurate_enough = (nlp_->actual_gradient_accuracy() <=
                       accuracy_*roundoff_error_factor);

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
  
  // update the hessian
  if (update_.nonnull()) {
    update_->update(hessian_,nlp_,xcurrent,gcurrent);
  }

  // just take the Newton-Raphson step first iteration
  RefSCVector xdisp = -1.0*(hessian_.i() * gcurrent);
  // try steepest descent
  // RefSCVector xdisp = -1.0*gcurrent;
  RefSCVector xnext = xcurrent + xdisp;
  nlp_->set_x(xnext);
    
  // compute the convergence criteria
  double con_crit1 = fabs(xdisp.scalar_product(gcurrent));
  double con_crit2 = maxabs_gradient;
  double con_crit3 = xdisp.maxabs();

  return ((con_crit1 <= convergence_)
          && (con_crit2 <= convergence_)
          && (con_crit3 <= convergence_));
}

RefNLP0
EFCOpt::nlp()
{
  return nlp_;
}
