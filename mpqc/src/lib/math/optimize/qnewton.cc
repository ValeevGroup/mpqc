
#ifdef __GNUC__
#pragma implementation
#endif

extern "C" {
#  include <math.h>
}

#include <math/optimize/qnewton.h>
#include <util/keyval/keyval.h>

#define CLASSNAME QNewtonOpt
#define PARENTS public Optimize
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>

/////////////////////////////////////////////////////////////////////////
// QNewtonOpt

void *
QNewtonOpt::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = Optimize::_castdown(cd);
  return do_castdowns(casts,cd);
}

QNewtonOpt::QNewtonOpt(const RefKeyVal&keyval):
  Optimize(keyval),
  maxabs_gradient(-1.0)
{
  nlp_ = keyval->describedclassvalue("function");
  update_ = keyval->describedclassvalue("update");
  if (update_.nonnull()) update_->set_inverse();
  convergence_ = keyval->doublevalue("convergence");
  if (keyval->error() != KeyVal::OK) convergence_ = 1.0e-6;
  lineopt_ = keyval->describedclassvalue("lineopt");
  accuracy_ = keyval->doublevalue("accuracy");
  if (keyval->error() != KeyVal::OK) accuracy_ = 0.0001;

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
  ihessian_ = nlp_->inverse_hessian(hessian);
}

QNewtonOpt::QNewtonOpt(StateIn&s):
  Optimize(s)
  maybe_SavableState(s)
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

  cout.flush();
    
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
  // scale the displacement vector if it's too large
  double tot = sqrt(xdisp.scalar_product(xdisp));
  double maxstepsize=0.3;
  if (tot > maxstepsize) {
    double scal = maxstepsize/tot;
    printf("\n stepsize of %f is too big, scaling by %f\n",tot,scal);
    xdisp.scale(scal);
    tot *= scal;
  }
  printf("\n taking step of size %f\n",tot);
  fflush(stdout);
  
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
