
#ifdef __GNUC__
#pragma implementation
#endif

extern "C" {
#  include <math.h>
}

#include <math/optimize/efc.h>
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
  update_ = keyval->describedclassvalue("update");
  
  convergence_ = keyval->doublevalue("convergence");
  if (keyval->error() != KeyVal::OK) convergence_ = 1.0e-6;
  
  accuracy_ = keyval->doublevalue("accuracy");
  if (keyval->error() != KeyVal::OK) accuracy_ = 0.0001;

  tstate = keyval->booleanvalue("transition_state");
  if (keyval->error() != KeyVal::OK) tstate = 0;

  modef = keyval->booleanvalue("mode_following");
  if (keyval->error() != KeyVal::OK) modef = 0;

  if (tstate) {
    printf("\n  performing a transition state search\n\n");
  }
  
  RefSymmSCMatrix hessian(dimension(),matrixkit());
  // get a guess hessian from function
  function()->guess_hessian(hessian);
  
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
  last_mode_ = 0;
}

EFCOpt::EFCOpt(StateIn&s):
  Optimize(s)
  maybe_SavableState(s)
{
  s.get(tstate);
  s.get(modef);
  hessian_ = matrixkit()->symmmatrix(dimension());
  hessian_.restore(s);
  update_.restore_state(s);
  last_mode_ = matrixkit()->vector(dimension());
  last_mode_.restore(s);
  s.get(convergence_);
  s.get(accuracy_);
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
  s.put(modef);
  hessian_.save(s);
  update_.save_state(s);
  last_mode_.save(s);
  s.put(convergence_);
  s.put(accuracy_);
  s.put(maxabs_gradient);
}

void
EFCOpt::init()
{
  Optimize::init();
  maxabs_gradient = -1.0;
}

int
EFCOpt::update()
{
  int i,j;
  
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
    function()->set_desired_gradient_accuracy(accuracy_);
    
    xcurrent = function()->get_x().copy();
    gcurrent = function()->gradient().copy();

    // compute the gradient convergence criterion now so i can see if
    // the accuracy needs to be tighter
    maxabs_gradient = gcurrent.maxabs();
    // compute the required accuracy
    accuracy_ = maxabs_gradient * maxabs_gradient_to_desired_accuracy;

    // The roundoff_error_factor is thrown in to allow for round off making
    // the current gcurrent.maxabs() a bit smaller than the previous,
    // which would make the current required accuracy less than the
    // gradient's actual accuracy and cause everything to be recomputed.
    accurate_enough = (function()->actual_gradient_accuracy() <=
                       accuracy_*roundoff_error_factor);

    if (!accurate_enough) {
      printf("NOTICE: function()->actual_gradient_accuracy() > accuracy_:\n");
      printf("  function()->actual_gradient_accuracy() = %15.8f\n",
             function()->actual_gradient_accuracy());
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
    update_->update(hessian_,function(),xcurrent,gcurrent);
  }

  // begin efc junk
  // first diagonalize hessian
  RefSCMatrix evecs(dimension(),dimension(),matrixkit());
  RefDiagSCMatrix evals(dimension(),matrixkit());

  hessian_.diagonalize(evals,evecs);
  evals.print("hessian eigenvalues");
  evecs.print("hessian eigenvectors");

  // form gradient to local hessian modes F = Ug
  RefSCVector F = evecs.t() * gcurrent;
  F.print("F");

  // figure out if hessian has the right number of negative eigenvalues
  int ncoord = evals.n();
  int npos=0,nneg=0;
  for (i=0; i < ncoord; i++) {
    if (evals.get_element(i) >= 0.0) npos++;
    else nneg++;
  }

  RefSCVector xdisp(dimension(),matrixkit());
  xdisp.assign(0.0);
  
  // for now, we always take the P-RFO for tstate (could take NR if
  // nneg==1, but we won't make that an option yet)
  if (tstate) {
    int mode = 0;

    if (modef) {
      // which mode are we following.  find mode with maximum overlap with
      // last mode followed
      if (last_mode_.nonnull()) {
        double overlap=0;
        for (i=0; i < ncoord; i++) {
          double S=0;
          for (j=0; j < ncoord; j++) {
            S += last_mode_.get_element(j)*evecs.get_element(j,i);
          }
          S = fabs(S);
          if (S > overlap) {
            mode = i;
            overlap = S;
          }
        }
      } else {
        last_mode_ = matrixkit()->vector(dimension());
      
        // find mode with max component = coord 0 which should be the
        // mode being followed
        double comp=0;
        for (i=0; i < ncoord; i++) {
          double S = fabs(evecs.get_element(0,i));
          if (S>comp) {
            mode=i;
            comp=S;
          }
        }
      }
    
      for (i=0; i < ncoord; i++)
        last_mode_(i) = evecs(i,mode);

      printf("\n following mode %d\n",mode);
    }
    
    double bk = evals(mode);
    double Fk = F(mode);
    double lambda_p = 0.5*bk + 0.5*sqrt(bk*bk + 4*Fk*Fk);
    
    double lambda_n;
    double nlambda=1.0;
    do {
      lambda_n=nlambda;
      nlambda=0;
      for (i=0; i < ncoord; i++) {
        if (i==mode) continue;
        
        nlambda += F.get_element(i)*F.get_element(i) /
                    (lambda_n - evals.get_element(i));
      }
    } while(fabs(nlambda-lambda_n) > 1.0e-8);

    printf("lambda_p = %g\n",lambda_p);
    printf("lambda_n = %g\n",lambda_n);

    // form Xk
    double Fkobkl = F(mode)/(evals(mode)-lambda_p);
    for (j=0; j < F.n(); j++)
      xdisp(j) = xdisp(j) - evecs(j,mode) * Fkobkl;
    
    // form displacement x = sum -Fi*Vi/(bi-lam)
    for (i=0; i < F.n(); i++) {
      if (i==mode) continue;
      
      double Fiobil = F(i) / (evals(i)-lambda_n);
      for (j=0; j < F.n(); j++) {
        xdisp(j) = xdisp(j) - evecs(j,i) * Fiobil;
      }
    }
    
 // minimum search
  } else {
    // evaluate lambda
    double lambda;
    double nlambda=1.0;
    do {
      lambda=nlambda;
      nlambda=0;
      for (i=0; i < F.n(); i++) {
        double Fi = F(i);
        nlambda += Fi*Fi / (lambda - evals.get_element(i));
      }
    } while(fabs(nlambda-lambda) > 1.0e-8);
    printf("lambda = %g\n",lambda);

  // form displacement x = sum -Fi*Vi/(bi-lam)
    for (i=0; i < F.n(); i++) {
      double Fiobil = F(i) / (evals(i)-lambda);
      for (j=0; j < F.n(); j++) {
        xdisp(j) = xdisp(j) - evecs(j,i) * Fiobil;
      }
    }
  }

  // scale the displacement vector if it's too large
  double tot = sqrt(xdisp.scalar_product(xdisp));
  double maxstepsize=0.6;
  if (tot > maxstepsize) {
    double scal = maxstepsize/tot;
    printf("\n stepsize of %f is too big, scaling by %f\n",tot,scal);
    xdisp.scale(scal);
    tot *= scal;
  }
  printf("\n taking step of size %f\n",tot);
  fflush(stdout);
                    
  xdisp.print("xdisp");

  // try steepest descent
  // RefSCVector xdisp = -1.0*gcurrent;
  RefSCVector xnext = xcurrent + xdisp;
  function()->set_x(xnext);
    
  // compute the convergence criteria
  double con_crit1 = fabs(xdisp.scalar_product(gcurrent));
  double con_crit2 = maxabs_gradient;
  double con_crit3 = xdisp.maxabs();

  return ((con_crit1 <= convergence_)
          && (con_crit2 <= convergence_)
          && (con_crit3 <= convergence_));
}
