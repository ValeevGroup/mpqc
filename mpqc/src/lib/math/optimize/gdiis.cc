
#ifdef __GNUC__
#pragma implementation
#endif

extern "C" {
#  include <math.h>
}

#include "gdiis.h"
#include <util/keyval/keyval.h>
#include <math/scmat/local.h>

#define CLASSNAME GDIISOpt
#define PARENTS public Optimize
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>

/////////////////////////////////////////////////////////////////////////
// GDIISOpt

void *
GDIISOpt::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = Optimize::_castdown(cd);
  return do_castdowns(casts,cd);
}

GDIISOpt::GDIISOpt(const RefKeyVal&keyval):
  Optimize(keyval),
  maxabs_gradient(-1.0),
  diis_iter(0)
{
  nsave = keyval->intvalue("ngdiis");
  if (keyval->error() != KeyVal::OK) nsave = 5;
  
  nlp_ = keyval->describedclassvalue("function");
  update_ = keyval->describedclassvalue("update");
  update_->set_inverse();
  
  convergence_ = keyval->doublevalue("convergence");
  if (keyval->error() != KeyVal::OK) convergence_ = 1.0e-6;
  
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
  ihessian_ = hessian.gi();

  dim_ = nlp_->dimension();
  coords_ = new RefSCVector[nsave];
  grad_ = new RefSCVector[nsave];
  error_ = new RefSCVector[nsave];

//   LocalSCDimension *ldim =
//     LocalSCDimension::require_castdown(dim_.pointer(),
//                                     "GDIISOpt::GDIISOpt(const RefKeyVal&)");
//   
//   for (int i=0; i < nsave; i++) {
//     coords_[i] = new LocalSCVector(ldim);  coords_[i]->assign(0.0);
//     grad_[i] = new LocalSCVector(ldim);    grad_[i]->assign(0.0);
//     error_[i] = new LocalSCVector(ldim);   error_[i]->assign(0.0);
//   }

  for (int i=0; i < nsave; i++) {
    coords_[i] = dim_->create_vector(); coords_[i]->assign(0.0);
    grad_[i] = dim_->create_vector(); grad_[i]->assign(0.0);
    error_[i] = dim_->create_vector(); error_[i]->assign(0.0);
  }
}

GDIISOpt::GDIISOpt(StateIn&s):
  Optimize(s)
  maybe_SavableState(s)
{
  s.get(nsave);
  s.get(diis_iter);
  dim_.restore_state(s);
  nlp_.restore_state(s);
  ihessian_.restore_state(s);
  update_.restore_state(s);
  s.get(convergence_);
  s.get(accuracy_);
  s.get(maxabs_gradient);
  coords_ = new RefSCVector[nsave];
  grad_ = new RefSCVector[nsave];
  error_ = new RefSCVector[nsave];
  for (int i=0; i < nsave; i++) {
    coords_[i].restore_state(s);
    grad_[i].restore_state(s);
    error_[i].restore_state(s);
  }
}

GDIISOpt::~GDIISOpt()
{
}

void
GDIISOpt::save_data_state(StateOut&s)
{
  Optimize::save_data_state(s);
  s.put(nsave);
  s.put(diis_iter);
  dim_.save_state(s);
  nlp_.save_state(s);
  ihessian_.save_state(s);
  update_.save_state(s);
  s.put(convergence_);
  s.put(accuracy_);
  s.put(maxabs_gradient);
  for (int i=0; i < nsave; i++) {
    coords_[i].save_state(s);
    grad_[i].save_state(s);
    error_[i].save_state(s);
  }
}

void
GDIISOpt::init()
{
  Optimize::init();
  maxabs_gradient = -1.0;
}

int
GDIISOpt::update()
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
    update_->update(ihessian_,nlp_,xcurrent,gcurrent);
  }

  diis_iter++;

  int howmany = (diis_iter < nsave) ? diis_iter : nsave;

  if (diis_iter <= nsave) {
    coords_[diis_iter-1] = xcurrent;
    grad_[diis_iter-1] = gcurrent;
  } else {
    for (i=0; i < nsave-1; i++) {
      coords_[i] = coords_[i+1];
      grad_[i] = grad_[i+1];
    }
    coords_[nsave-1] = xcurrent;
    grad_[nsave-1] = gcurrent;
  }
  
  // take the step
  if (diis_iter==1 || maxabs_gradient > 0.005) {
    // just take the Newton-Raphson step first iteration
    // RefSCVector xdisp = -1.0*(ihessian_ * gcurrent);
    // try steepest descent
    RefSCVector xdisp = -1.0*gcurrent;
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

  // form the error vectors
  for (i=0; i < howmany; i++)
    error_[i] = -1.0*(ihessian_ * grad_[i]);

  // and form the A matrix
  RefSCMatrix A;
  RefSCVector coeff;
  double determ;
  int ntry=0;

  do {
    int num = howmany-ntry;

    RefLocalSCDimension size = new LocalSCDimension(num+1);
    A = new LocalSCMatrix(size.pointer(),size.pointer());
    coeff = new LocalSCVector(size.pointer());

    for (ii=0, i=ntry; i < howmany; i++,ii++) {
      coeff(ii) = 0;
      for (j=ntry,jj=0; j <= i; j++,jj++) {
        A(ii,jj) = error_[i].scalar_product(error_[j]);
        A(jj,ii) = A(ii,jj);
      }
    }

    A->scale(1.0/A(0,0));

    coeff(num) = 1.0;
    for (i=0; i < num; i++)
      A(num,i) = A(i,num) = 1.0;

    A(num,num) = 0;

    ntry++;

  } while ((determ=fabs(A.solve_lin(coeff))) < 1.0e-12);

  RefSCVector xstar = xcurrent.dim()->create_vector();
  RefSCVector delstar = gcurrent.dim()->create_vector();

  xstar.assign(0.0);
  delstar.assign(0.0);

  for (i=0,ii=ntry-1; ii < howmany; i++,ii++) {
    RefSCVector tmp = grad_[ii].copy(); tmp.scale(coeff[i]);
    delstar.accumulate(tmp);
    tmp = coords_[ii].copy(); tmp.scale(coeff[i]);
    xstar.accumulate(tmp);
  }

  RefSCVector xdisp = xstar - xcurrent - ihessian_*delstar;
  RefSCVector xnext = xcurrent + xdisp;
  nlp_->set_x(xnext);
  
  // compute the convergence criteria
  double con_crit1 = fabs(xdisp.scalar_product(gcurrent));
  double con_crit2 = maxabs_gradient;
  double con_crit3 = xdisp.maxabs();

  return   (con_crit1 <= convergence_)
            && (con_crit2 <= convergence_)
            && (con_crit3 <= convergence_);
}

RefNLP0
GDIISOpt::nlp()
{
  return nlp_;
}
