//
// gdiis.cc
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

#ifdef __GNUC__
#pragma implementation
#endif

#include <math.h>

#include <math/optimize/gdiis.h>
#include <util/keyval/keyval.h>
#include <math/scmat/local.h>
#include <util/misc/formio.h>

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
  
  update_ = keyval->describedclassvalue("update");
  update_->set_inverse();
  
  convergence_ = keyval->doublevalue("convergence");
  if (keyval->error() != KeyVal::OK) convergence_ = 1.0e-6;
  
  accuracy_ = keyval->doublevalue("accuracy");
  if (keyval->error() != KeyVal::OK) accuracy_ = 0.0001;

  RefSymmSCMatrix hessian(dimension(),matrixkit());
  // get a guess hessian from the function
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
  ihessian_ = function()->inverse_hessian(hessian);

  coords_ = new RefSCVector[nsave];
  grad_ = new RefSCVector[nsave];
  error_ = new RefSCVector[nsave];

  for (int i=0; i < nsave; i++) {
    coords_[i] = matrixkit()->vector(dimension()); coords_[i]->assign(0.0);
    grad_[i] = matrixkit()->vector(dimension()); grad_[i]->assign(0.0);
    error_[i] = matrixkit()->vector(dimension()); error_[i]->assign(0.0);
  }
}

GDIISOpt::GDIISOpt(StateIn&s):
  Optimize(s)
  maybe_SavableState(s)
{
  s.get(nsave);
  s.get(diis_iter);
  ihessian_ = matrixkit()->symmmatrix(dimension());
  ihessian_.restore(s);
  update_.restore_state(s);
  s.get(convergence_);
  s.get(accuracy_);
  s.get(maxabs_gradient);
  coords_ = new RefSCVector[nsave];
  grad_ = new RefSCVector[nsave];
  error_ = new RefSCVector[nsave];
  for (int i=0; i < nsave; i++) {
    coords_[i] = matrixkit()->vector(dimension());
    grad_[i] = matrixkit()->vector(dimension());
    error_[i] = matrixkit()->vector(dimension());
    coords_[i].restore(s);
    grad_[i].restore(s);
    error_[i].restore(s);
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
  ihessian_.save(s);
  update_.save_state(s);
  s.put(convergence_);
  s.put(accuracy_);
  s.put(maxabs_gradient);
  for (int i=0; i < nsave; i++) {
    coords_[i].save(s);
    grad_[i].save(s);
    error_[i].save(s);
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
      cout << node0 << indent
           << "NOTICE: function()->actual_gradient_accuracy() > accuracy_:\n"
           << indent
           << scprintf(
             "        function()->actual_gradient_accuracy() = %15.8e",
             function()->actual_gradient_accuracy()) << endl << indent
           << scprintf(
             "                                     accuracy_ = %15.8e",
             accuracy_) << endl;
    }
  } while(!accurate_enough);

  if (old_maxabs_gradient >= 0.0 && old_maxabs_gradient < maxabs_gradient) {
    cout << node0 << indent
         << scprintf("NOTICE: maxabs_gradient increased from %8.4e to %8.4e",
                     old_maxabs_gradient, maxabs_gradient) << endl;
  }

  // make the next gradient computed more accurate, since it will
  // be smaller
  accuracy_ = maxabs_gradient * maxabs_gradient_to_next_desired_accuracy;
  
  // update the hessian
  if (update_.nonnull()) {
    update_->update(ihessian_,function(),xcurrent,gcurrent);
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
  if (diis_iter==1 || maxabs_gradient > 0.05) {
    // just take the Newton-Raphson step first iteration
    RefSCVector xdisp = -1.0*(ihessian_ * gcurrent);
    // try steepest descent
    // RefSCVector xdisp = -1.0*gcurrent;
    
    // scale displacement vector if it's too large
    double tot = sqrt(xdisp.scalar_product(xdisp));
    if (tot > max_stepsize_) {
      double scal = max_stepsize_/tot;
      cout << node0 << endl << indent
           << scprintf("stepsize of %f is too big, scaling by %f",tot,scal)
           << endl;
      xdisp.scale(scal);
      tot *= scal;
    }
                    
    cout << node0 << endl << indent
         << scprintf("taking step of size %f", tot) << endl;
    
    RefSCVector xnext = xcurrent + xdisp;
    
    conv_->reset();
    conv_->get_grad(function());
    conv_->get_x(function());

    function()->set_x(xnext);
  
    conv_->get_nextx(function());

    return conv_->converged();
  }

  // form the error vectors
  for (i=0; i < howmany; i++)
    error_[i] = -1.0*(ihessian_ * grad_[i]);

  // and form the A matrix
  RefSCMatrix A;
  RefSCVector coeff;
  int ntry=0;

  do {
    int num = howmany-ntry;

    RefSCDimension size = new SCDimension(num+1);
    RefSCMatrixKit lkit = new LocalSCMatrixKit;
    A = lkit->matrix(size,size);
    coeff = lkit->vector(size);

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

  } while (fabs(A.solve_lin(coeff)) < 1.0e-12);

  RefSCVector xstar = matrixkit()->vector(dimension());
  RefSCVector delstar = matrixkit()->vector(dimension());

  xstar.assign(0.0);
  delstar.assign(0.0);

  for (i=0,ii=ntry-1; ii < howmany; i++,ii++) {
    RefSCVector tmp = grad_[ii].copy(); tmp.scale(coeff[i]);
    delstar.accumulate(tmp);
    tmp = coords_[ii].copy(); tmp.scale(coeff[i]);
    xstar.accumulate(tmp);
  }

  RefSCVector xdisp = xstar - xcurrent - ihessian_*delstar;
  // scale displacement vector if it's too large
  double tot = sqrt(xdisp.scalar_product(xdisp));
  if (tot > max_stepsize_) {
    double scal = max_stepsize_/tot;
    cout << node0 << endl << indent
         << scprintf("stepsize of %f is too big, scaling by %f",tot,scal)
         << endl;
    xdisp.scale(scal);
    tot *= scal;
  }

  cout << node0 << endl << indent
       << scprintf("taking step of size %f", tot) << endl;
    
  RefSCVector xnext = xcurrent + xdisp;

  conv_->reset();
  conv_->get_grad(function());
  conv_->get_x(function());

  function()->set_x(xnext);
  
  conv_->get_nextx(function());

  return conv_->converged();
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
