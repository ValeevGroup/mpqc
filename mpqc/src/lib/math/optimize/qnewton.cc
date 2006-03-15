//
// qnewton.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
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
#include <float.h>

#include <util/state/stateio.h>
#include <math/optimize/qnewton.h>
#include <util/keyval/keyval.h>
#include <util/misc/formio.h>

using namespace std;
using namespace sc;

/////////////////////////////////////////////////////////////////////////
// QNewtonOpt

static ClassDesc QNewtonOpt_cd(
  typeid(QNewtonOpt),"QNewtonOpt",2,"public Optimize",
  0, create<QNewtonOpt>, create<QNewtonOpt>);

QNewtonOpt::QNewtonOpt(const Ref<KeyVal>&keyval):
  Optimize(keyval)
{

  if (function_.null()) {
      ExEnv::err0() << "QNewtonOpt requires a function keyword" << endl;
      abort();
  }

  init();

  update_ << keyval->describedclassvalue("update");
  if (update_.nonnull()) update_->set_inverse();
  lineopt_ << keyval->describedclassvalue("lineopt");
  accuracy_ = keyval->doublevalue("accuracy");
  if (keyval->error() != KeyVal::OK) accuracy_ = 0.0001;
  print_x_ = keyval->booleanvalue("print_x");
  print_hessian_ = keyval->booleanvalue("print_hessian");
  print_gradient_ = keyval->booleanvalue("print_gradient");
  linear_ = keyval->booleanvalue("linear");
  if (keyval->error() != KeyVal::OK) linear_ = 0;
  restrict_ = keyval->booleanvalue("restrict");
  if (keyval->error() != KeyVal::OK) restrict_ = 1;
  dynamic_grad_acc_ = keyval->booleanvalue("dynamic_grad_acc");
  if (keyval->error() != KeyVal::OK) dynamic_grad_acc_ = 1;
  force_search_ = keyval->booleanvalue("force_search");
  if (keyval->error() != KeyVal::OK) force_search_ = 0;
  restart_ = keyval->booleanvalue("restart");
  if (keyval->error() != KeyVal::OK) restart_ = 1;

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
}

QNewtonOpt::QNewtonOpt(StateIn&s):
  SavableState(s),
  Optimize(s)
{
  ihessian_ = matrixkit()->symmmatrix(dimension());
  ihessian_.restore(s);
  update_ << SavableState::restore_state(s);
  s.get(accuracy_);
  s.get(take_newton_step_);
  s.get(maxabs_gradient);
  if (s.version(::class_desc<QNewtonOpt>()) > 1) {
    s.get(print_hessian_);
    s.get(print_x_);
    s.get(print_gradient_);
  }
  else {
    print_hessian_ = 0;
    print_x_ = 0;
    print_gradient_ = 0;
  }
  lineopt_ << SavableState::restore_state(s);
}

QNewtonOpt::~QNewtonOpt()
{
}

void
QNewtonOpt::save_data_state(StateOut&s)
{
  Optimize::save_data_state(s);
  ihessian_.save(s);
  SavableState::save_state(update_.pointer(),s);
  s.put(accuracy_);
  s.put(take_newton_step_);
  s.put(maxabs_gradient);
  s.put(print_hessian_);
  s.put(print_x_);
  s.put(print_gradient_);
  SavableState::save_state(lineopt_.pointer(),s);
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
  const double maxabs_gradient_to_next_desired_accuracy = 0.005;
  const double roundoff_error_factor = 1.1;
  
  // the gradient convergence criterion.
  double old_maxabs_gradient = maxabs_gradient;
  RefSCVector xcurrent;
  RefSCVector gcurrent;

  if( dynamic_grad_acc_ ) {
    // get the next gradient at the required level of accuracy.
    // usually only one pass is needed, unless we happen to find
    // that the accuracy was set too low.
    int accurate_enough;
    do {
        // compute the current point
        function()->set_desired_gradient_accuracy(accuracy_);
        xcurrent = function()->get_x();
        gcurrent = function()->gradient().copy();
      
        // compute the gradient convergence criterion now so i can see if
        // the accuracy needs to be tighter
        maxabs_gradient = gcurrent.maxabs();
        // compute the required accuracy
        accuracy_ = maxabs_gradient * maxabs_gradient_to_desired_accuracy;
      
        if (accuracy_ < DBL_EPSILON) accuracy_ = DBL_EPSILON;

        // The roundoff_error_factor is thrown in to allow for round off making
        // the current gcurrent.maxabs() a bit smaller than the previous,
        // which would make the current required accuracy less than the
        // gradient's actual accuracy and cause everything to be recomputed.
        accurate_enough = (
            function()->actual_gradient_accuracy()
            <= accuracy_*roundoff_error_factor);

        if (!accurate_enough) {
          ExEnv::out0().unsetf(ios::fixed);
          ExEnv::out0() << indent << 
	       "NOTICE: function()->actual_gradient_accuracy() > accuracy_:\n"
               << indent
               << scprintf(
                 "        function()->actual_gradient_accuracy() = %15.8e",
                 function()->actual_gradient_accuracy()) << endl << indent
               << scprintf(
                 "                                     accuracy_ = %15.8e",
                 accuracy_) << endl;
        }
     } while(!accurate_enough);
    // increase accuracy, since the next gradient will be smaller
    accuracy_ = maxabs_gradient * maxabs_gradient_to_next_desired_accuracy;
  }
  else {
    xcurrent = function()->get_x();
    gcurrent = function()->gradient().copy();
  }

  if (old_maxabs_gradient >= 0.0 && old_maxabs_gradient < maxabs_gradient) {
    ExEnv::out0() << indent
	 << scprintf("NOTICE: maxabs_gradient increased from %8.4e to %8.4e",
                     old_maxabs_gradient, maxabs_gradient) << endl;
  }

  // update the hessian
  if(update_.nonnull())
    update_->update(ihessian_,function(),xcurrent,gcurrent);

  conv_->reset();
  conv_->get_grad(function());
  conv_->get_x(function());

  // compute the quadratic step
  RefSCVector xdisp = -1.0*(ihessian_ * gcurrent);
  RefSCVector xnext = xcurrent + xdisp;

  // either do a lineopt or check stepsize
  double tot;
  if(lineopt_.nonnull()) {
    if (dynamic_cast<Backtrack*>(lineopt_.pointer()) != 0) {
      // The Backtrack line search is a special case.

      // perform a search
      double factor;
      if( n_iterations_ == 0 && force_search_ ) 
        factor = lineopt_->set_decrease_factor(1.0);
      lineopt_->init(xdisp,function());
      // reset value acc here so line search "precomputes" are 
      // accurate enough for subsequent gradient evals
      function()->set_desired_value_accuracy(accuracy_/100);
      int acceptable = lineopt_->update();
      if( n_iterations_ == 0 && force_search_ )
        lineopt_->set_decrease_factor( factor );

      if( !acceptable ) {
        if( force_search_ ) factor = lineopt_->set_decrease_factor(1.0);

        // try a new guess hessian
        if( restart_ ) {
          ExEnv::out0() << endl << indent << 
            "Restarting Hessian approximation" << endl;
          RefSymmSCMatrix hessian(dimension(),matrixkit());
          function()->guess_hessian(hessian);
          ihessian_ = function()->inverse_hessian(hessian);
          xdisp = -1.0 * (ihessian_ * gcurrent);
          lineopt_->init(xdisp,function());
          acceptable = lineopt_->update();
        }
      
        // try steepest descent direction
        if( !acceptable ) {
          ExEnv::out0() << endl << indent << 
            "Trying steepest descent direction." << endl;
          xdisp = -1.0 * gcurrent;
          lineopt_->init(xdisp,function());
          acceptable = lineopt_->update();
        }

        // give up and use steepest descent step
        if( !acceptable ) {
          ExEnv::out0() << endl << indent << 
            "Resorting to unscaled steepest descent step." << endl;
          function()->set_x(xcurrent + xdisp);
          Ref<NonlinearTransform> t = function()->change_coordinates();
          apply_transform(t);
        }

        if( force_search_ ) lineopt_->set_decrease_factor( factor );
      }
    }
    else {
      // All line searches other than Backtrack use this
      ExEnv::out0() << indent
                    << "......................................."
                    << endl
                    << indent
                    << "Starting line optimization."
                    << endl;
      lineopt_->init(xdisp,function());
      int nlineopt = 0;
      int maxlineopt = 3;
      for (int ilineopt=0; ilineopt<maxlineopt; ilineopt++) {
        double maxabs_gradient = function()->gradient()->maxabs();

        int converged = lineopt_->update();

        ExEnv::out0() << indent
                      << "Completed line optimization step " << ilineopt+1
                      << (converged?" (converged)":" (not converged)")
                      << endl
                      << indent
                      << "......................................."
                      << endl;

        if (converged) break;

        // Improve accuracy, since we might be able to reuse the next
        // gradient for the next quasi-Newton step.
        if (dynamic_grad_acc_)  {
          accuracy_ = maxabs_gradient*maxabs_gradient_to_next_desired_accuracy;
          function()->set_desired_gradient_accuracy(accuracy_);
        }
      }
    }
    xnext = function()->get_x();
    xdisp = xnext - xcurrent;
    tot = sqrt(xdisp.scalar_product(xdisp));
  }
  else {

    tot = sqrt(xdisp.scalar_product(xdisp));

    if ( tot > max_stepsize_ ) {
      if( restrict_ ) {
	double scal = max_stepsize_/tot;
	ExEnv::out0() << endl << indent <<
	  scprintf("stepsize of %f is too big, scaling by %f",tot,scal)
	  << endl;
	xdisp.scale(scal);
	tot *= scal;
      }
      else {
	ExEnv::out0() << endl << indent <<
	  scprintf("stepsize of %f is too big, but scaling is disabled",tot)
	  << endl;
      }
    }
    xnext = xcurrent + xdisp;
  }

  if (print_hessian_) {
    RefSymmSCMatrix hessian = ihessian_.gi();
    ExEnv::out0() << indent << "hessian = [" << endl;
    ExEnv::out0() << incindent;
    int n = hessian.n();
    for (int i=0; i<n; i++) {
      ExEnv::out0() << indent << "[";
      for (int j=0; j<=i; j++) {
        ExEnv::out0() << scprintf(" % 10.6f",double(hessian(i,j)));
      }
      ExEnv::out0() << " ]" << endl;
    }
    ExEnv::out0() << decindent;
    ExEnv::out0() << indent << "]" << endl;
  }
  if (print_x_) {
    int n = xcurrent.n();
    ExEnv::out0() << indent << "x = [";
    for (int i=0; i<n; i++) {
      ExEnv::out0() << scprintf(" % 16.12f",double(xcurrent(i)));
    }
    ExEnv::out0() << " ]" << endl;
  }
  if (print_gradient_) {
    int n = gcurrent.n();
    ExEnv::out0() << indent << "gradient = [";
    for (int i=0; i<n; i++) {
      ExEnv::out0() << scprintf(" % 16.12f",double(gcurrent(i)));
    }
    ExEnv::out0() << " ]" << endl;
  }

  // check for convergence
  conv_->set_nextx(xnext);
  int converged = conv_->converged();
  if (converged) return converged;

  ExEnv::out0() << indent
		<< scprintf("taking step of size %f", tot) << endl;
  ExEnv::out0() << indent << "Optimization iteration " 
		<< n_iterations_ + 1 << " complete" << endl;
  ExEnv::out0() << indent
    << "//////////////////////////////////////////////////////////////////////"
    << endl;

  if( lineopt_.null() ) {
    function()->set_x(xnext);
    Ref<NonlinearTransform> t = function()->change_coordinates();
    apply_transform(t);
  }

  if( dynamic_grad_acc_ ) 
    function()->set_desired_gradient_accuracy(accuracy_);

  return converged;
}

void
QNewtonOpt::apply_transform(const Ref<NonlinearTransform> &t)
{
  if (t.null()) return;
  Optimize::apply_transform(t);
  if (lineopt_.nonnull()) lineopt_->apply_transform(t);
  if (ihessian_.nonnull()) t->transform_ihessian(ihessian_);
  if (update_.nonnull()) update_->apply_transform(t);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
