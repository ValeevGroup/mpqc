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

#define CLASSNAME QNewtonOpt
#define VERSION 2
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
  Optimize(keyval)
{
  init();

  update_ = keyval->describedclassvalue("update");
  if (update_.nonnull()) update_->set_inverse();
  lineopt_ = keyval->describedclassvalue("lineopt");
  accuracy_ = keyval->doublevalue("accuracy");
  if (keyval->error() != KeyVal::OK) accuracy_ = 0.0001;
  print_x_ = keyval->booleanvalue("print_x");
  print_hessian_ = keyval->booleanvalue("print_hessian");
  print_gradient_ = keyval->booleanvalue("print_gradient");

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
  update_.restore_state(s);
  s.get(accuracy_);
  s.get(take_newton_step_);
  s.get(maxabs_gradient);
  if (s.version(static_class_desc()) > 1) {
    s.get(print_hessian_);
    s.get(print_x_);
    s.get(print_gradient_);
  }
  else {
    print_hessian_ = 0;
    print_x_ = 0;
    print_gradient_ = 0;
  }
  lineopt_.restore_state(s);
}

QNewtonOpt::~QNewtonOpt()
{
}

void
QNewtonOpt::save_data_state(StateOut&s)
{
  Optimize::save_data_state(s);
  ihessian_.save(s);
  update_.save_state(s);
  s.put(accuracy_);
  s.put(take_newton_step_);
  s.put(maxabs_gradient);
  s.put(print_hessian_);
  s.put(print_x_);
  s.put(print_gradient_);
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

  ExEnv::out().flush();
    
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
        ExEnv::out().unsetf(ios::fixed);
        ExEnv::out() << node0 << indent
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
    ExEnv::out() << node0 << indent
         << scprintf("NOTICE: maxabs_gradient increased from %8.4e to %8.4e",
                     old_maxabs_gradient, maxabs_gradient) << endl;
  }

  if (!take_newton_step_ && lineopt_.nonnull()) {
      // see if the line min step is really needed
//       if (line min really needed) {
//           take_newton_step_ = (lineopt_->update() == 1);
//           // maybe check for convergence and return here
//         }
    }

  // update the hessian
  if (update_.nonnull()) {
      update_->update(ihessian_,function(),xcurrent,gcurrent);
    }

  if (print_hessian_) {
    RefSymmSCMatrix hessian = ihessian_.gi();
    ExEnv::out() << node0 << indent << "hessian = [" << endl;
    ExEnv::out() << incindent;
    int n = hessian.n();
    for (int i=0; i<n; i++) {
      ExEnv::out() << node0 << indent << "[";
      for (int j=0; j<=i; j++) {
        ExEnv::out() << node0 << scprintf(" % 10.6f",double(hessian(i,j)));
      }
      ExEnv::out() << node0 << " ]" << endl;
    }
    ExEnv::out() << decindent;
    ExEnv::out() << node0 << indent << "]" << endl;
  }
  if (print_x_) {
    int n = xcurrent.n();
    ExEnv::out() << node0 << indent << "x = [";
    for (int i=0; i<n; i++) {
      ExEnv::out() << node0 << scprintf(" % 16.12f",double(xcurrent(i)));
    }
    ExEnv::out() << node0 << " ]" << endl;
  }
  if (print_gradient_) {
    int n = gcurrent.n();
    ExEnv::out() << node0 << indent << "gradient = [";
    for (int i=0; i<n; i++) {
      ExEnv::out() << node0 << scprintf(" % 16.12f",double(gcurrent(i)));
    }
    ExEnv::out() << node0 << " ]" << endl;
  }

  // take the step
  RefSCVector xdisp = -1.0*(ihessian_ * gcurrent);
  // scale the displacement vector if it's too large
  double tot = sqrt(xdisp.scalar_product(xdisp));
  if (tot > max_stepsize_) {
    double scal = max_stepsize_/tot;
    ExEnv::out() << node0 << endl << indent
         << scprintf("stepsize of %f is too big, scaling by %f",tot,scal)
         << endl;
    xdisp.scale(scal);
    tot *= scal;
  }

  RefSCVector xnext = xcurrent + xdisp;

  conv_->reset();
  conv_->get_grad(function());
  conv_->get_x(function());
  conv_->set_nextx(xnext);

  // check for convergence before resetting the geometry
  int converged = conv_->converged();
  if (converged)
    return converged;

  ExEnv::out() << node0 << endl << indent
       << scprintf("taking step of size %f", tot) << endl;
  
  function()->set_x(xnext);
  RefNonlinearTransform t = function()->change_coordinates();
  apply_transform(t);

  // do a line min step next time
  if (lineopt_.nonnull()) take_newton_step_ = 0;

  // make the next gradient computed more accurate, since it will
  // be smaller
  accuracy_ = maxabs_gradient * maxabs_gradient_to_next_desired_accuracy;
  
  return converged;
}

void
QNewtonOpt::apply_transform(const RefNonlinearTransform &t)
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
