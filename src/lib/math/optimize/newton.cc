//
// newton.cc
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

#include <math/optimize/newton.h>
#include <util/keyval/keyval.h>
#include <util/misc/formio.h>
#include <util/state/stateio.h>

using namespace std;
using namespace sc;

/////////////////////////////////////////////////////////////////////////
// NewtonOpt

static ClassDesc NewtonOpt_cd(
  typeid(NewtonOpt),"NewtonOpt",1,"public Optimize",
  0, create<NewtonOpt>, create<NewtonOpt>);

NewtonOpt::NewtonOpt(const Ref<KeyVal>&keyval):
  Optimize(keyval)
{
  init();

  accuracy_ = keyval->doublevalue("accuracy",KeyValValuedouble(0.0001));
  print_x_ = keyval->booleanvalue("print_x");
  print_hessian_ = keyval->booleanvalue("print_hessian");
  print_gradient_ = keyval->booleanvalue("print_gradient");
}

NewtonOpt::NewtonOpt(StateIn&s):
  SavableState(s),
  Optimize(s)
{
  s.get(accuracy_);
  s.get(maxabs_gradient);
  s.get(print_hessian_);
  s.get(print_x_);
  s.get(print_gradient_);
}

NewtonOpt::~NewtonOpt()
{
}

void
NewtonOpt::save_data_state(StateOut&s)
{
  Optimize::save_data_state(s);
  s.put(accuracy_);
  s.put(maxabs_gradient);
  s.put(print_hessian_);
  s.put(print_x_);
  s.put(print_gradient_);
}

void
NewtonOpt::init()
{
  Optimize::init();
  maxabs_gradient = -1.0;
}

int
NewtonOpt::update()
{
  // these are good candidates to be input options
  const double maxabs_gradient_to_desired_accuracy = 0.05;
  const double maxabs_gradient_to_next_desired_accuracy = 0.005;
  const double roundoff_error_factor = 1.1;

  // the gradient convergence criterion.
  double old_maxabs_gradient = maxabs_gradient;
  RefSCVector xcurrent;
  RefSCVector gcurrent;
    
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
        ExEnv::out0() << indent
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
    ExEnv::out0() << indent
         << scprintf("NOTICE: maxabs_gradient increased from %8.4e to %8.4e",
                     old_maxabs_gradient, maxabs_gradient) << endl;
  }

  RefSymmSCMatrix hessian = function()->hessian();
  RefSymmSCMatrix ihessian = function()->inverse_hessian(hessian);

  if (print_hessian_) {
    hessian.print("hessian");
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

  // take the step
  RefSCVector xdisp = -1.0*(ihessian * gcurrent);
  // scale the displacement vector if it's too large
  double tot = sqrt(xdisp.scalar_product(xdisp));
  if (tot > max_stepsize_) {
    double scal = max_stepsize_/tot;
    ExEnv::out0() << endl << indent
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

  ExEnv::out0() << endl << indent
       << scprintf("taking step of size %f", tot) << endl;
  
  function()->set_x(xnext);

  // make the next gradient computed more accurate, since it will
  // be smaller
  accuracy_ = maxabs_gradient * maxabs_gradient_to_next_desired_accuracy;
  
  return converged;
}

void
NewtonOpt::print(std::ostream&o) const
{
  o << indent
    << "NewtonOpt:"
    << std::endl
    << incindent
    << indent << "accuracy         = " << accuracy_
    << std::endl
    << indent << "print_x          = " << print_x_
    << std::endl
    << indent << "print_hessian    = " << print_hessian_
    << std::endl
    << indent << "print_gradient   = " << print_gradient_
    << std::endl;

  Optimize::print(o);
  o << decindent;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
