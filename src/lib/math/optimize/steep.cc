//
// steep.cc --- implementation of steepest descent
//
// Copyright (C) 1997 Limit Point Systems, Inc.
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

#include <math.h>
#include <float.h>

#include <util/state/stateio.h>
#include <math/optimize/steep.h>
#include <util/keyval/keyval.h>
#include <util/misc/formio.h>

using namespace std;
using namespace sc;

/////////////////////////////////////////////////////////////////////////
// SteepestDescentOpt

static ClassDesc SteepestDescentOpt_cd(
  typeid(SteepestDescentOpt),"SteepestDescentOpt",2,"public Optimize",
  0, create<SteepestDescentOpt>, create<SteepestDescentOpt>);

SteepestDescentOpt::SteepestDescentOpt(const Ref<KeyVal>&keyval):
  Optimize(keyval),
  maxabs_gradient(-1.0)
{
  lineopt_ << keyval->describedclassvalue("lineopt");
  accuracy_ = keyval->doublevalue("accuracy");
  if (keyval->error() != KeyVal::OK) accuracy_ = 0.0001;
  print_x_ = keyval->booleanvalue("print_x");
  print_gradient_ = keyval->booleanvalue("print_gradient");
}

SteepestDescentOpt::SteepestDescentOpt(StateIn&s):
  SavableState(s),
  Optimize(s)
{
  s.get(accuracy_);
  s.get(take_newton_step_);
  s.get(maxabs_gradient);
  s.get(print_x_);
  s.get(print_gradient_);
  lineopt_ << SavableState::restore_state(s);
}

SteepestDescentOpt::~SteepestDescentOpt()
{
}

void
SteepestDescentOpt::save_data_state(StateOut&s)
{
  Optimize::save_data_state(s);
  s.put(accuracy_);
  s.put(take_newton_step_);
  s.put(maxabs_gradient);
  s.put(print_x_);
  s.put(print_gradient_);
  SavableState::save_state(lineopt_.pointer(),s);
}

void
SteepestDescentOpt::init()
{
  Optimize::init();
  take_newton_step_ = 1;
  maxabs_gradient = -1.0;
}

int
SteepestDescentOpt::update()
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

  // make the next gradient computed more accurate, since it will
  // be smaller
  accuracy_ = maxabs_gradient * maxabs_gradient_to_next_desired_accuracy;
  
  if (!take_newton_step_ && lineopt_) {
      // see if the line min step is really needed
//       if (line min really needed) {
//           take_newton_step_ = (lineopt_->update() == 1);
//           // maybe check for convergence and return here
//         }
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
  RefSCVector xdisp = -1.0*gcurrent;

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

  ExEnv::out0() << endl << indent
       << scprintf("taking step of size %f", tot) << endl;
  
  RefSCVector xnext = xcurrent + xdisp;

  conv_->reset();
  conv_->get_grad(function());
  conv_->get_x(function());
  conv_->set_nextx(xnext);

  function()->set_x(xnext);

  // do a line min step next time
  if (lineopt_) take_newton_step_ = 0;

  return conv_->converged();
}

void
SteepestDescentOpt::print(std::ostream&o) const
{
  o << indent
    << "SteepestDescentOpt:"
    << std::endl
    << incindent
    << indent << "accuracy         = " << accuracy_
    << std::endl
    << indent << "print_x          = " << print_x_
    << std::endl
    << indent << "print_gradient   = " << print_gradient_
    << std::endl;

  if (lineopt_.null()) {
    o << indent << "lineopt          = 0 (line optimization will not be performed)"
      << std::endl;
  }
  else {
    o << indent << "lineopt          =" << std::endl;
    o << incindent;
    lineopt_->print(o);
    o << decindent;
  }

  Optimize::print(o);
  o << decindent;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
