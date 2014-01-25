//
// opt.cc
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

#include <math.h>
#include <deque>

#include <math/optimize/opt.h>
#include <util/keyval/keyval.h>
#include <util/misc/formio.h>
#include <util/misc/regtime.h>
#include <util/state/stateio.h>
#include <util/state/state_bin.h>

using namespace std;
using namespace sc;

/////////////////////////////////////////////////////////////////////////
// Optimize

static ClassDesc Optimize_cd(
  typeid(Optimize),"Optimize",2,"virtual public SavableState",
  0, 0, 0);

Optimize::Optimize() :
  ckpt_(0), print_timings_(0), ckpt_file_("opt_ckpt.dat"),
  max_iterations_(20), max_stepsize_(0.6), conv_(new Convergence())
{
  msg_ = MessageGrp::get_default_messagegrp();
}

Optimize::Optimize(StateIn&s):
  SavableState(s)
{
  msg_ = MessageGrp::get_default_messagegrp();

  s.get(ckpt_,"checkpoint");
  s.get(ckpt_file_);
  s.get(max_iterations_,"max_iterations");
  s.get(max_stepsize_,"max_stepsize");
  if (s.version(::class_desc<Optimize>()) > 1) {
      s.get(print_timings_,"print_timings");
    }
  n_iterations_ = 0;
  conv_ << SavableState::restore_state(s);
  function_ << SavableState::key_restore_state(s,"function");
}

Optimize::Optimize(const Ref<KeyVal>&keyval)
{
  msg_ = MessageGrp::get_default_messagegrp();

  print_timings_ = keyval->booleanvalue("print_timings");
  if (keyval->error() != KeyVal::OK) print_timings_ = 0;
  ckpt_ = keyval->booleanvalue("checkpoint");
  if (keyval->error() != KeyVal::OK) ckpt_ = 0;
  ckpt_file_ = keyval->stringvalue("checkpoint_file");
  if (keyval->error() != KeyVal::OK) {
    ckpt_file_ = "opt_ckpt.dat";
  }

  max_iterations_ = keyval->intvalue("max_iterations", KeyValValueint(20));
  n_iterations_ = 0;

  max_stepsize_ = keyval->doublevalue("max_stepsize");
  if (keyval->error() != KeyVal::OK) max_stepsize_ = 0.6;

  function_ << keyval->describedclassvalue("function");
//  if (function_.null()) {
//      ExEnv::err0() << "Optimize requires a function keyword" << endl;
//      ExEnv::err0() << "which is an object of type Function" << endl;
//      abort();
//    }
// can't assume lineopt's have a function keyword

  conv_ << keyval->describedclassvalue("convergence");
  if (conv_.null()) {
      double convergence = keyval->doublevalue("convergence");
      if (keyval->error() == KeyVal::OK) {
          conv_ = new Convergence(convergence);
        }
    }
  if (conv_.null()) conv_ = new Convergence();
}

Optimize::~Optimize()
{
}

void
Optimize::save_data_state(StateOut&s)
{
  s.put(ckpt_);
  s.put(ckpt_file_);
  s.put(max_iterations_);
  s.put(max_stepsize_);
  s.put(print_timings_);
  SavableState::save_state(conv_.pointer(),s);
  SavableState::save_state(function_.pointer(),s);
}

void
Optimize::init()
{
  n_iterations_ = 0;
}

void
Optimize::set_checkpoint()
{
  ckpt_=1;
}

void
Optimize::set_max_iterations(int mi)
{
  max_iterations_ = mi;
}

void
Optimize::set_checkpoint_file(const char *path)
{
  if (path) {
    ckpt_file_ = path;
    }
  else {
    ckpt_file_.resize(0);
    }
}
  
void
Optimize::set_function(const Ref<Function>& f)
{
  function_ = f;
}

#ifndef OPTSTATEOUT
#define OPTSTATEOUT StateOutBin
#endif

int
Optimize::optimize()
{
  int result=0;
  while((n_iterations_ < max_iterations_) && (!(result = update()))) {
      ++n_iterations_;
      if (ckpt_) {
        OPTSTATEOUT so(ckpt_file_.c_str());
        this->save_state(so);
      }
      if (print_timings_) {
          Timer t;
          t.print();
        }
    }
  return result;
}

void
Optimize::apply_transform(const Ref<NonlinearTransform> &t)
{
}

static inline const char *bp(bool b) { return b?"yes":"no"; };

void
Optimize::print(std::ostream& o) const
{
  o << indent << "Optimize";
  if (::class_desc<Optimize>() != class_desc()) {
      o << " (base class of " << class_desc()->name() << ")";
    }
  o << ":" << std::endl;
  o << incindent;

  o << indent << "print_timings  = " << bp(print_timings_) << std::endl;
  o << indent << "checkpoint     = " << bp(ckpt_) << std::endl;
  o << indent << "max_iterations = " << max_iterations_ << std::endl;
  o << indent << "max_stepsize   = " << max_stepsize_ << std::endl;

  if (conv_.nonnull()) {
      o << indent << "convergence    = " << std::endl;
      o << incindent;
      conv_->print(o);
      o << decindent;
    }
  else {
      o << indent << "convergence    = 0 (no Convergence object given)"
        << std::endl;
    }

  o << decindent;
}

/////////////////////////////////////////////////////////////////////////
// LineOpt

static ClassDesc LineOpt_cd(
  typeid(LineOpt),"LineOpt",1,"public Optimize",
  0, 0, 0);

LineOpt::LineOpt(StateIn&s):
  Optimize(s), SavableState(s)
{
  search_direction_ = matrixkit()->vector(dimension());
  search_direction_.restore(s);
}

LineOpt::LineOpt(const Ref<KeyVal>&keyval)
{
}

LineOpt::LineOpt()
{
}

LineOpt::~LineOpt()
{
}

void
LineOpt::save_data_state(StateOut&s)
{
  search_direction_.save(s);
}

void
LineOpt::init(RefSCVector& direction)
{
  if (function().null()) {
      ExEnv::err0() << "LineOpt requires a function object through" << endl;
      ExEnv::err0() << "constructor or init method" << endl;
      abort();
  }
  search_direction_ = direction.copy();
  initial_x_ = function()->get_x();
  initial_value_ = function()->value();
  initial_grad_ = function()->gradient();
  Optimize::init();
}

void
LineOpt::init(RefSCVector& direction, Ref<Function> function )
{
  set_function(function);
  init(direction);
}

void
LineOpt::apply_transform(const Ref<NonlinearTransform> &t)
{
  if (t.null()) return;
  apply_transform(t);
  t->transform_gradient(search_direction_);
}

void
LineOpt::print(std::ostream&o) const
{
  Optimize::print(o);
}

/////////////////////////////////////////////////////////////////////////
// Backtrack

static ClassDesc Backtrack_cd(
  typeid(Backtrack),"Backtrack",2,"public LineOpt",
  create<Backtrack>, create<Backtrack>, create<Backtrack>);

Backtrack::Backtrack()
  : LineOpt(), backtrack_factor_(0.1),
    decrease_factor_(0.1), force_search_(0)
{
}

Backtrack::Backtrack(const Ref<KeyVal>& keyval) 
  : LineOpt(keyval)
{ 
  backtrack_factor_ = keyval->doublevalue("backtrack_factor",
                                          KeyValValuedouble(0.1));
  decrease_factor_ = keyval->doublevalue("decrease_factor",
                                         KeyValValuedouble(0.1));
  force_search_ = keyval->booleanvalue("force_search",
                                       KeyValValueboolean(false));
}

Backtrack::Backtrack(StateIn&s):
  LineOpt(s), SavableState(s)
{
  if (s.version(::class_desc<Backtrack>()) > 1) {
      s.get(backtrack_factor_);
      s.get(decrease_factor_);
      s.get(force_search_);
    }
}

void
Backtrack::save_data_state(StateOut&s)
{
  s.put(backtrack_factor_);
  s.put(decrease_factor_);
  s.put(force_search_);
}

Backtrack::~Backtrack()
{
}

int
Backtrack::sufficient_decrease(RefSCVector& step) {

  const double ftarget = initial_value_ + decrease_factor_ *
      initial_grad_.scalar_product(step);
  
  RefSCVector xnext = initial_x_ + step;
  function()->set_x(xnext);
  Ref<NonlinearTransform> t = function()->change_coordinates();
  apply_transform(t);

  return function()->value() <= ftarget;
}

int
Backtrack::update() {
  
  deque<double> values;
  int acceptable=0;
  int descent=1;
  int took_step=0;
  int using_step;

  RefSCVector backtrack = -1.0 * backtrack_factor_ * search_direction_;
  RefSCVector step = search_direction_.copy();
  
  // check if sufficient decrease is possible within the required accuracy
  if (function()->desired_value_accuracy() >= fabs(decrease_factor_ *
                                                   initial_grad_.scalar_product(step))) {
    ExEnv::out0() << endl << indent <<
      "Sufficient decrease cannot be achieved for the requested value accuracy. Skipping line search." << endl;
    RefSCVector xnext = initial_x_ + step;
    function()->set_x(xnext);
    Ref<NonlinearTransform> t = function()->change_coordinates();
    apply_transform(t);
    return 1;
  }

  // check if line search is needed
  if( sufficient_decrease(step) ) {
    ExEnv::out0() << endl << indent << 
      "Unscaled initial step yields sufficient decrease." << endl;
    return 1;
  }

  ExEnv::out0() << endl << indent 
    << "Unscaled initial step does not yield a sufficient decrease."
    << endl << indent
    << "Initiating backtracking line search." << endl; 

  // perform a simple backtrack
  values.push_back(function()->value());
  for (int i = 0; i < max_iterations_ && !acceptable && descent; ++i) {

    step = step + backtrack;

    if (sqrt(step.scalar_product(step))
        >= 0.1 * sqrt(search_direction_.scalar_product(search_direction_))) {

      ++took_step;
      if (sufficient_decrease(step)) {
        ExEnv::out0() << endl << indent << "Backtrack " << i + 1
            << " yields a sufficient decrease." << endl;
        acceptable = 1;
        using_step = i + 1;
      }

      else if (values.back() < function()->value()) {
        ExEnv::out0() << endl << indent << "Backtrack " << i + 1
            << " increases value; terminating search." << endl;
        acceptable = 1;
        using_step = i;
      }

      else {
        ExEnv::out0() << endl << indent << "Backtrack " << i + 1
            << " does not yield a sufficient decrease." << endl;
        using_step = i + 1;
      }
      
      values.push_back(function()->value());
    }

    else {
      ExEnv::out0() << indent
          << "Search direction does not appear to be a descent direction;"
          << " terminating search." << endl;
      descent = 0;
    }
  }
	     

  if ( !acceptable && descent ) {
    ExEnv::out0() << indent << 
      "Maximum number of backtrack iterations has been exceeded." << endl;
    acceptable = 1;
  }
 
  for(int i=0; i <= took_step; ++i) {
    if(i==0) ExEnv::out0() << indent << "initial step    " << " value: ";
    else ExEnv::out0() << indent << "backtrack step " << i << " value: ";
    ExEnv::out0() << scprintf("%15.10lf", values.front()) << endl;
    values.pop_front();
  }

  if(descent) { 
    ExEnv::out0() << indent << "Using step " << using_step << endl;
   
    // use next to last step if value went up
    if( using_step != took_step ) {
      function()->set_x( function()->get_x() - backtrack );
      Ref<NonlinearTransform> t = function()->change_coordinates();
      apply_transform(t);
    }
  }
  else {
    function()->set_x( initial_x_ );
    Ref<NonlinearTransform> t = function()->change_coordinates();
    apply_transform(t);
  }

  // returning 0 only if search direction is not descent direction
  return acceptable;
}

void
Backtrack::print(std::ostream&o) const
{
  o << indent
    << "Backtrack:"
    << std::endl
    << incindent
    << indent << "backtrack_factor = " << backtrack_factor_
    << std::endl
    << indent << "decrease_factor  = " << decrease_factor_
    << std::endl
    << indent << "force_search     = " << (force_search_?"yes":"no")
    << std::endl;

  // Parent class parameters are not relevant to Backtrack
  //LineOpt::print(o);

  o << decindent;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
