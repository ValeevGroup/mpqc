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

#ifdef __GNUC__
#pragma implementation
#endif

#include <math.h>

#include <math/optimize/opt.h>
#include <util/keyval/keyval.h>
#include <util/misc/formio.h>
#include <util/misc/timer.h>
#include <util/state/state_bin.h>

SavableState_REF_def(Optimize);
#define CLASSNAME Optimize
#define VERSION 2
#define PARENTS virtual_base public SavableState
#include <util/state/statei.h>
#include <util/class/classia.h>

SavableState_REF_def(LineOpt);
#define CLASSNAME LineOpt
#define PARENTS public Optimize
#include <util/state/statei.h>
#include <util/class/classia.h>

/////////////////////////////////////////////////////////////////////////
// Optimize

void *
Optimize::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

Optimize::Optimize() :
  ckpt_(0), ckpt_file(0)
{
}

Optimize::Optimize(StateIn&s):
  SavableState(s)
{
  s.get(ckpt_,"checkpoint");
  s.getstring(ckpt_file);
  s.get(max_iterations_,"max_iterations");
  s.get(max_stepsize_,"max_stepsize");
  if (s.version(static_class_desc()) > 1) {
      s.get(print_timings_,"print_timings");
    }
  n_iterations_ = 0;
  conv_.restore_state(s);
  function_.restore_state(s);
}

Optimize::Optimize(const RefKeyVal&keyval)
{
  print_timings_ = keyval->booleanvalue("print_timings");
  if (keyval->error() != KeyVal::OK) print_timings_ = 0;
  ckpt_ = keyval->booleanvalue("checkpoint");
  if (keyval->error() != KeyVal::OK) ckpt_ = 0;
  ckpt_file = keyval->pcharvalue("checkpoint_file");
  if (keyval->error() != KeyVal::OK) {
    ckpt_file = new char[13];
    strcpy(ckpt_file,"opt_ckpt.dat");
  }

  max_iterations_ = keyval->intvalue("max_iterations");
  if (keyval->error() != KeyVal::OK) max_iterations_ = 10;
  n_iterations_ = 0;

  max_stepsize_ = keyval->doublevalue("max_stepsize");
  if (keyval->error() != KeyVal::OK) max_stepsize_ = 0.6;

  function_ = keyval->describedclassvalue("function");
  if (function_.null()) {
      cerr << node0 << "Optimize requires a function keyword" << endl;
      cerr << node0 << "which is an object of type Function" << endl;
      abort();
    }

  conv_ = keyval->describedclassvalue("convergence");
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
  if (ckpt_file) delete[] ckpt_file;
  ckpt_file=0;
}

void
Optimize::save_data_state(StateOut&s)
{
  s.put(ckpt_);
  s.putstring(ckpt_file);
  s.put(max_iterations_);
  s.put(max_stepsize_);
  s.put(print_timings_);
  conv_.save_state(s);
  function_.save_state(s);
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
  if (ckpt_file) delete[] ckpt_file;
  if (path) {
    ckpt_file = new char[strlen(path)+1];
    strcpy(ckpt_file,path);
  } else
    ckpt_file=0;
}
  
void
Optimize::set_function(const RefFunction& f)
{
  function_ = f;
}

#ifndef OPTSTATEOUT
#define OPTSTATEOUT StateOutBin
#endif

int
Optimize::optimize()
{
  int result;
  while((n_iterations_ < max_iterations_) && (!(result = update()))) {
      n_iterations_++;
      if (ckpt_) {
        OPTSTATEOUT so(ckpt_file);
        this->save_state(so);
      }
      if (print_timings_) {
          tim_print(0);
        }
    }
  return result;
}

void
Optimize::apply_transform(const RefNonlinearTransform &t)
{
}

/////////////////////////////////////////////////////////////////////////
// LineOpt

void *
LineOpt::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = Optimize::_castdown(cd);
  return do_castdowns(casts,cd);
}

LineOpt::LineOpt()
{
}

LineOpt::LineOpt(StateIn&s):
  Optimize(s)
  maybe_SavableState(s)
{
  search_direction_ = matrixkit()->vector(dimension());
  search_direction_.restore(s);
}

LineOpt::LineOpt(const RefKeyVal&keyval):
  Optimize(keyval)
{
}

LineOpt::~LineOpt()
{
}

void
LineOpt::save_data_state(StateOut&s)
{
  Optimize::save_data_state(s);
  search_direction_.save(s);
}

void
LineOpt::set_search_direction(RefSCVector&s)
{
  search_direction_ = s.copy();
}

void
LineOpt::apply_tranform(const RefNonlinearTransform &t)
{
  if (t.null()) return;
  Optimize::apply_transform(t);
  t->transform_gradient(search_direction_);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
