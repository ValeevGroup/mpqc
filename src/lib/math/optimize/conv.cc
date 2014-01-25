//
// conv.cc
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

#include <util/misc/formio.h>
#include <util/keyval/keyval.h>
#include <util/state/stateio.h>
#include <math/optimize/conv.h>

using namespace std;
using namespace sc;

/////////////////////////////////////////////////////////////////////////
// Convergence

static ClassDesc Convergence_cd(
  typeid(Convergence),"Convergence",1,"virtual public SavableState",
  0, create<Convergence>, create<Convergence>);

Convergence::Convergence()
{
  set_defaults();
}

Convergence::Convergence(double tolerance)
{
  set_defaults();
  max_disp_ = tolerance;
  max_grad_ = tolerance;
  rms_disp_ = tolerance;
  rms_grad_ = tolerance;
  graddisp_ = tolerance;
}

Convergence::Convergence(StateIn&s):
  SavableState(s)
{
  s.get(use_max_disp_);
  s.get(use_max_grad_);
  s.get(use_rms_disp_);
  s.get(use_rms_grad_);
  s.get(use_graddisp_);
  s.get(max_disp_);
  s.get(max_grad_);
  s.get(rms_disp_);
  s.get(rms_grad_);
  s.get(graddisp_);
}

Convergence::Convergence(const Ref<KeyVal>&keyval)
{
  use_max_disp_ = keyval->exists("max_disp");
  use_max_grad_ = keyval->exists("max_grad");
  use_rms_disp_ = keyval->exists("rms_disp");
  use_rms_grad_ = keyval->exists("rms_grad");
  use_graddisp_ = keyval->exists("graddisp");
  if (use_max_disp_) max_disp_ = keyval->doublevalue("max_disp");
  if (use_max_grad_) max_grad_ = keyval->doublevalue("max_grad");
  if (use_rms_disp_) rms_disp_ = keyval->doublevalue("rms_disp");
  if (use_rms_grad_) rms_grad_ = keyval->doublevalue("rms_grad");
  if (use_graddisp_) graddisp_ = keyval->doublevalue("graddisp");

  if (!use_max_disp_ && !use_max_grad_
      && !use_rms_disp_ && !use_rms_grad_
      && !use_graddisp_) {
      set_defaults();
    }
}

Convergence::~Convergence()
{
}

void
Convergence::save_data_state(StateOut&s)
{
  s.put(use_max_disp_);
  s.put(use_max_grad_);
  s.put(use_rms_disp_);
  s.put(use_rms_grad_);
  s.put(use_graddisp_);
  s.put(max_disp_);
  s.put(max_grad_);
  s.put(rms_disp_);
  s.put(rms_grad_);
  s.put(graddisp_);
}

void
Convergence::set_defaults()
{
  use_max_disp_ = 0;
  use_max_grad_ = 1;
  use_rms_disp_ = 0;
  use_rms_grad_ = 1;
  use_graddisp_ = 0;
  max_grad_ = 4.0e-6;
  rms_grad_ = 1.0e-6;
}

void
Convergence::get_x(const Ref<Function> &f)
{
  x_ = f->get_x();
}

void
Convergence::set_nextx(const RefSCVector &x)
{
  nextx_ = x->copy();
}

void
Convergence::get_grad(const Ref<Function> &f)
{
  grad_ = f->gradient();
}

int
Convergence::converged()
{
  int fail = 0;
  int pass = 0;

  RefSCVector disp;
  if (x_.nonnull() && nextx_.nonnull()) disp = nextx_ - x_;

  ExEnv::out0() << endl;
  
  if (use_max_grad_ && grad_.nonnull()) {
      check_conv("Max Gradient     ", grad_.maxabs(), max_grad_, pass, fail);
    }
  if (use_rms_grad_ && grad_.nonnull()) {
      check_conv("RMS Gradient     ",
                 sqrt(grad_.scalar_product(grad_)/grad_.n()),
                 rms_grad_, pass, fail);
    }
  if (use_max_disp_ && disp.nonnull()) {
      check_conv("Max Displacement ", disp.maxabs(), max_disp_, pass, fail);
    }
  if (use_rms_disp_ && disp.nonnull()) {
      check_conv("RMS Displacement ",
                 sqrt(disp.scalar_product(disp)/disp.n()),
                 rms_disp_, pass, fail);
    }
  if (use_graddisp_ && disp.nonnull() && grad_.nonnull()) {
      check_conv("Gradient*Displace", fabs(disp.scalar_product(grad_)),
                 graddisp_, pass, fail);
    }
  if (fail + pass == 0) {
      ExEnv::errn() << "ERROR: Convergence::converged: no applicable convergence tests"
           << endl;
      abort();
    }
  if (!fail) {
      ExEnv::out0() << endl
           << indent << "All convergence criteria have been met."
           << endl;
    }
  return !fail;
}

void
Convergence::check_conv(const char *heading,
                        double val, double bound,
                        int &pass, int &fail)
{
  int converged = val <= bound;
  ExEnv::out0() << indent << heading << ": "
       << scprintf("%14.10f ", val)
       << scprintf("%14.10f  ", bound)
       << (converged?"yes":"no")
       << endl;
  if (converged) pass++;
  else fail++;
}

void
Convergence::reset()
{
  grad_ = 0;
  x_ = 0;
  nextx_ = 0;
}

void
Convergence::print(std::ostream&o) const
{
  o << indent << "Convergence";
  if (::class_desc<Convergence>() != class_desc()) {
      o << " (base class of " << class_desc()->name() << ")";
    }
  o << ":" << std::endl;
  o << incindent;
  o << indent << "The following criteria must be simultaneously satified:"
    << std::endl;
  if (use_max_disp_) {
      o << indent << "max_disp         = " << max_disp_ << std::endl;
    }
  if (use_max_grad_) {
      o << indent << "max_grad         = " << max_grad_ << std::endl;
    }
  if (use_rms_disp_) {
      o << indent << "rms_disp         = " << rms_disp_ << std::endl;
    }
  if (use_rms_grad_) {
      o << indent << "rms_grad         = " << rms_grad_ << std::endl;
    }
  if (use_graddisp_) {
      o << indent << "graddisp         = " << graddisp_ << std::endl;
    }
  o << decindent;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
