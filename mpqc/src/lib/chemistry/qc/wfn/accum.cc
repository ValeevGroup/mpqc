//
// accum.cc
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

#include <util/state/stateio.h>
#include <chemistry/qc/wfn/accum.h>
#include <chemistry/qc/basis/integral.h>

///////////////////////////////////////////////////////////////////////////
// AccumH

#define CLASSNAME AccumH
#define PARENTS public SavableState
#include <util/class/classia.h>

void *
AccumH::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

AccumH::AccumH()
{
}

AccumH::AccumH(StateIn&s) :
  SavableState(s)
{
  wfn_.restore_state(s);
}

AccumH::AccumH(const RefKeyVal& keyval)
{
  wfn_ = keyval->describedclassvalue("wavefunction");
}

AccumH::~AccumH()
{
}

void
AccumH::save_data_state(StateOut& s)
{
  wfn_.save_state(s);
}

void
AccumH::init(const RefWavefunction& w)
{
  wfn_ = w;
}

void
AccumH::done()
{
  wfn_ = 0;
}

double
AccumH::e()
{
  return 0.0;
}

///////////////////////////////////////////////////////////////////////////
// AccumHNull

#define CLASSNAME AccumHNull
#define PARENTS public AccumH
#define HAVE_CTOR
#define HAVE_STATEIN_CTOR
#define HAVE_KEYVAL_CTOR
#include <util/class/classi.h>

void *
AccumHNull::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = AccumH::_castdown(cd);
  return do_castdowns(casts,cd);
}

AccumHNull::AccumHNull()
{
}

AccumHNull::AccumHNull(StateIn&s) :
  maybe_SavableState(s)
  AccumH(s)
{
}

AccumHNull::AccumHNull(const RefKeyVal& keyval) :
  AccumH(keyval)
{
}

AccumHNull::~AccumHNull()
{
}

void
AccumHNull::save_data_state(StateOut& s)
{
  AccumH::save_data_state(s);
}

void
AccumHNull::accum(const RefSymmSCMatrix& h)
{
}

/////////////////////////////////////////////////////////////////////////////
// SumAccumH

#define CLASSNAME SumAccumH
#define PARENTS public AccumH
#define HAVE_STATEIN_CTOR
#define HAVE_KEYVAL_CTOR
#include <util/class/classi.h>

void *
SumAccumH::_castdown(const ClassDesc* cd)
{
  void *casts[1];
  casts[0] = AccumH::_castdown(cd);
  return do_castdowns(casts, cd);
}

SumAccumH::SumAccumH(StateIn& s) :
  maybe_SavableState(s)
  AccumH(s)
{
  s.get(n_);
  accums_ = new RefAccumH[n_];
  for (int i=0; i < n_; i++)
    accums_[i].restore_state(s);
}

SumAccumH::SumAccumH(const RefKeyVal& keyval) :
  AccumH(keyval)
{
  n_ = keyval->count("accums");
  accums_ = new RefAccumH[n_];
  for (int i=0; i < n_; i++)
    accums_[i] = keyval->describedclassvalue("accums", i);
}

SumAccumH::~SumAccumH()
{
  if (accums_) {
    delete[] accums_;
    accums_=0;
  }
  n_=0;
}

void
SumAccumH::save_data_state(StateOut& s)
{
  AccumH::save_data_state(s);
  s.put(n_);
  for (int i=0; i < n_; i++)
    accums_[i].save_state(s);
}

void
SumAccumH::init(const RefWavefunction& w)
{
  for (int i=0; i < n_; i++)
    accums_[i]->init(w);
}

void
SumAccumH::accum(const RefSymmSCMatrix& h)
{
  for (int i=0; i < n_; i++)
    accums_[i]->accum(h);
}

void
SumAccumH::done()
{
  for (int i=0; i < n_; i++)
    accums_[i]->done();
}

double
SumAccumH::e()
{
  double te = 0.0;

  for (int i=0; i < n_; i++) {
    te += accums_[i]->e();
  }
    
  return te;
}
  
/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
