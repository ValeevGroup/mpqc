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

#include <util/state/stateio.h>
#include <chemistry/qc/wfn/accum.h>
#include <chemistry/qc/basis/integral.h>

using namespace sc;

///////////////////////////////////////////////////////////////////////////
// AccumH

static ClassDesc AccumH_cd(
  typeid(AccumH),"AccumH",1,"public SavableState",
  0, 0, 0);

AccumH::AccumH()
{
}

AccumH::AccumH(StateIn&s) :
  SavableState(s)
{
  wfn_ << SavableState::restore_state(s);
}

AccumH::AccumH(const Ref<KeyVal>& keyval)
{
  wfn_ << keyval->describedclassvalue("wavefunction");
}

AccumH::~AccumH()
{
}

void
AccumH::save_data_state(StateOut& s)
{
  SavableState::save_state(wfn_.pointer(),s);
}

void
AccumH::init(const Ref<Wavefunction>& w)
{
  wfn_ = w;
}

void
AccumH::done()
{
  wfn_ = 0;
}

void
AccumH::print_summary()
{
}

double
AccumH::e()
{
  return 0.0;
}

///////////////////////////////////////////////////////////////////////////
// AccumHNull

static ClassDesc AccumHNull_cd(
  typeid(AccumHNull),"AccumHNull",1,"public AccumH",
  create<AccumHNull>, create<AccumHNull>, create<AccumHNull>);

AccumHNull::AccumHNull()
{
}

AccumHNull::AccumHNull(StateIn&s) :
  SavableState(s),
  AccumH(s)
{
}

AccumHNull::AccumHNull(const Ref<KeyVal>& keyval) :
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

static ClassDesc SumAccumH_cd(
  typeid(SumAccumH),"SumAccumH",1,"public AccumH",
  0, create<SumAccumH>, create<SumAccumH>);

SumAccumH::SumAccumH(StateIn& s) :
  SavableState(s),
  AccumH(s)
{
  s.get(n_);
  accums_ = new Ref<AccumH>[n_];
  for (int i=0; i < n_; i++)
    accums_[i] << SavableState::restore_state(s);
}

SumAccumH::SumAccumH(const Ref<KeyVal>& keyval) :
  AccumH(keyval)
{
  n_ = keyval->count("accums");
  accums_ = new Ref<AccumH>[n_];
  for (int i=0; i < n_; i++)
    accums_[i] << keyval->describedclassvalue("accums", i);
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
    SavableState::save_state(accums_[i].pointer(),s);
}

void
SumAccumH::init(const Ref<Wavefunction>& w)
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
