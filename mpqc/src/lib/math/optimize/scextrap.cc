//
// scextrap.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@ca.sandia.gov>
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

#include <math/optimize/scextrap.h>

SavableState_REF_def(SCExtrapData);

#define CLASSNAME SCExtrapData
#define PARENTS public SavableState
#include <util/class/classia.h>

void *
SCExtrapData::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

SCExtrapData::SCExtrapData()
{
}

SCExtrapData::SCExtrapData(StateIn& s) :
  SavableState(s)
{
}

SCExtrapData::~SCExtrapData()
{
}

void
SCExtrapData::save_data_state(StateOut& s)
{
}

////////////////////////////////////////////////////////////////////////////

SavableState_REF_def(SCExtrapError);

#define CLASSNAME SCExtrapError
#define PARENTS public SavableState
#include <util/class/classia.h>
void *
SCExtrapError::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

SCExtrapError::SCExtrapError()
{
}

SCExtrapError::SCExtrapError(StateIn& s) :
  SavableState(s)
{
}

SCExtrapError::~SCExtrapError()
{
}

void
SCExtrapError::save_data_state(StateOut& s)
{
}

////////////////////////////////////////////////////////////////////////////

SavableState_REF_def(SelfConsistentExtrapolation);

#define CLASSNAME SelfConsistentExtrapolation
#define PARENTS public SavableState
#include <util/class/classia.h>
void *
SelfConsistentExtrapolation::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

SelfConsistentExtrapolation::SelfConsistentExtrapolation()
{
  errorset_ = 0;
  tolerance_ = 1.0e-8;
}

SelfConsistentExtrapolation::SelfConsistentExtrapolation(StateIn& s) :
  SavableState(s)
{
  s.get(error_);
  s.get(errorset_);
  s.get(tolerance_);
}

SelfConsistentExtrapolation::SelfConsistentExtrapolation(
    const RefKeyVal&keyval)
{
  errorset_ = 0;
  tolerance_ = keyval->doublevalue("tolerance");
  if (keyval->error() != KeyVal::OK) tolerance_ = 1.0e-8;
}

SelfConsistentExtrapolation::~SelfConsistentExtrapolation()
{
}

void
SelfConsistentExtrapolation::save_data_state(StateOut& s)
{
  s.put(error_);
  s.put(errorset_);
  s.put(tolerance_);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")
// End:
