//
// cartiter.cc
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

#include <stdlib.h>
#include <math.h>

#include <util/misc/formio.h>
#include <util/misc/exenv.h>
#include <chemistry/qc/basis/cartiter.h>

////////////////////////////////////////////////////////////////////////
// CartianIter

CartesianIter::CartesianIter(int l) :
  l_(l)
{
}

CartesianIter::~CartesianIter()
{
}

////////////////////////////////////////////////////////////////////////
// RedundantCartianIter

RedundantCartesianIter::RedundantCartesianIter(int l)
{
  l_ = l;
  axis_ = new int[l_];
}

RedundantCartesianIter::~RedundantCartesianIter()
{
  if (axis_) {
    delete[] axis_;
    axis_=0;
  }
}

////////////////////////////////////////////////////////////////////////
// RedundantCartianIter

RedundantCartesianSubIter::RedundantCartesianSubIter(int l)
{
  l_ = l;
  axis_ = new int[l_];
}

RedundantCartesianSubIter::~RedundantCartesianSubIter()
{
  if (axis_) {
    delete[] axis_;
    axis_=0;
  }
}

void
RedundantCartesianSubIter::start(int a, int b, int c)
{
  if (l_ != a + b + c) {
    ExEnv::err() << node0 << indent
         << "RedundantCartesianSubIter::start: bad args\n";
    abort();
  }

  if (l_==0) {
    done_ = 1;
    return;
  } else {
    done_ = 0;
  }

  e_[0] = a;
  e_[1] = b;
  e_[2] = c;

  for (int i=0; i<l_; i++)
    axis_[i] = 0;

  while (!done_ && !valid()) {
    advance();
  }
}

void
RedundantCartesianSubIter::next()
{
  advance();
  while (!done_ && !valid()) {
    advance();
  }
}

void
RedundantCartesianSubIter::advance()
{
  int i;
  for (i=0; i<l_; i++) {
    if (axis_[i] == 2)
      axis_[i] = 0;
    else {
      axis_[i]++;
      return;
    }
  }
  done_ = 1;
}

int
RedundantCartesianSubIter::valid()
{
  int t[3];
  int i;

  for (i=0; i<3; i++)
    t[i] = 0;

  for (i=0; i<l_; i++)
    t[axis_[i]]++;

  return t[0] == e_[0] && t[1] == e_[1] && t[2] == e_[2];
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
