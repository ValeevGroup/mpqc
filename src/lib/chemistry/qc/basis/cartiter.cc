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

using namespace sc;

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
  delete[] axis_;
}

////////////////////////////////////////////////////////////////////////
// RedundantCartianIter

RedundantCartesianSubIter::RedundantCartesianSubIter(int l)
{
  l_ = l;
  axis_ = new int[l_];
  zloc_ = new int[l_];
  yloc_ = new int[l_];
}

RedundantCartesianSubIter::~RedundantCartesianSubIter()
{
  delete[] axis_;
  delete[] zloc_;
  delete[] yloc_;
}

void
RedundantCartesianSubIter::start(int a, int b, int c)
{
  if (l_ != a + b + c) {
    ExEnv::err0() << indent
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

  int ii=0;
  for (int i=0; i<c; i++,ii++) { axis_[ii] = 2; zloc_[i] = c-i-1; }
  for (int i=0; i<b; i++,ii++) { axis_[ii] = 1; yloc_[i] = b-i-1; }
  for (int i=0; i<a; i++,ii++) axis_[ii] = 0;
}

static bool
advance(int l, int *loc, int n)
{
  int maxloc = l-1;
  for (int i=0; i<n; i++) {
    if (loc[i] < maxloc) {
      loc[i]++;
      for (int j=i-1; j>=0; j--) loc[j] = loc[j+1] + 1;
      return true;
    }
    else {
      maxloc = loc[i]-1;
    }
  }
  return false;
}

// This loops through all unique axis vectors that have a
// given total a, b, and c.  It is done by looping through
// all possible positions for z, then y, leaving x to be
// filled in.
void
RedundantCartesianSubIter::next()
{
  int currentz = 0;
  int currenty = 0;
  int nz = c();
  int ny = b();

  if (!::advance(l(),zloc_,nz)) {
    if (!::advance(l()-nz,yloc_,ny)) {
      done_ = 1;
      return;
    }
    else {
      for (int i=0; i<nz; i++) { zloc_[i] = nz-i-1; }
    }
  }

  int nonz = l()-nz-1;
  for (int i = l()-1; i>=0; i--) {
    if (currentz<nz && zloc_[currentz]==i) {
      axis_[i] = 2;
      currentz++;
    }
    else if (currenty<ny && yloc_[currenty]==nonz) {
      axis_[i] = 1;
      currenty++;
      nonz--;
    }
    else {
      axis_[i] = 0;
      nonz--;
    }
  }
//    for (int i=0; i<3; i++) cout << " " << e_[i];
//    cout << ": ";
//    for (int i=0; i<l(); i++) cout << " " << axis_[i];
//    cout << endl;
//    if (!valid()) {
//      cout << "ERROR: invalid" << endl;
//      cout << "z: ";
//      for (int i=0; i<c(); i++) cout << " " << zloc_[i];
//      cout << endl;
//      cout << "y: ";
//      for (int i=0; i<b(); i++) cout << " " << yloc_[i];
//      cout << endl;
//    }
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
