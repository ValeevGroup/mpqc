//
// pairiter.cc
//
// Copyright (C) 2004 Edward Valeev
//
// Author: Edward Valeev <edward.valeev@chemistry.gatech.edu>
// Maintainer: EV
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

#ifdef __GNUG__
#pragma implementation
#endif

#include <chemistry/qc/mbptr12/pairiter.h>

using namespace std;
using namespace sc;

inline int max(int a,int b) { return (a > b) ? a : b;}

MOPairIter::MOPairIter(const Ref<MOIndexSpace>& space_i, const Ref<MOIndexSpace>& space_j)
{
  i_eq_j_ = (space_i == space_j);
  ni_ = space_i->rank();
  nj_ = space_j->rank();
  
  i_ = -1;
  j_ = -1;
}

MOPairIter::~MOPairIter()
{
}

MOPairIter_SD::MOPairIter_SD(const Ref<MOIndexSpace>& space) :
  MOPairIter(space,space)
{
  nij_ = ni_*(ni_+1)/2;
  ij_ = 0;

  nij_aa_ = ni_*(ni_-1)/2;
  ij_aa_ = -1;
  nij_ab_ = ni_*nj_;
  ij_ab_ = 0;
  ij_ab_ = 0;
}

MOPairIter_SD::~MOPairIter_SD()
{
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
