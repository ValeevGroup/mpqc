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

#include <util/class/scexception.h>
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

SpatialMOPairIter_eq::SpatialMOPairIter_eq(const Ref<MOIndexSpace>& space) :
  SpatialMOPairIter(space,space)
{
  nij_ = ni_*(ni_+1)/2;
  ij_ = 0;

  nij_aa_ = ni_*(ni_-1)/2;
  ij_aa_ = -1;
  nij_ab_ = ni_*nj_;
  ij_ab_ = 0;
  ij_ab_ = 0;
}

SpatialMOPairIter_eq::~SpatialMOPairIter_eq()
{
}

SpatialMOPairIter_neq::SpatialMOPairIter_neq(const Ref<MOIndexSpace>& space1, const Ref<MOIndexSpace>& space2) :
SpatialMOPairIter(space1,space2)
{
  /// If debugging, check if spaces are the same
  if (classdebug() > 0) {
    if (space1 == space2)
      throw ProgrammingError("SpatialMOPairIter_neq::SpatialMOPairIter_neq() -- space1 == space2",__FILE__,__LINE__);
  }
  nij_ = ni_*nj_;
  ij_ = 0;
  IJ_ = 0;
}

SpatialMOPairIter_neq::~SpatialMOPairIter_neq()
{
}

///////

SpinMOPairIter::SpinMOPairIter(const Ref<MOIndexSpace>& space1,
                               const Ref<MOIndexSpace>& space2,
                               const SpinCase2& S) :
  MOPairIter(space1, space2), IJ_(0)
{
  i_eq_j_ = (S!=AlphaBeta && (space1 == space2));
  if (i_eq_j_) {
    if (ni_ <= 0)
      throw ProgrammingError("SpinMOPairIter() initialized with rank-0 space",__FILE__,__LINE__);
    nij_ = (ni_ * (ni_-1))/2;
    ij_ = 0;
    i_ = 1;
    j_ = 0;
  }
  else {
    nij_ = ni_ * nj_;
    ij_ = 0;
    i_ = 0;
    j_ = 0;
  }
}

SpinMOPairIter::~SpinMOPairIter() {}

void
SpinMOPairIter::start(const int first_ij)
{
  IJ_ = 0;
  if (nij_ > 0) {
    ij_ = first_ij%nij_;
    if (i_eq_j_) {
      i_ = (int)floor((sqrt(1.0+8.0*ij_) - 1.0)/2.0) + 1;
      const int i_off = i_*(i_-1)/2;
      j_ = ij_ - i_off;
    }
    else {
      i_ = ij_/nj_;
      j_ = ij_ - i_*nj_;
    }
  }
}

void
SpinMOPairIter::next()
{
  IJ_++;
  ij_++;
  if (ij_ == nij_)
    ij_ = 0;
  if (i_eq_j_) {
    ++j_;
    if (j_ >= i_) {
      i_++;
      j_ = 0;
      if (i_ >= ni_) {
        i_ = 1;
      }
    }
  }
  else {
    ++j_;
    if (j_ == nj_) {
      j_ = 0;
      ++i_;
      if (i_ >= ni_) {
        i_ = 0;
        j_ = 0;
      }
    }
  }
}

SpinMOPairIter::operator int() const { return nij_ > IJ_; }

///////

PureSpinPairIter::PureSpinPairIter(const Ref<MOIndexSpace>& space,
				   const PureSpinCase2& S) :
  MOPairIter(space, space), spin_(S), IJ_(0)
{
  if (spin_ == Triplet) {
    nij_ = (ni_ * (ni_-1))/2;
    ij_ = 0;
    i_ = 1;
    j_ = 0;
  }
  else {
    nij_ = (ni_ * (ni_+1))/2;
    ij_ = 0;
    i_ = 0;
    j_ = 0;
  }
}

PureSpinPairIter::~PureSpinPairIter() {}

void
PureSpinPairIter::start(const int first_ij)
{
  IJ_ = 0;
  if (nij_ > 0) {
    ij_ = first_ij%nij_;
    if (spin_ == Triplet) {
      i_ = (int)floor((sqrt(1.0+8.0*ij_) - 1.0)/2.0) + 1;
      const int i_off = i_*(i_-1)/2;
      j_ = ij_ - i_off;
    }
    else {
      i_ = (int)floor((sqrt(1.0+8.0*ij_) + 1.0)/2.0) - 1;
      const int i_off = i_*(i_+1)/2;
      j_ = ij_ - i_off;
    }
  }
}

void
PureSpinPairIter::next()
{
  IJ_++;
  ij_++;
  if (ij_ == nij_)
    ij_ = 0;
  if (spin_ == Triplet) {
    ++j_;
    if (j_ >= i_) {
      i_++;
      j_ = 0;
      if (i_ >= ni_) {
        i_ = 1;
      }
    }
  }
  else {
    ++j_;
    if (j_ > i_) {
      i_++;
      j_ = 0;
      if (i_ >= ni_) {
        i_ = 0;
      }
    }
  }
}

PureSpinPairIter::operator int() const { return nij_ > IJ_; }

///////

Ref<SpatialMOPairIter>
MOPairIterFactory::mopairiter(const Ref<MOIndexSpace>& space1, const Ref<MOIndexSpace>& space2)
{
  if (space1 == space2)
    return new SpatialMOPairIter_eq(space1);
  else
    return new SpatialMOPairIter_neq(space1, space2);
}

RefSCDimension
MOPairIterFactory::scdim_aa(const Ref<MOIndexSpace>& space1, const Ref<MOIndexSpace>& space2)
{
  if (space1 != space2)
    return scdim_ab(space1,space2);
  else {
    const int n = space1->rank();
    const int npair_aa = n*(n-1)/2;
    std::string name = "Alpha-alpha pair (" + space1->name() + "," + space1->name() + ")";
    return new SCDimension(npair_aa,name.c_str());
  }
}

RefSCDimension
MOPairIterFactory::scdim_ab(const Ref<MOIndexSpace>& space1, const Ref<MOIndexSpace>& space2)
{
  Ref<MOPairIter> piter = mopairiter(space1,space2);
  int npair_ab = space1->rank() * space2->rank();
  std::string name = "Alpha-beta pair (" + space1->name() + "," + space2->name() + ")";
  return new SCDimension(npair_ab,name.c_str());
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
