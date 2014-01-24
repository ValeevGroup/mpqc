//
// pairiter.cc
//
// Copyright (C) 2004 Edward Valeev
//
// Author: Edward Valeev <evaleev@vt.edu>
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

#include <cassert>
#include <util/misc/scexception.h>
#include <math/mmisc/pairiter.h>
#include <math/mmisc/pairiter.impl.h>

using namespace std;
using namespace sc;

inline int max(int a,int b) { return (a > b) ? a : b;}

MOPairIter::MOPairIter(unsigned int n_i, unsigned int n_j) : ni_(n_i), nj_(n_j)
{
  i_ = -1;
  j_ = -1;
}

MOPairIter::~MOPairIter()
{
}

SpatialMOPairIter_eq::SpatialMOPairIter_eq(unsigned int n) :
  SpatialMOPairIter(n, n)
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

SpatialMOPairIter_neq::SpatialMOPairIter_neq(unsigned int n_i, unsigned int n_j) :
SpatialMOPairIter(n_i, n_j)
{
  MPQC_ASSERT(n_i == n_j);
  nij_ = ni_*nj_;
  ij_ = 0;
  IJ_ = 0;
}

SpatialMOPairIter_neq::~SpatialMOPairIter_neq()
{
}

///////

SpinMOPairIter::SpinMOPairIter(unsigned int n_i, unsigned int n_j, bool i_eq_j) :
  MOPairIter(n_i, n_j), IJ_(0), i_eq_j_(i_eq_j)
{
  if (i_eq_j_) {
    MPQC_ASSERT(n_i == n_j);
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
#if 0
PureSpinPairIter::PureSpinPairIter(const Ref<OrbitalSpace>& space,
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
#endif

/////////

namespace sc { namespace fastpairiter {

  template <>
  void MOPairIter<Symm>::init() {
    I_ = 0; J_ = 0; IJ_ = 0;
  }

  template <>
  void MOPairIter<AntiSymm>::init() {
    I_ = 1; J_ = 0; IJ_ = 0;
  }

  template <>
  void MOPairIter<ASymm>::init() {
    I_ = 0; J_ = 0; IJ_ = 0;
  }

  template<> MOPairIter<Symm>::MOPairIter(int nI, int nJ) :
    nI_(nI), nJ_(nJ), IJ_(0), nIJ_(0), I_(0), J_(0) {
    MPQC_ASSERT(nI == nJ);
    nIJ_ = nI_ * (nI_ + 1)/2;
    init();
  }
      
  template<> MOPairIter<AntiSymm>::MOPairIter(int nI, int nJ) :
      nI_(nI), nJ_(nJ), IJ_(0), nIJ_(0), I_(0), J_(0) {
      MPQC_ASSERT(nI == nJ);
      nIJ_ = nI_ * (nI_ - 1)/2;
    }
    
    template<> MOPairIter<ASymm>::MOPairIter(int nI, int nJ) :
      nI_(nI), nJ_(nJ), IJ_(0), nIJ_(0), I_(0), J_(0) {
      nIJ_ = nI_ * nJ_;
    }
    
    template<> void MOPairIter<AntiSymm>::next() {
      ++J_;
      if (J_ == I_) {
        ++I_;
        J_ = 0;
      }
      ++IJ_;
    }
    
    template<> void MOPairIter<ASymm>::next() {
      ++J_;
      if (J_ == nJ_) {
        ++I_;
        J_ = 0;
      }
      ++IJ_;
    }
    
    template<> void MOPairIter<Symm>::next() {
      ++J_;
      if (J_ > I_) {
        ++I_;
        J_ = 0;
      }
      ++IJ_;
    }
    
    template<> int MOPairIter<AntiSymm>::ij(int i, int j) const {
      return i*(i-1)/2 + j;
    }
    
    template<> int MOPairIter<Symm>::ij(int i, int j) const {
      return i*(i+1)/2 + j;
    }
    
    template<> int MOPairIter<ASymm>::ij(int i, int j) const {
      return i*nJ_ + j;
    }
    
    template<> std::size_t npair<AntiSymm>(unsigned int nI, unsigned int nJ) {
      MPQC_ASSERT(nI == nJ);
      return nI * (std::size_t)(nI-1)/2;
    }
    
    template<> std::size_t npair<ASymm>(unsigned int nI, unsigned int nJ) {
      return nI * (std::size_t)nJ;
    }
    
    template<> std::size_t npair<Symm>(unsigned int nI, unsigned int nJ) {
      MPQC_ASSERT(nI == nJ);
      return nI * (std::size_t)(nI+1)/2;
    }

}} // namespace sc::fastpairiter

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
