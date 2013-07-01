//
// cartit.h
//
// Copyright (C) 2001 Edward Valeev
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

#ifndef _chemistry_qc_basis_cartit_h
#define _chemistry_qc_basis_cartit_h

#include <chemistry/qc/basis/cartiter.h>

namespace sc {

class CartesianIterCCA : public CartesianIter {

  public:
    CartesianIterCCA(int l) : CartesianIter(l) {}

    void start() {
      bfn_=b_=c_=0;
      a_=l_;
    }

    void next() {
      if (c_ < l_ - a_) {
	b_--;
	c_++;
      }
      else {
	a_--;
	c_ = 0;
	b_ = l_ - a_;
      }
      bfn_++;
    }
    
    operator int() {
      return (a_ >= 0);
    }
};

class RedundantCartesianIterCCA : public RedundantCartesianIter {
  public:
    RedundantCartesianIterCCA(int l) : RedundantCartesianIter(l) {}

    int bfn() {
      int i = a();
      int am = l();
      if (am == i)
	return 0;
      else {
	int j = b();
	int c = am - i;
	return ((((c+1)*c)>>1)+c-j);
      }
    }
};

class RedundantCartesianSubIterCCA : public RedundantCartesianSubIter {

    int bfn_;

  public:
    RedundantCartesianSubIterCCA(int l) : RedundantCartesianSubIter(l) {}

    void start(int aa, int bb, int cc) {
      RedundantCartesianSubIter::start(aa, bb, cc);
      int i = a();
      int am = l();
      if (am == i)
        bfn_ = 0;
      else {
        int j = b();
        int c = am - i;
        bfn_ = ((((c + 1) * c) >> 1) + c - j);
      }
    }

    int bfn() const { return bfn_; }
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
