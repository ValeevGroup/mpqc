//
// cartitv3.h
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

#ifndef _chemistry_qc_intv3_cartitv3_h
#define _chemistry_qc_intv3_cartitv3_h

#include <chemistry/qc/basis/cartiter.h>

namespace sc {

class CartesianIterV3 : public CartesianIter {
  public:
    CartesianIterV3(int l) : CartesianIter(l) {}

    void start() {
      bfn_=a_=c_=0;
      b_=l_;
    }

    void next() {
      if (c_<l_-a_)
        c_++;
      else {
        c_=0;
        a_++;
      }
      bfn_++;
      b_ = l_-a_-c_;
    }
    
    operator int() {
      return (a_ <= l_);
    }
};

class RedundantCartesianIterV3 : public RedundantCartesianIter {
  public:
    RedundantCartesianIterV3(int l) : RedundantCartesianIter(l) {}

    int bfn() {
      int i = a();
      int j = b();
      int am = l();
      return (((((((am)+1)<<1)-(i))*((i)+1))>>1)-(j)-1);
    }
};

class RedundantCartesianSubIterV3 : public RedundantCartesianSubIter {
  public:
    RedundantCartesianSubIterV3(int l) : RedundantCartesianSubIter(l) {}

    int bfn() {
      int i = a();
      int j = b();
      int am = l();
      return (((((((am)+1)<<1)-(i))*((i)+1))>>1)-(j)-1);
    }
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
