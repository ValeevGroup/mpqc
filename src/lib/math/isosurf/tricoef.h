//
// tricoef.h
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

#ifndef _math_isosurf_tricoef_h
#define _math_isosurf_tricoef_h

#ifdef __GNUC__
#pragma interface
#endif

#include <util/ref/ref.h>

namespace sc {

class TriInterpCoefKey {
  private:
    int order_;
    double L2_;
    double L3_;
  public:
    TriInterpCoefKey(int order, double L2, double L3):
      order_(order), L2_(L2), L3_(L3) {}
    int order() const { return order_; }
    double L1() const { return 1.0 - L2_ - L3_; }
    double L2() const { return L2_; }
    double L3() const { return L3_; }
    int cmp(const TriInterpCoefKey&t) const {
        if (order_ < t.order_) return -1;
        if (order_ > t.order_) return 1;
        if (L2_ < t.L2_) return -1;
        if (L2_ > t.L2_) return 1;
        if (L3_ < t.L3_) return -1;
        if (L3_ > t.L3_) return 1;
        return 0;
      }
};

#define TriInterpCoefKeyEQ(k1,k2) ((k1).cmp(k2)==0)
#define TriInterpCoefKeyCMP(k1,k2) ((k1).cmp(k2))

class TriInterpCoef: public RefCount {
    double *coef_;
    double *rderiv_;
    double *sderiv_;
  public:
    TriInterpCoef(const TriInterpCoefKey& key);
    ~TriInterpCoef();
    double& coef(int i, int j, int k) {return coef_[ijk_to_index(i,j,k)];}
    double& rderiv(int i, int j, int k) {return rderiv_[ijk_to_index(i,j,k)];}
    double& sderiv(int i, int j, int k) {return sderiv_[ijk_to_index(i,j,k)];}

    static int
    ijk_to_index(int i, int j, int k)
    {
      int n = i + j + k;
      int ir = n - i;
      return (ir*(ir+1)>>1) + j;
    }

    static int
    order_to_nvertex(int order)
    {
      return ((order+1)*(order+2)>>1);
    }
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
