//
// tricoef.cc
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

#include <math/isosurf/triangle.h>
#include <math/isosurf/tricoef.h>

using namespace sc;

/////////////////////////////////////////////////////////////////////////
// Utility functions

static inline void
init_coef_deriv(double L, int order, double *Lcoef, double *Lcoefderiv)
{
  int i;
  Lcoef[0] = 1.0;
  Lcoefderiv[0] = 0.0;
  double spacing = 1.0/order;
  for (i=1; i<=order; i++) {
      Lcoef[i] = Lcoef[i-1] * (L - (i-1)*spacing)/(i*spacing);
      Lcoefderiv[i] = Lcoefderiv[i-1] * (L - (i-1)*spacing)/(i*spacing)
                      + Lcoef[i-1]/(i*spacing);
    }
}


/////////////////////////////////////////////////////////////////////////
// The TriInterpCoef Utility Class

TriInterpCoef::TriInterpCoef(const TriInterpCoefKey& key)
{
  int i,j,k;

  int order = key.order();
  double L1 = key.L1();
  double L2 = key.L2();
  double L3 = key.L3();
  int n = order_to_nvertex(order);
  coef_ = new double[n];
  rderiv_ = new double[n];
  sderiv_ = new double[n];

  double L1coef[Triangle::max_order+1];
  double L2coef[Triangle::max_order+1];
  double L3coef[Triangle::max_order+1];

  double L1coefderiv[Triangle::max_order+1];
  double L2coefderiv[Triangle::max_order+1];
  double L3coefderiv[Triangle::max_order+1];

  init_coef_deriv(L1, order, L1coef, L1coefderiv);
  init_coef_deriv(L2, order, L2coef, L2coefderiv);
  init_coef_deriv(L3, order, L3coef, L3coefderiv);

  // the r derivatives
  double L1coef_r[Triangle::max_order+1];
  double L2coef_r[Triangle::max_order+1];
  double L3coef_r[Triangle::max_order+1];

  // the s derivatives
  double L1coef_s[Triangle::max_order+1];
  double L2coef_s[Triangle::max_order+1];
  double L3coef_s[Triangle::max_order+1];

  // convert into r and s derivatives
  for (i=0; i<=order; i++) {
      L1coef_r[i] = -L1coefderiv[i];
      L1coef_s[i] = -L1coefderiv[i];
      L2coef_r[i] =  L2coefderiv[i];
      L2coef_s[i] =  0.0;
      L3coef_r[i] =  0.0;
      L3coef_s[i] =  L3coefderiv[i];
    }

  for (i=0; i<=order; i++) {
      for (j=0; j <= order-i; j++) {
          k = order - i - j;
          coef(i,j,k) = L1coef[i]*L2coef[j]*L3coef[k];
          sderiv(i,j,k) = L1coef_s[i]*L2coef[j]*L3coef[k]
                         +L1coef[i]*L2coef_s[j]*L3coef[k]
                         +L1coef[i]*L2coef[j]*L3coef_s[k];
          rderiv(i,j,k) = L1coef_r[i]*L2coef[j]*L3coef[k]
                         +L1coef[i]*L2coef_r[j]*L3coef[k]
                         +L1coef[i]*L2coef[j]*L3coef_r[k];
        }
    }
}

TriInterpCoef::~TriInterpCoef()
{
  delete[] coef_;
  delete[] rderiv_;
  delete[] sderiv_;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
