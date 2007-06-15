//
// bounds.cc
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

#include <cmath>
#include <chemistry/qc/libint2/bounds.h>

using namespace std;
using namespace sc;

/* Compute the partial bound arrays, either Q or R can be computed
 * with appropiate choice of flag. */
Log2Bounds::int_bound_t
Log2Bounds::bound_cast(double value)
{
  static const double loginv = 1.0/log(2.0);
  double tol = pow(2.0,double(int_bound_min));
  int_bound_t res;

  if (value > tol) res = (int_bound_t) std::ceil(std::log(value)*loginv);
  else res = int_bound_min;
  return res;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
