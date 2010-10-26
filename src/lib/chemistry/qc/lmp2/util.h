
/*
 * Copyright 2009 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 *
 * This file is a part of the MPQC LMP2 library.
 *
 * The MPQC LMP2 library is free software: you can redistribute it
 * and/or modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/>.
 *
 */

#ifndef _chemistry_qc_lmp2_util_h
#define _chemistry_qc_lmp2_util_h

#include <stdlib.h>
#include <iostream>

namespace sc {

namespace sma2 {

/// Stores a triplet of data.
template <class T1, class T2, class T3>
struct triplet {
  /// The first member of the triplet.
  T1 a;
  /// The second member of the triplet.
  T2 b;
  /// The third member of the triplet.
  T3 c;

  triplet() {}
  triplet(const sma2::triplet<T1, T2, T3>&t):a(t.a),b(t.b),c(t.c) {}
  triplet(const T1&a_,const T2&b_,const T3&c_):a(a_),b(b_),c(c_) {}
  bool operator == (const sma2::triplet<T1, T2, T3> &t)const {
    return a==t.a && b==t.b && c=t.c;
  }
  bool operator < (const sma2::triplet<T1, T2, T3> &t)const {
    if (a<t.a) return true;
    else if (a>t.a) return false;
    if (b<t.b) return true;
    else if (b>t.b) return false;
    if (c<t.c) return true;
    return false;
  }
    
};

template <class T1, class T2, class T3>
sma2::triplet<T1, T2, T3> make_triplet(const T1&a,const T2&b, const T3&c)
{
  return triplet<T1, T2, T3>(a,b,c);
}

}

}

#endif
