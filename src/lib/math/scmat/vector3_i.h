//
// vector3_i.h
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

#ifdef INLINE_FUNCTIONS
#define INLINE inline
#else
#define INLINE
#endif

namespace sc {

INLINE void
SCVector3::spherical_coord(double theta, double phi,
                           double r)
{
  double rsin_theta = r*sin(theta);
  _v[0]=rsin_theta*cos(phi);
  _v[1]=rsin_theta*sin(phi);
  _v[2]=r*cos(theta);
}

INLINE  double
SCVector3::dist(const SCVector3 &s) const
{
  double x=_v[0]-s._v[0],y=_v[1]-s._v[1],z=_v[2]-s._v[2];
  return sqrt(x*x + y*y + z*z);
}

}

#undef INLINE

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
