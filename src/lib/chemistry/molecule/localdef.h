//
// localdef.h
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Edward Seidl <seidl@janed.com>
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

// some inline functions for dealing with 3 dimensional vectors

#ifndef _localdef_h
#define _localdef_h

#include <math.h>

namespace sc {

static const double pi=M_PI;
static const double pih=M_PI_2;
static const double tpi=2.0*pi;

static const double bohr = 0.52917706;

// /////////////////////////////////////////////////////////

static inline void
delta(double u[], const double a[], const double b[])
{
  u[0]=a[0]-b[0];
  u[1]=a[1]-b[1];
  u[2]=a[2]-b[2];
}

// /////////////////////////////////////////////////////////

inline static double
dot(double v[3], double w[3])
{
  return v[0]*w[0] + v[1]*w[1] + v[2]*w[2];
}

// /////////////////////////////////////////////////////////

// returns the distance between two points
static inline double
dist(const double a[], const double b[])
{
  double ab[3];  delta(ab, a, b);
  return sqrt(dot(ab,ab));
}

// /////////////////////////////////////////////////////////

// given sin(x) returns cos(x) 
static inline double
s2(double x)
{
  double tmp = 1.0 - x*x;
  if (tmp < 0.0) tmp = 0.0;
  return sqrt(tmp);
}

// /////////////////////////////////////////////////////////

// returns the dot product for two vectors
static inline double
scalar(const double a[], const double b[])
{
  double x = a[0]*b[0];
  double x1 = a[1]*b[1];
  x += a[2]*b[2];
  return x+x1;
}

// /////////////////////////////////////////////////////////

// given vectors a and b, returns a unit vector directed along the difference
// of the two vectors
static inline void
norm(double u[], const double a[], const double b[])
{
  delta(u,a,b);
  double x = 1.0/sqrt(scalar(u,u));
  u[0] *= x; u[1] *= x; u[2] *= x;
}

// /////////////////////////////////////////////////////////

// given two vectors, returns the normalized cross product of those vectors
static inline void
normal(const double a[], const double b[], double w[])
{
  w[0] = a[1]*b[2]-a[2]*b[1];
  w[1] = a[2]*b[0]-a[0]*b[2];
  w[2] = a[0]*b[1]-a[1]*b[0];
  double x = 1.0/sqrt(scalar(w,w));
  w[0] *= x; w[1] *= x; w[2] *= x;
}

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
