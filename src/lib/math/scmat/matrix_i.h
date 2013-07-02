//
// matrix_i.h
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

#ifndef _math_scmat_matrix_i_h
#define _math_scmat_matrix_i_h

#include <math/scmat/matrix.h>

// These are the inline candidates for the members defined in matrix.h.

#ifdef INLINE_FUNCTIONS
#define INLINE inline
#else
#define INLINE
#endif

namespace sc {

// /////////////////////////////////////////////////////////////////////////
// SCMatrixdouble inline candidates

INLINE
SCMatrixdouble::SCMatrixdouble(SCMatrix*a,int b,int c):
  matrix(a),i(b),j(c)
{
}
INLINE
SCMatrixdouble::~SCMatrixdouble()
{
}
INLINE double
SCMatrixdouble::operator=(double a)
{
  matrix.set_element(i,j,a);
  return a;
}
INLINE double
SCMatrixdouble::operator=(const SCMatrixdouble& md)
{
  double a = md.val();
  matrix.set_element(i,j,a);
  return a;
}
INLINE
SCMatrixdouble::operator double() const
{
  return matrix.get_element(i,j);
}
INLINE double
SCMatrixdouble::val() const
{
  return matrix.get_element(i,j);
}

// /////////////////////////////////////////////////////////////////////////
// SymmSCMatrixdouble inline candidates

INLINE
SymmSCMatrixdouble::SymmSCMatrixdouble(SymmSCMatrix*a,int b,int c):
  matrix(a),i(b),j(c)
{
}
INLINE
SymmSCMatrixdouble::~SymmSCMatrixdouble()
{
}
INLINE double
SymmSCMatrixdouble::operator=(double a)
{
  matrix.set_element(i,j,a);
  return a;
}
INLINE double
SymmSCMatrixdouble::operator=(const SymmSCMatrixdouble& md)
{
  double a = md.val();
  matrix.set_element(i,j,a);
  return a;
}
INLINE
SymmSCMatrixdouble::operator double()
{
  return matrix.get_element(i,j);
}
INLINE double
SymmSCMatrixdouble::val() const
{
  return matrix.get_element(i,j);
}

// /////////////////////////////////////////////////////////////////////////
// DiagSCMatrixdouble inline candidates

INLINE
DiagSCMatrixdouble::DiagSCMatrixdouble(DiagSCMatrix*a,int b,int c):
  matrix(a),i(b),j(c)
{
}
INLINE
DiagSCMatrixdouble::~DiagSCMatrixdouble()
{
}
INLINE double
DiagSCMatrixdouble::operator=(double a)
{
  matrix.set_element(i,a);
  return a;
}
INLINE double
DiagSCMatrixdouble::operator=(const DiagSCMatrixdouble& md)
{
  double a = md.val();
  matrix.set_element(i,a);
  return a;
}
INLINE
DiagSCMatrixdouble::operator double()
{
  return matrix.get_element(i);
}
INLINE double
DiagSCMatrixdouble::val() const
{
  return matrix.get_element(i);
}

// /////////////////////////////////////////////////////////////////////////
// SCVectordouble inline candidates

INLINE
SCVectordouble::SCVectordouble(SCVector*a,int b):
  vector(a),i(b)
{
}
INLINE
SCVectordouble::~SCVectordouble()
{
}
INLINE double
SCVectordouble::operator=(double a)
{
  vector.set_element(i,a);
  return a;
}
INLINE double
SCVectordouble::operator=(const SCVectordouble& vd)
{
  double a = vd.val();
  vector.set_element(i,a);
  return a;
}
INLINE
SCVectordouble::operator double()
{
  return vector.get_element(i);
}
INLINE double
SCVectordouble::val() const
{
  return vector.get_element(i);
}

}

#undef INLINE

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
