//
// matrix3.cc
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

#ifdef __GNUC__
#pragma implementation
#endif

#include <iostream>
#include <iomanip>
#include <math.h>

#include <util/misc/formio.h>
#include <math/scmat/matrix3.h>
#include <math/scmat/vector3.h>

using namespace std;
using namespace sc;

namespace sc {

////////////////////////////////////////////////////////////////////////
// DMatrix3

// Commented out for debugging symmetry class
#if 0
SCMatrix3::SCMatrix3(const RefSCMatrix&x)
{
  if (x.dim().n() != 3) {
      ExEnv::errn() << indent "SCMatrix3::SCMatrix3(RefSCMatrix&): bad length\n";
      abort();
    }
  _v[0] = x.get_element(0);
  _v[1] = x.get_element(1);
  _v[2] = x.get_element(2);
};
#endif 

SCMatrix3::SCMatrix3(double x[9])
{
    _m[0] = x[0];
    _m[1] = x[1];
    _m[2] = x[2];
    _m[3] = x[3];
    _m[4] = x[4];
    _m[5] = x[5];
    _m[6] = x[6];
    _m[7] = x[7];
    _m[8] = x[8];
};

SCMatrix3::SCMatrix3(const SCVector3& c0,
                     const SCVector3& c1,
                     const SCVector3& c2)
{
    operator()(0,0)=c0[0];
    operator()(1,0)=c0[1];
    operator()(2,0)=c0[2];
    operator()(0,1)=c1[0];
    operator()(1,1)=c1[1];
    operator()(2,1)=c1[2];
    operator()(0,2)=c2[0];
    operator()(1,2)=c2[1];
    operator()(2,2)=c2[2];
};

SCMatrix3::SCMatrix3(const SCMatrix3&p)
{
    _m[0] = p[0];
    _m[1] = p[1];
    _m[2] = p[2];
    _m[3] = p[3];
    _m[4] = p[4];
    _m[5] = p[5];
    _m[6] = p[6];
    _m[7] = p[7];
    _m[8] = p[8];
};

SCMatrix3& SCMatrix3::operator=(const SCMatrix3&p)
{
    _m[0] = p[0];
    _m[1] = p[1];
    _m[2] = p[2];
    _m[3] = p[3];
    _m[4] = p[4];
    _m[5] = p[5];
    _m[6] = p[6];
    _m[7] = p[7];
    _m[8] = p[8];

    return *this;
};

// This function builds a rotation matrix that rotates clockwise
// around the given axis
SCMatrix3 rotation_mat(const SCVector3& inaxis, double theta)
{

    // Normalize the rotation axis
    SCVector3 axis=inaxis;
    axis.normalize();
    
    // Calculate the e0-e3  (Following formulae in Goldstein's Classical
    // Mechanics eqn 4-67
    double e0=cos(theta/2.0);
    double e1=axis.x()*sin(theta/2.0);
    double e2=axis.y()*sin(theta/2.0);
    double e3=axis.z()*sin(theta/2.0);
    
    SCMatrix3 result;

    result(0,0)=e0*e0+e1*e1-e2*e2-e3*e3;
    result(1,0)=2.*(e1*e2-e0*e3);
    result(2,0)=2.*(e1*e3+e0*e2);
    result(0,1)=2.*(e1*e2+e0*e3);
    result(1,1)=e0*e0-e1*e1+e2*e2-e3*e3;
    result(2,1)=2.*(e2*e3-e0*e1);
    result(0,2)=2.*(e1*e3-e0*e2);
    result(1,2)=2.*(e2*e3+e0*e1);
    result(2,2)=e0*e0-e1*e1-e2*e2+e3*e3;

    return result;
}

SCMatrix3 rotation_mat(const SCVector3& v1,  const SCVector3& v2, double theta)
{
    return rotation_mat(v1.cross(v2), theta);
}

// This function builds the rotation matrix that will rotate the vector
// ref to the vector target, through an axis that is the cross product
// of the two.
SCMatrix3 rotation_mat(const SCVector3& ref,  const SCVector3& target)
{
    return rotation_mat(target.perp_unit(ref),
                        acos(ref.dot(target)/(ref.norm()*target.norm())));
}

// This function builds a reflection matrix, that reflects the 
// coordinates though a plane perpendicular with unit normal n
// and intersecting the origin
SCMatrix3 reflection_mat(const SCVector3& innormal)
{
    // Normalize the reflection normal
    SCVector3 n = innormal;
    n.normalize();
    SCMatrix3 result;
    
    for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
            result(i,j)=delta(i,j)-2.0*n[i]*n[j];

    return result;
}

SCMatrix3 SCMatrix3::operator*(const SCMatrix3& v) const
{
  SCMatrix3 result;
  for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
      {
          result(i,j) = 0;
          for (int k=0; k<3; k++)
              result(i,j)+=operator()(i,k)*v(k,j);
      }
  return result;
}

SCMatrix3
operator*(double d, const SCMatrix3& v)
{
  SCMatrix3 result;
  for (int i=0; i<9; i++) result[i] = d * v[i];
  return result;
}

SCMatrix3 SCMatrix3::operator*(double d) const
{
  return d*(*this);
}

SCMatrix3 SCMatrix3::operator+(const SCMatrix3&v) const
{
  SCMatrix3 result;
  for (int i=0; i<9; i++) result[i] = _m[i] + v[i];
  return result;
}

SCMatrix3 SCMatrix3::operator-(const SCMatrix3&v) const
{
  SCMatrix3 result;
  for (int i=0; i<9; i++) result[i] = _m[i] - v[i];
  return result;
}


void SCMatrix3::print(ostream& os) const
{
  os << indent << "{"
     << setw(8) << setprecision(5) << operator()(0,0) << " "
     << setw(8) << setprecision(5) << operator()(0,1) << " "
     << setw(8) << setprecision(5) << operator()(0,2) << "}\n";
  os << indent << "{"
     << setw(8) << setprecision(5) << operator()(1,0) << " "
     << setw(8) << setprecision(5) << operator()(1,1) << " "
     << setw(8) << setprecision(5) << operator()(1,2) << "}\n";
  os << indent << "{"
     << setw(8) << setprecision(5) << operator()(2,0) << " "
     << setw(8) << setprecision(5) << operator()(2,1) << " "
     << setw(8) << setprecision(5) << operator()(2,2) << "}\n";
}

}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
