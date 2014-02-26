//
// vector3.cc
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

#include <iostream>
#include <iomanip>

#include <math/scmat/matrix.h>
#include <math/scmat/vector3.h>
#include <math.h>

#include <util/misc/formio.h>
#include <util/keyval/keyval.h>

using namespace std;
using namespace sc;

namespace sc {

////////////////////////////////////////////////////////////////////////
// DVector3

SCVector3::SCVector3(const Ref<KeyVal>&keyval)
{
  _v[0] = keyval->doublevalue(0);
  _v[1] = keyval->doublevalue(1);
  _v[2] = keyval->doublevalue(2);
}

SCVector3::SCVector3(const RefSCVector&x)
{
  if (x.dim().n() != 3) {
      ExEnv::errn() << indent << "SCVector3::SCVector3(RefSCVEctor&): bad length\n";
      abort();
    }
  _v[0] = x.get_element(0);
  _v[1] = x.get_element(1);
  _v[2] = x.get_element(2);
};

SCVector3
operator*(double d,const SCVector3& v)
{
  SCVector3 result;
  for (int i=0; i<3; i++) result[i] = d * v[i];
  return result;
}

SCVector3 SCVector3::operator*(double d) const
{
  return d*(*this);
}

SCVector3 SCVector3::cross(const SCVector3&v) const
{
  SCVector3 result(_v[1]*v._v[2]-_v[2]*v._v[1],
                _v[2]*v._v[0]-_v[0]*v._v[2],
                _v[0]*v._v[1]-_v[1]*v._v[0]);
  return result;
}

SCVector3 SCVector3::perp_unit(const SCVector3&v) const
{
  // try the cross product
  SCVector3 result(_v[1]*v._v[2]-_v[2]*v._v[1],
                   _v[2]*v._v[0]-_v[0]*v._v[2],
                   _v[0]*v._v[1]-_v[1]*v._v[0]);
  double resultdotresult = result.dot(result);
  if (resultdotresult < 1.e-16) {
      // the cross product is too small to normalize

      // find the largest of this and v
      double dotprodt = this->dot(*this);
      double dotprodv = v.dot(v);
      const SCVector3 *d;
      double dotprodd;
      if (dotprodt < dotprodv) {
          d = &v;
          dotprodd = dotprodv;
        }
      else {
          d = this;
          dotprodd = dotprodt;
        }
      // see if d is big enough
      if (dotprodd < 1.e-16) {
          // choose an arbitrary vector, since the biggest vector is small
          result[0] = 1.0;
          result[1] = 0.0;
          result[2] = 0.0;
          return result;
        }
      else {
          // choose a vector perpendicular to d
          // choose it in one of the planes xy, xz, yz
          // choose the plane to be that which contains the two largest
          // components of d
          double absd[3];
          absd[0] = fabs(d->_v[0]);
          absd[1] = fabs(d->_v[1]);
          absd[2] = fabs(d->_v[2]);
          int axis0, axis1;
          if (absd[0] < absd[1]) {
              axis0 = 1;
              if (absd[0] < absd[2]) {
                  axis1 = 2;
                }
              else {
                  axis1 = 0;
                }
            }
          else {
              axis0 = 0;
              if (absd[1] < absd[2]) {
                  axis1 = 2;
                }
              else {
                  axis1 = 1;
                }
            }
          result[0] = 0.0;
          result[1] = 0.0;
          result[2] = 0.0;
          // do the pi/2 rotation in the plane
          result[axis0] = d->_v[axis1];
          result[axis1] = -d->_v[axis0];
        }
      result.normalize();
      return result;
    }
  else {
      // normalize the cross product and return the result
      result *= 1.0/sqrt(resultdotresult);
      return result;
    }
}

void SCVector3::rotate(double theta,SCVector3&axis)
{
  SCVector3 result;
  SCVector3 unitaxis = axis;
  unitaxis.normalize();

  // split this into parallel and perpendicular components along axis
  SCVector3 parallel = axis * (this->dot(axis) / axis.dot(axis));
  SCVector3 perpendicular = (*this) - parallel;

  // form a unit vector perpendicular to parallel and perpendicular
  SCVector3 third_axis = axis.perp_unit(perpendicular);
  third_axis = third_axis * perpendicular.norm();

  result = parallel + cos(theta) * perpendicular + sin(theta) * third_axis;
  (*this) = result;
}

void SCVector3::normalize()
{
  double tmp=0.0;
  int i;
  for (i=0; i<3; i++) tmp += _v[i]*_v[i];
  tmp = 1.0/sqrt(tmp);
  for (i=0; i<3; i++) _v[i] *= tmp;
}

double
SCVector3::maxabs() const
{
  double result = fabs(_v[0]);
  double tmp;
  if ((tmp = fabs(_v[1])) > result) result = tmp;
  if ((tmp = fabs(_v[2])) > result) result = tmp;
  return result;
}

void
SCVector3::spherical_to_cartesian(SCVector3&cart) const
{
  cart.spherical_coord(theta(), phi(), r());
}

void SCVector3::print(ostream& os) const
{
  os << indent << "{"
     << setw(8) << setprecision(5) << x() << " "
     << setw(8) << setprecision(5) << y() << " "
     << setw(8) << setprecision(5) << z() << "}"
     << endl;
}

ostream &
operator<<(ostream&o, const SCVector3 &v)
{
  o << scprintf("{% 8.5f % 8.5f % 8.5f}", v.x(), v.y(), v.z());
  return o;
}

bool operator==(const SCVector3& a, const SCVector3& b) {
  for(int xyz=0; xyz!=3; ++xyz)
    if (a.elem(xyz) != b.elem(xyz))
      return false;
  return true;
}

}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
