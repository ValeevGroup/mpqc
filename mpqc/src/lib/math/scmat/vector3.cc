
#ifdef __GNUC__
#pragma implementation
#endif

#include <stdio.h>
#include "matrix.h"
#include "vector3.h"
#include <math.h>
#include <util/keyval/keyval.h>

////////////////////////////////////////////////////////////////////////
// DVector3

SCVector3::SCVector3(const RefKeyVal&keyval)
{
  _v[0] = keyval->doublevalue(0);
  _v[1] = keyval->doublevalue(1);
  _v[2] = keyval->doublevalue(2);
}

SCVector3::SCVector3(const RefSCVector&x)
{
  if (x.dim().n() != 3) {
      fprintf(stderr,"SCVector3::SCVector3(RefSCVEctor&): bad length\n");
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

void SCVector3::print(FILE*fp) const
{
  fprintf(fp,"{%8.5f %8.5f %8.5f}\n",x(),y(),z());
}


