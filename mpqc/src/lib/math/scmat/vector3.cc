
#include <stdio.h>
#include "matrix.h"
#include "vector3.h"
#include <math.h>

////////////////////////////////////////////////////////////////////////
// DVector3

SCVector3::SCVector3(RefSCVector&x)
{
  if (v.dim().n() != 3) {
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

SCVector3 SCVector3::operator+(const SCVector3&v) const
{
  SCVector3 result;
  for (int i=0; i<3; i++) result[i] = _v[i] + v[i];
  return result;
}

SCVector3 SCVector3::operator-(const SCVector3&v) const
{
  SCVector3 result;
  for (int i=0; i<3; i++) result[i] = _v[i] - v[i];
  return result;
}

double SCVector3::dot(const SCVector3&v) const
{
  return _v[0]*v[0] + _v[1]*v[1] + _v[2]*v[2];
}

SCVector3 SCVector3::cross(const SCVector3&v) const
{
  SCVector3 result(_v[1]*v._v[2]-_v[2]*v._v[1],
                _v[2]*v._v[0]-_v[0]*v._v[2],
                _v[0]*v._v[1]-_v[1]*v._v[0]);
  return result;
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
  SCVector3 third_axis = axis.cross(perpendicular);
  third_axis.normalize();
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

void SCVector3::print(FILE*fp) const
{
  fprintf(fp,"{%8.5f %8.5f %8.5f}\n",x(),y(),z());
}


