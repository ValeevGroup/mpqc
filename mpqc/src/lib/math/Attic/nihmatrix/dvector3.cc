
#include <stdio.h>
#include "dvector3.h"
#include <math.h>

////////////////////////////////////////////////////////////////////////
// DVector3

DVector3 operator*(double d,const DVector3& v)
{
  DVector3 result;
  for (int i=0; i<3; i++) result[i] = d * v[i];
  return result;
}

DVector3 DVector3::operator*(double d) const
{
  return d*(*this);
}

DVector3 DVector3::operator+(const DVector3&v) const
{
  DVector3 result;
  for (int i=0; i<3; i++) result[i] = _v[i] + v[i];
  return result;
}

DVector3 DVector3::operator-(const DVector3&v) const
{
  DVector3 result;
  for (int i=0; i<3; i++) result[i] = _v[i] - v[i];
  return result;
}

double DVector3::dot(const DVector3&v) const
{
  return _v[0]*v[0] + _v[1]*v[1] + _v[2]*v[2];
}

DVector3 DVector3::cross(const DVector3&v) const
{
  DVector3 result(_v[1]*v._v[2]-_v[2]*v._v[1],
                _v[2]*v._v[0]-_v[0]*v._v[2],
                _v[0]*v._v[1]-_v[1]*v._v[0]);
  return result;
}

void DVector3::rotate(double theta,DVector3&axis)
{
  DVector3 result;
  DVector3 unitaxis = axis;
  unitaxis.normalize();

  // split this into parallel and perpendicular components along axis
  DVector3 parallel = axis * (this->dot(axis) / axis.dot(axis));
  DVector3 perpendicular = (*this) - parallel;

  // form a unit vector perpendicular to parallel and perpendicular
  DVector3 third_axis = axis.cross(perpendicular);
  third_axis.normalize();
  third_axis = third_axis * perpendicular.norm();

  result = parallel + cos(theta) * perpendicular + sin(theta) * third_axis;
  (*this) = result;
}

void DVector3::normalize()
{
  double tmp=0.0;
  int i;
  for (i=0; i<3; i++) tmp += _v[i]*_v[i];
  tmp = 1.0/sqrt(tmp);
  for (i=0; i<3; i++) _v[i] *= tmp;
}

void DVector3::print(FILE*fp) const
{
  fprintf(fp,"{%8.5f %8.5f %8.5f}\n",x(),y(),z());
}


