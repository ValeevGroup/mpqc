
#ifndef _math_nihmatrix_dvector3_h
#define _math_nihmatrix_dvector3_h

#include <stdio.h>
#include <math.h>
#include <math/topology/point.h>

class DVector3
{
 private:
  double _v[3];
 public:
  inline DVector3() {};
  inline DVector3(double x[3]) { _v[0] = x[0]; _v[1] = x[1]; _v[2] = x[2]; };
  inline DVector3(double x,double y,double z)
    { _v[0] = x; _v[1] = y; _v[2] = z; };
  inline DVector3(const Point3&p)
    { _v[0] = p[0]; _v[1] = p[1]; _v[2] = p[2]; };
  inline DVector3(const DVector3&p)
    { _v[0] = p[0]; _v[1] = p[1]; _v[2] = p[2]; };
  void normalize();
  DVector3 operator*(double) const;
  DVector3 operator+(const DVector3&) const;
  DVector3 operator-(const DVector3&) const;
  double dot(const DVector3&) const;
  DVector3 cross(const DVector3&) const;
  void rotate(double theta,DVector3 &v);
  inline double norm() { return sqrt(this->dot(*this)); };
  inline double& elem(int xyz) { return _v[xyz]; }
  inline const double& elem(int xyz) const { return _v[xyz]; }
  inline double& operator [] (int i) { return _v[i]; }
  inline const double& operator [] (int i) const { return _v[i]; }
  inline double& x() const { return _v[0]; }
  inline double& y() const { return _v[1]; }
  inline double& z() const { return _v[2]; }
  void print(FILE*fp=stdout) const;
  operator Point3() const { Point3 r(_v); return r; };
};
DVector3 operator*(double,const DVector3&);


#endif
