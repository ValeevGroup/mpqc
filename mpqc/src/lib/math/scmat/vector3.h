
#ifndef _math_scmat_vector3_h
#define _math_scmat_vector3_h
#ifdef __GNUC__
#pragma interface
#endif

#include <stdio.h>
#include <math.h>

class RefSCVector;

class SCVector3
{
  private:
    double _v[3];
  public:
    SCVector3();
    SCVector3(double x[3]);
    SCVector3(double x,double y,double z);
    SCVector3(const SCVector3&p);
    SCVector3(RefSCVector&);
    void normalize();
    SCVector3 operator*(double) const;
    SCVector3 operator+(const SCVector3&) const;
    SCVector3 operator-(const SCVector3&) const;
    double dot(const SCVector3&) const;
    SCVector3 cross(const SCVector3&) const;
    void rotate(double theta,SCVector3 &v);
    double norm() const;
    double& elem(int xyz);
    const double& elem(int xyz) const;
    double& operator [] (int i);
    const double& operator [] (int i) const;
    double& operator () (int i);
    const double& operator () (int i) const;
    const double* data() const;
    double& x();
    double& y();
    double& z();
    const double& x() const;
    const double& y() const;
    const double& z() const;
    void print(FILE*fp=stdout) const;
};
SCVector3 operator*(double,const SCVector3&);

#ifdef INLINE_FUNCTIONS
#include <math/scmat/vector3_i.h>
#endif


#endif
