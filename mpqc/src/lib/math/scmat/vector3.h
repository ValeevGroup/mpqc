
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
    SCVector3() {}
    SCVector3(double x[3]);
    SCVector3(double x,double y,double z);
    SCVector3(const SCVector3&p);
    SCVector3(const RefSCVector&);
    void normalize();
    SCVector3 operator*(double) const;
    void operator *= (double m) { _v[0] *= m; _v[1] *= m; _v[2] *= m; }
    SCVector3 operator+(const SCVector3&) const;
    SCVector3 operator-(const SCVector3&) const;
    double dot(const SCVector3&v) const {
        return _v[0]*v[0] + _v[1]*v[1] + _v[2]*v[2]; }
    SCVector3 cross(const SCVector3&) const;
    // returns a unit vector that is perpendicular to the two vectors
    SCVector3 perp_unit(const SCVector3&) const;
    void spherical_coord(double theta, double phi, 
                         double r);
    // this returns the length of the difference vector
    double dist(const SCVector3&) const;
    void rotate(double theta,SCVector3 &v);
    double norm() const { return sqrt(this->dot(*this)); }
    double& elem(int xyz) { return _v[xyz]; }
    const double& elem(int xyz) const { return _v[xyz]; }
    double& operator [] (int i) { return _v[i]; }
    const double& operator [] (int i) const { return _v[i]; }
    double& operator () (int i) { return _v[i]; }
    const double& operator () (int i) const { return _v[i]; }
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
