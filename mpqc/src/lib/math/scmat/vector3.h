
#ifndef _math_scmat_vector3_h
#define _math_scmat_vector3_h
#ifdef __GNUC__
#pragma interface
#endif

#include <stdio.h>
#include <math.h>

class RefSCVector;
class SCMatrix3;

class SCVector3
{
    friend class SCMatrix3;
  private:
    double _v[3];
  public:
    SCVector3() {}
    SCVector3(double p[3]) {
        _v[0] = p[0]; _v[1] = p[1]; _v[2] = p[2];
      }
    SCVector3(double x,double y,double z) {
        _v[0] = x; _v[1] = y; _v[2] = z;
      }
    SCVector3(const SCVector3&p) {
        _v[0] = p._v[0]; _v[1] = p._v[1]; _v[2] = p._v[2];
      }
    SCVector3(const RefSCVector&);
    void normalize();
    SCVector3 operator*(double) const;
    void operator -= (const SCVector3& v) {
        _v[0] -= v._v[0];
        _v[1] -= v._v[1];
        _v[2] -= v._v[2];
      }
    void operator += (const SCVector3& v) {
        _v[0] += v._v[0];
        _v[1] += v._v[1];
        _v[2] += v._v[2];
      }
    void operator *= (double m) { _v[0] *= m; _v[1] *= m; _v[2] *= m; }
    SCVector3 operator+(const SCVector3&v) const {
        SCVector3 result;
        result._v[0] = _v[0] + v._v[0];
        result._v[1] = _v[1] + v._v[1];
        result._v[2] = _v[2] + v._v[2];
        return result;
      }
    SCVector3 operator-(const SCVector3&v) const {
        SCVector3 result;
        result._v[0] = _v[0] - v._v[0];
        result._v[1] = _v[1] - v._v[1];
        result._v[2] = _v[2] - v._v[2];
        return result;
      }
    double dot(const SCVector3&v) const {
        return _v[0]*v._v[0] + _v[1]*v._v[1] + _v[2]*v._v[2]; }
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
    const double* data() const { return _v; }
    double& x() { return _v[0]; }
    double& y() { return _v[1]; }
    double& z() { return _v[2]; }
    const double& x() const { return _v[0]; }
    const double& y() const { return _v[1]; }
    const double& z() const { return _v[2]; }
    void print(FILE*fp=stdout) const;
};
SCVector3 operator*(double,const SCVector3&);

#ifdef INLINE_FUNCTIONS
#include <math/scmat/vector3_i.h>
#endif


#endif
