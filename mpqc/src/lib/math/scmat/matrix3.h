#ifndef _math_scmat_matrix3_h
#define _math_scmat_matrix3_h
#ifdef __GNUC__
#pragma interface
#endif

#include <stdio.h>
#include <math.h>

#include "vector3.h"
class RefSCmatrix;

class SCMatrix3
{
  private:
    double _m[9];
  public:
    SCMatrix3() {}
    SCMatrix3(const SCMatrix3&);
#if 0
    SCMatrix3(const RefSCMatrix3&);
#endif
    SCMatrix3(double x[9]);
    SCMatrix3(const SCVector3&p1, const SCVector3&p2, const SCVector3&p3);
    ~SCMatrix3() {}
    SCMatrix3& operator=(const SCMatrix3&);
    SCMatrix3 operator*(double) const;
    SCMatrix3 operator*(const SCMatrix3&) const;
    SCVector3 operator*(const SCVector3&v) const {
        SCVector3 result;
        result._v[0] = _m[0+3*0]*v._v[0]+_m[0+3*1]*v._v[1]+_m[0+3*2]*v._v[2];
        result._v[1] = _m[1+3*0]*v._v[0]+_m[1+3*1]*v._v[1]+_m[1+3*2]*v._v[2];
        result._v[2] = _m[2+3*0]*v._v[0]+_m[2+3*1]*v._v[1]+_m[2+3*2]*v._v[2];
        return result;
      }
    SCMatrix3 operator+(const SCMatrix3&) const;
    SCMatrix3 operator-(const SCMatrix3&) const;
    double& elem(int i, int j) { return _m[i+3*j]; }
    const double& elem(int i, int j) const { return _m[i+3*j]; }
    double& elem(int i) { return _m[i]; }
    const double& elem(int i) const { return _m[i]; }
    double& operator[] (int i) { return _m[i]; }
    const double& operator[] (int i) const { return _m[i]; }
    double& operator() (int i, int j) { return _m[i+3*j]; }
    const double& operator() (int i, int j) const { return _m[i+3*j]; }
    const double* data() const { return _m; }
    void print(FILE*fp=stdout) const;
};
SCMatrix3 operator*(double,const SCMatrix3&);
SCMatrix3 rotation_mat(const SCVector3&, const SCVector3&, double theta);
SCMatrix3 rotation_mat(const SCVector3&, const SCVector3&);
SCMatrix3 rotation_mat(const SCVector3&, double theta);
SCMatrix3 reflection_mat(const SCVector3&);
inline int delta(int i, int j) { return i==j; }

#endif
