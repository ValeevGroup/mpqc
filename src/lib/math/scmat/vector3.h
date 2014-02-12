//
// vector3.h
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

#ifndef _math_scmat_vector3_h
#define _math_scmat_vector3_h

#include <iostream>
#include <math.h>

#include <util/misc/exenv.h>
#include <util/keyval/keyval.h>

namespace sc {

class RefSCVector;
class SCMatrix3;

/// a 3-element version of SCVector
class SCVector3
{
    friend class SCMatrix3;
  private:
    double _v[3];
  public:
    SCVector3() {}
    SCVector3(const double p[3]) {
        _v[0] = p[0]; _v[1] = p[1]; _v[2] = p[2];
      }
    SCVector3(double d) { _v[0] = d; _v[1] = d; _v[2] = d; }
    SCVector3(double x,double y,double z) {
        _v[0] = x; _v[1] = y; _v[2] = z;
      }
    SCVector3(const SCVector3&p) {
        _v[0] = p._v[0]; _v[1] = p._v[1]; _v[2] = p._v[2];
      }
    SCVector3(const RefSCVector&);
    SCVector3(const Ref<KeyVal>&);
    ~SCVector3() {}
    void normalize();
    SCVector3 operator -() { return SCVector3(-_v[0],-_v[1],-_v[2]); }
    SCVector3 operator*(double) const;
    void operator = (const double *x) {
        _v[0] = x[0];
        _v[1] = x[1];
        _v[2] = x[2];
      }
    void operator = (const SCVector3& x) {
        _v[0] = x._v[0];
        _v[1] = x._v[1];
        _v[2] = x._v[2];
      }
    void operator = (double d) { _v[0] = d; _v[1] = d; _v[2] = d; }
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
    void spherical_to_cartesian(SCVector3&cart) const;
    double maxabs() const;
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
    double* data() { return _v; }
    double& x() { return _v[0]; }
    double& y() { return _v[1]; }
    double& z() { return _v[2]; }
    const double& x() const { return _v[0]; }
    const double& y() const { return _v[1]; }
    const double& z() const { return _v[2]; }
    double& r() { return _v[0]; }
    double& theta() { return _v[1]; }
    double& phi() { return _v[2]; }
    const double& r() const { return _v[0]; }
    const double& theta() const { return _v[1]; }
    const double& phi() const { return _v[2]; }
    void print(std::ostream& =ExEnv::out0()) const;
};
SCVector3 operator*(double,const SCVector3&);
std::ostream &operator<<(std::ostream&, const SCVector3 &);

/// @return true if \c a and \c b are \b exactly identical
bool operator==(const SCVector3& a, const SCVector3& b);
/// @return false if \c a and \c b are \b exactly identical
inline bool operator!=(const SCVector3& a, const SCVector3& b) { return not (a==b); }
}

#ifdef INLINE_FUNCTIONS
#include <math/scmat/vector3_i.h>
#endif

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
