//
// triangle.h
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

#ifndef _math_isosurf_triangle_h
#define _math_isosurf_triangle_h

#ifdef __GNUC__
#pragma interface
#endif

#include <math/isosurf/tricoef.h>
#include <math/isosurf/edge.h>

class Triangle: public VRefCount {
  protected:
    // these break gcc 2.5.8
    //unsigned int _order:5;
    //unsigned int _orientation0:1;
    //unsigned int _orientation1:1;
    //unsigned int _orientation2:1;
    //unsigned char _order;
    //unsigned char _orientation0;
    //unsigned char _orientation1;
    //unsigned char _orientation2;
    unsigned int _order;
    unsigned int _orientation0;
    unsigned int _orientation1;
    unsigned int _orientation2;
    RefEdge _edges[3];
    RefVertex *_vertices;
  public:
    enum {max_order = 10};

    Triangle(const RefEdge& v1, const RefEdge& v2, const RefEdge& v3,
             unsigned int orient0 = 0);
    RefEdge edge(int i) { return _edges[i]; };
    int contains(const RefEdge&) const;
    unsigned int orientation(int i) const
    {
      return i==0?_orientation0:i==1?_orientation1:_orientation2;
    }
    unsigned int orientation(const RefEdge&) const;
    ~Triangle();
    void add_edges(AVLSet<RefEdge>&);
    void add_vertices(AVLSet<RefVertex>&);

    // returns the surface area element
    // 0<=r<=1, 0<=s<=1, 0<=r+s<=1
    // RefVertex is the intepolated vertex (both point and normal)
    void interpolate(const RefTriInterpCoef&,
                     double r,double s,const RefVertex&v, SCVector3& dA);
    void interpolate(double r,double s,const RefVertex&v, SCVector3& dA);
    void interpolate(double r,double s,const RefVertex&v, SCVector3& dA,
                     const RefVolume &vol, double isovalue);

    // returns a corner vertex from the triangle
    // i = 0 is the (0,0) vertex (or L1 = 1, L2 = 0, L3 = 0)
    // i = 1 is the (r=1,s=0) vertex (or L1 = 0, L2 = 1, L3 = 0)
    // i = 2 is the (r=0,s=1) vertex (or L1 = 0, L2 = 0, L3 = 1)
    RefVertex vertex(int i);

    double flat_area();

    // flip the orientation
    void flip();

    unsigned int order() const { return _order; }

    void set_order(int order, const RefVolume&vol,double isovalue);
};

REF_dec(Triangle);

class TriangleIntegrator: public DescribedClass {
#   define CLASSNAME TriangleIntegrator
#   define HAVE_KEYVAL_CTOR
//#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    int _n;
    double* _r;
    double* _s;
    double* _w;
    // precomputed interpolation coefficients for triangles of various orders
    RefTriInterpCoef **coef_; // max_order by _n
  protected:
    void set_r(int i,double r);
    void set_s(int i,double s);
    void set_w(int i,double w);
    void init_coef();
    void clear_coef();
  public:
    TriangleIntegrator(const RefKeyVal&);
    TriangleIntegrator(int n);
    virtual ~TriangleIntegrator();
    inline double w(int i) { return _w[i]; }
    inline double r(int i) { return _r[i]; }
    inline double s(int i) { return _s[i]; }
    inline int n() { return _n; }
    virtual void set_n(int n);
    RefTriInterpCoef coef(int order, int i) { return coef_[order-1][i]; }
};
DescribedClass_REF_dec(TriangleIntegrator);

class GaussTriangleIntegrator: public TriangleIntegrator {
#   define CLASSNAME GaussTriangleIntegrator
#   define HAVE_KEYVAL_CTOR
//#   include <util/state/stated.h>
#   include <util/class/classd.h>
  private:
    void init_rw(int order);
  public:
    GaussTriangleIntegrator(const RefKeyVal&);
    GaussTriangleIntegrator(int order);
    ~GaussTriangleIntegrator();
    void set_n(int n);
};

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
