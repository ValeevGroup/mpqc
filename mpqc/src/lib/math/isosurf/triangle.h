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

namespace sc {

class Triangle: public RefCount {
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
    Ref<Edge> _edges[3];
    Ref<Vertex> *_vertices;
  public:
    enum {max_order = 10};

    Triangle(const Ref<Edge>& v1, const Ref<Edge>& v2, const Ref<Edge>& v3,
             unsigned int orient0 = 0);
    Ref<Edge> edge(int i) { return _edges[i]; };
    int contains(const Ref<Edge>&) const;
    unsigned int orientation(int i) const
    {
      return i==0?_orientation0:i==1?_orientation1:_orientation2;
    }
    unsigned int orientation(const Ref<Edge>&) const;
    ~Triangle();
    void add_edges(std::set<Ref<Edge> >&);
    void add_vertices(std::set<Ref<Vertex> >&);

    // returns the surface area element
    // 0<=r<=1, 0<=s<=1, 0<=r+s<=1
    // Ref<Vertex> is the intepolated vertex (both point and normal)
    void interpolate(const Ref<TriInterpCoef>&,
                     double r,double s,const Ref<Vertex>&v, SCVector3& dA);
    void interpolate(double r,double s,const Ref<Vertex>&v, SCVector3& dA);
    void interpolate(double r,double s,const Ref<Vertex>&v, SCVector3& dA,
                     const Ref<Volume> &vol, double isovalue);

    // returns a corner vertex from the triangle
    // i = 0 is the (0,0) vertex (or L1 = 1, L2 = 0, L3 = 0)
    // i = 1 is the (r=1,s=0) vertex (or L1 = 0, L2 = 1, L3 = 0)
    // i = 2 is the (r=0,s=1) vertex (or L1 = 0, L2 = 0, L3 = 1)
    Ref<Vertex> vertex(int i);

    double flat_area();

    // flip the orientation
    void flip();

    unsigned int order() const { return _order; }

    void set_order(int order, const Ref<Volume>&vol,double isovalue);
};



class TriangleIntegrator: public DescribedClass {
  private:
    int _n;
    double* _r;
    double* _s;
    double* _w;
    // precomputed interpolation coefficients for triangles of various orders
    Ref<TriInterpCoef> **coef_; // max_order by _n
  protected:
    void set_r(int i,double r);
    void set_s(int i,double s);
    void set_w(int i,double w);
    void init_coef();
    void clear_coef();
  public:
    TriangleIntegrator(const Ref<KeyVal>&);
    TriangleIntegrator(int n);
    virtual ~TriangleIntegrator();
    inline double w(int i) { return _w[i]; }
    inline double r(int i) { return _r[i]; }
    inline double s(int i) { return _s[i]; }
    inline int n() { return _n; }
    virtual void set_n(int n);
    Ref<TriInterpCoef> coef(int order, int i) { return coef_[order-1][i]; }
};


class GaussTriangleIntegrator: public TriangleIntegrator {
  private:
    void init_rw(int order);
  public:
    GaussTriangleIntegrator(const Ref<KeyVal>&);
    GaussTriangleIntegrator(int order);
    ~GaussTriangleIntegrator();
    void set_n(int n);
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
