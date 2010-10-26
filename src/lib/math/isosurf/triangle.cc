//
// triangle.cc
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

#ifdef __GNUC__
#pragma implementation
#endif

#include <util/misc/formio.h>
#include <util/keyval/keyval.h>
#include <math/isosurf/triangle.h>
#include <math/scmat/vector3.h>

using namespace std;
using namespace sc;

/////////////////////////////////////////////////////////////////////////
// Triangle
///////////////////////////////////////////////////////////////////////////
// Here is the layout of the vertices and edges in the triangles:         |
//                                                                        |
//                       vertex(1) (r=0, s=1)                             |
//                          +                                             |
//                         / \  _edges[1].vertex(_orientation1)           |
//                        /   \                                           |
//                       /     \                                          |
//                      /       \                                         |
//                     /         \                                        |
//                    /           \                                       |
//         _edges[0] /             \  _edges[1]                           |
//          (r = 0) /               \   (1-r-s = 0)                       |
//                 /                 \                                    |
//                /                   \                                   |
//               /                     \                                  |
//              /                       \ _edges[1].vertex(!_orientation1)|
//             /                         \                                |
//   vertex(0)+---------------------------+                               |
// (r=0, s=0)        _edges[2] (s = 0)      vertex(2) (r=1, s=0)          |
//                                                                        |
//  Zienkiewicz and Taylor, "The Finite Element Method", 4th Ed, Vol 1,   |
//  use                                                                   |
//      L1 = 1 - r - s                                                    |
//      L2 = r,                                                           |
//      L3 = s.                                                           |
//  I also use these below.                                               |
///////////////////////////////////////////////////////////////////////////

Triangle::Triangle(const Ref<Edge>& v1, const Ref<Edge>& v2, const Ref<Edge>& v3,
                   unsigned int orientation0)
{
  _orientation0 = orientation0;

  _edges[0] = v1;
  _edges[1] = v2;
  _edges[2] = v3;

  // edge 0 corresponds to r = 0
  // edge 1 corresponds to (1-r-s) = 0
  // edge 2 corresponds to s = 0
  // edge 0 vertex _orientation0 is (r=0,s=0)
  // edge 1 vertex _orientation1 is (r=0,s=1)
  // edge 2 vertex _orientation2 is (r=1,s=0)
  // edge 0 vertex (1-_orientation0) is edge 1 vertex _orientation1
  // edge 1 vertex (1-_orientation1) is edge 2 vertex _orientation2
  // edge 2 vertex (1-_orientation2) is edge 0 vertex _orientation0

  Ref<Edge> *e = _edges;

  // swap edges 1 and 2 if necessary
  if (e[0]->vertex(1-_orientation0) != e[1]->vertex(0)) {
      if (e[0]->vertex(1-_orientation0) != e[1]->vertex(1)) {
          e[1] = v3;
          e[2] = v2;
        }
    }

  // compute the orientation of _edge[1]
  if (e[0]->vertex(1-_orientation0) == e[1]->vertex(0)) {
      _orientation1 = 0;
    }
  else {
      _orientation1 = 1;
    }

  // compute the orientation of _edge[2]
  if (e[1]->vertex(1-_orientation1) == e[2]->vertex(0)) {
      _orientation2 = 0;
    }
  else {
      _orientation2 = 1;
    }

  // make sure that the triangle is really a triangle
  if ( e[0]->vertex(1-_orientation0) != e[1]->vertex(_orientation1)
       || e[1]->vertex(1-_orientation1) != e[2]->vertex(_orientation2)
       || e[2]->vertex(1-_orientation2) != e[0]->vertex(_orientation0))
    {
      ExEnv::errn() << "Triangle: given edges that don't form a triangle" << endl;
      abort();
    }

  _order = 1;
  _vertices = new Ref<Vertex>[3];

  _vertices[TriInterpCoef::ijk_to_index(_order, 0, 0)] = vertex(0);
  _vertices[TriInterpCoef::ijk_to_index(0, 0, _order)] = vertex(1);
  _vertices[TriInterpCoef::ijk_to_index(0, _order, 0)] = vertex(2);
}

Triangle::~Triangle()
{
  if (_vertices) delete[] _vertices;
}

Ref<Vertex>
Triangle::vertex(int i)
{
  return _edges[i]->vertex(orientation(i));
}

unsigned int
Triangle::orientation(const Ref<Edge>& e) const
{
  if (e == _edges[0]) return orientation(0);
  if (e == _edges[1]) return orientation(1);
  return orientation(2);
}

int
Triangle::contains(const Ref<Edge>& e) const
{
  if (_edges[0] == e) return 1;
  if (_edges[1] == e) return 1;
  if (_edges[2] == e) return 1;
  return 0;
}

double Triangle::flat_area()
{
  double a = _edges[0]->straight_length();
  double b = _edges[1]->straight_length();
  double c = _edges[2]->straight_length();
  double a2 = a*a;
  double b2 = b*b;
  double c2 = c*c;
  return 0.25 * sqrt( 2.0 * (c2*b2 + c2*a2 + a2*b2)
                      - c2*c2 - b2*b2 - a2*a2);
}

void Triangle::add_vertices(std::set<Ref<Vertex> >&set)
{
  for (int i=0; i<3; i++) set.insert(_edges[i]->vertex(orientation(i)));
}

void Triangle::add_edges(std::set<Ref<Edge> >&set)
{
  for (int i=0; i<3; i++) set.insert(_edges[i]);
}

void
Triangle::interpolate(double r,double s,const Ref<Vertex>&result,SCVector3&dA)
{
  TriInterpCoefKey key(_order, r, s);
  Ref<TriInterpCoef> coef = new TriInterpCoef(key);
  interpolate(coef, r, s, result, dA);
}

void
Triangle::interpolate(const Ref<TriInterpCoef>& coef,
                      double r, double s,
                      const Ref<Vertex>&result, SCVector3&dA)
{
  unsigned int i, j, k;

  //double L1 = 1 - r - s;
  //double L2 = r;
  //double L3 = s;

  SCVector3 tmp(0.0);
  SCVector3 x_s(0.0);
  SCVector3 x_r(0.0);
  for (i=0; i<=_order; i++) {
      for (j=0; j <= _order-i; j++) {
          k = _order - i - j;
          int index = TriInterpCoef::ijk_to_index(i,j,k);
          tmp += coef->coef(i,j,k)*_vertices[index]->point();
          x_s += coef->sderiv(i,j,k)*_vertices[index]->point();
          x_r += coef->rderiv(i,j,k)*_vertices[index]->point();
        }
    }
  result->set_point(tmp);

  if (result->has_normal()) {
      for (i=0; i<_order; i++) {
          for (j=0; j <= _order-i; j++) {
              k = _order - i - j;
              int index = TriInterpCoef::ijk_to_index(i,j,k);
              tmp += coef->coef(i,j,k)*_vertices[index]->point();
            }
        }
      result->set_normal(tmp);
    }

  // Find the surface element
  dA = x_r.cross(x_s);
}

void
Triangle::interpolate(double r, double s,
                      const Ref<Vertex>&result, SCVector3&dA,
                      const Ref<Volume> &vol, double isoval)
{
  // set up an initial dummy norm
  SCVector3 norm(0.0);
  result->set_normal(norm);

  // initial guess
  interpolate(r,s,result,dA);

  // now refine that guess
  SCVector3 trialpoint = result->point();
  SCVector3 trialnorm = result->normal();
  SCVector3 newpoint;
  vol->solve(trialpoint,trialnorm,isoval,newpoint);
  // compute the true normal
  vol->set_x(newpoint);
  if (vol->gradient_implemented()) {
      vol->get_gradient(trialnorm);
    }
  trialnorm.normalize();
  result->set_point(newpoint);
  result->set_normal(trialnorm);
}

void
Triangle::flip()
{
  _orientation0 = _orientation0?0:1;
  _orientation1 = _orientation1?0:1;
  _orientation2 = _orientation2?0:1;
}

void
Triangle::set_order(int order, const Ref<Volume>&vol, double isovalue)
{
  if (order > max_order) {
      ExEnv::errn() << scprintf("Triangle::set_order: max_order = %d but order = %d\n",
                       max_order, order);
      abort();
    }

  unsigned int i, j, k;

  if (edge(0)->order() != order
      ||edge(1)->order() != order
      ||edge(2)->order() != order) {
      ExEnv::errn() << "Triangle::set_order: edge order doesn't match" << endl;
      abort();
    }

  _order = order;
  delete[] _vertices;

  _vertices = new Ref<Vertex>[TriInterpCoef::order_to_nvertex(_order)];

  // fill in the corner vertices
  _vertices[TriInterpCoef::ijk_to_index(_order, 0, 0)] = vertex(0);
  _vertices[TriInterpCoef::ijk_to_index(0, 0, _order)] = vertex(1);
  _vertices[TriInterpCoef::ijk_to_index(0, _order, 0)] = vertex(2);

  // fill in the interior edge vertices
  for (i = 1; i < _order; i++) {
      j = _order - i;
      _vertices[TriInterpCoef::ijk_to_index(0, i, j)]
          = _edges[1]->interior_vertex(_orientation1?i:j);
      _vertices[TriInterpCoef::ijk_to_index(j, 0, i)]
          = _edges[0]->interior_vertex(_orientation0?i:j);
      _vertices[TriInterpCoef::ijk_to_index(i, j, 0)]
          = _edges[2]->interior_vertex(_orientation2?i:j);
    }

  const SCVector3& p0 = vertex(0)->point();
  const SCVector3& p1 = vertex(1)->point();
  const SCVector3& p2 = vertex(2)->point();
  const SCVector3& norm0 = vertex(0)->normal();
  const SCVector3& norm1 = vertex(1)->normal();
  const SCVector3& norm2 = vertex(2)->normal();

  for (i=0; i<=_order; i++) {
      double I = (1.0*i)/_order;
      for (j=0; j<=_order-i; j++) {
          SCVector3 trialpoint;
          SCVector3 trialnorm;
          SCVector3 newpoint;
          double J = (1.0*j)/_order;
          k = _order - i - j;
          if (!i || !j || !k) continue; // interior point check
          double K = (1.0*k)/_order;
          int index = TriInterpCoef::ijk_to_index(i,j,k);
          // this get approximate vertices and normals
          trialpoint = I*p0 + J*p1 + K*p2;
          trialnorm = I*norm0 + J*norm1 + K*norm2;
          // now refine that guess
          vol->solve(trialpoint,trialnorm,isovalue,newpoint);
          // compute the true normal
          vol->set_x(newpoint);
          if (vol->gradient_implemented()) {
              vol->get_gradient(trialnorm);
            }
          trialnorm.normalize();
          _vertices[index] = new Vertex(newpoint,trialnorm);
        }
    }
}

/////////////////////////////////////////////////////////////////////////
// TriangleIntegrator

static ClassDesc TriangleIntegrator_cd(
  typeid(TriangleIntegrator),"TriangleIntegrator",1,"public DescribedClass",
  0, create<TriangleIntegrator>, 0);

TriangleIntegrator::TriangleIntegrator(int order):
  _n(order)
{
  _r = new double [_n];
  _s = new double [_n];
  _w = new double [_n];
  coef_ = 0;
}

TriangleIntegrator::TriangleIntegrator(const Ref<KeyVal>& keyval)
{
  _n = keyval->intvalue("n");
  if (keyval->error() != KeyVal::OK) {
      _n = 7;
    }
  _r = new double [_n];
  _s = new double [_n];
  _w = new double [_n];
  coef_ = 0;
}

TriangleIntegrator::~TriangleIntegrator()
{
  delete[] _r;
  delete[] _s;
  delete[] _w;
  clear_coef();
}

void
TriangleIntegrator::set_n(int n)
{
  delete[] _r;
  delete[] _s;
  delete[] _w;
  _n = n;
  _r = new double [_n];
  _s = new double [_n];
  _w = new double [_n];
  clear_coef();
}

void
TriangleIntegrator::set_w(int i,double w)
{
  _w[i] = w;
}

void
TriangleIntegrator::set_r(int i,double r)
{
  _r[i] = r;
}

void
TriangleIntegrator::set_s(int i,double s)
{
  _s[i] = s;
}

void
TriangleIntegrator::init_coef()
{
  int i, j;

  clear_coef();

  coef_ = new Ref<TriInterpCoef>*[Triangle::max_order];
  for (i=0; i<Triangle::max_order; i++) {
      coef_[i] = new Ref<TriInterpCoef>[_n];
      for (j=0; j<_n; j++) {
          TriInterpCoefKey key(i+1,_r[j],_s[j]);
          coef_[i][j] = new TriInterpCoef(key);
        }
    }
}

void
TriangleIntegrator::clear_coef()
{
  int i, j;

  if (coef_) {
      for (i=0; i<Triangle::max_order; i++) {
          for (j=0; j<_n; j++) {
              coef_[i][j] = 0;
            }
          delete[] coef_[i];
        }
    }
  delete[] coef_;
  coef_ = 0;
}

/////////////////////////////////////////////////////////////////////////
// GaussTriangleIntegrator

static ClassDesc GaussTriangleIntegrator_cd(
  typeid(GaussTriangleIntegrator),"GaussTriangleIntegrator",1,"public TriangleIntegrator",
  0, create<GaussTriangleIntegrator>, 0);

GaussTriangleIntegrator::GaussTriangleIntegrator(const Ref<KeyVal>& keyval):
  TriangleIntegrator(keyval)
{
  ExEnv::out0() << "Created a GaussTriangleIntegrator with n = " << n() << endl;
  init_rw(n());
  init_coef();
}

GaussTriangleIntegrator::GaussTriangleIntegrator(int order):
  TriangleIntegrator(order)
{
  init_rw(n());
  init_coef();
}

void
GaussTriangleIntegrator::set_n(int n)
{
  TriangleIntegrator::set_n(n);
  init_rw(n);
  init_coef();
}

void
GaussTriangleIntegrator::init_rw(int order)
{
  if (order == 1) {
      set_r(0, 1.0/3.0);
      set_s(0, 1.0/3.0);
      set_w(0, 1.0);
    }
  else if (order == 3) {
      set_r(0, 1.0/6.0);
      set_r(1, 2.0/3.0);
      set_r(2, 1.0/6.0);
      set_s(0, 1.0/6.0);
      set_s(1, 1.0/6.0);
      set_s(2, 2.0/3.0);
      set_w(0, 1.0/3.0);
      set_w(1, 1.0/3.0);
      set_w(2, 1.0/3.0);
    }
  else if (order == 4) {
      set_r(0, 1.0/3.0);
      set_r(1, 1.0/5.0);
      set_r(2, 3.0/5.0);
      set_r(3, 1.0/5.0);
      set_s(0, 1.0/3.0);
      set_s(1, 1.0/5.0);
      set_s(2, 1.0/5.0);
      set_s(3, 3.0/5.0);
      set_w(0, -27.0/48.0);
      set_w(1, 25.0/48.0);
      set_w(2, 25.0/48.0);
      set_w(3, 25.0/48.0);
    }
  else if (order == 6) {
      set_r(0, 0.091576213509771);
      set_r(1, 0.816847572980459);
      set_r(2, r(0));
      set_r(3, 0.445948490915965);
      set_r(4, 0.108103018168070);
      set_r(5, r(3));
      set_s(0, r(0));
      set_s(1, r(0));
      set_s(2, r(1));
      set_s(3, r(3));
      set_s(4, r(3));
      set_s(5, r(4));
      set_w(0, 0.109951743655322);
      set_w(1, w(0));
      set_w(2, w(0));
      set_w(3, 0.223381589678011);
      set_w(4, w(3));
      set_w(5, w(3));
    }
  else if (order == 7) {
      set_r(0, 1.0/3.0);
      set_r(1, 0.101286507323456);
      set_r(2, 0.797426985353087);
      set_r(3, r(1));
      set_r(4, 0.470142064105115);
      set_r(5, 0.059715871789770);
      set_r(6, r(4));
      set_s(0, r(0));
      set_s(1, r(1));
      set_s(2, r(1));
      set_s(3, r(2));
      set_s(4, r(4));
      set_s(5, r(4));
      set_s(6, r(5));
      set_w(0, 0.225);
      set_w(1, 0.125939180544827);
      set_w(2, w(1));
      set_w(3, w(1));
      set_w(4, 0.132394152788506);
      set_w(5, w(4));
      set_w(6, w(4));
    }
  else {
      ExEnv::errn() << "GaussTriangleIntegrator: invalid order " << order << endl;
      abort();
    }

  // scale the weights by the area of the unit triangle
  for (int i=0; i<order; i++) {
      set_w(i, w(i)*0.5);
    }
}

GaussTriangleIntegrator::~GaussTriangleIntegrator()
{
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
