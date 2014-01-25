//
// surf.cc
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

#include <util/misc/formio.h>
#include <util/keyval/keyval.h>
#include <math/scmat/matrix.h>
#include <math/scmat/vector3.h>
#include <math/isosurf/surf.h>
#include <math/isosurf/isosurf.h>
#include <util/render/polygons.h>

using namespace std;
using namespace sc;

#ifndef WRITE_OOGL
#define WRITE_OOGL 1
#endif

#if WRITE_OOGL
#include <util/render/oogl.h>
#endif

/////////////////////////////////////////////////////////////////////////
// TriangulatedSurface
static ClassDesc TriangulatedSurface_cd(
  typeid(TriangulatedSurface),"TriangulatedSurface",1,"public DescribedClass",
  create<TriangulatedSurface>, create<TriangulatedSurface>, 0);

TriangulatedSurface::TriangulatedSurface():
  _verbose(0),
  _debug(0),
  _triangle_vertex(0),
  _triangle_edge(0),
  _edge_vertex(0),
  _integrator(new GaussTriangleIntegrator(1))
{
  clear();
}

TriangulatedSurface::TriangulatedSurface(const Ref<KeyVal>& keyval):
  _triangle_vertex(0),
  _triangle_edge(0),
  _edge_vertex(0)
{
  _verbose = keyval->booleanvalue("verbose");
  _debug = keyval->booleanvalue("debug");
  Ref<TriangleIntegrator> triint;
  triint << keyval->describedclassvalue("integrator");
  if (triint == 0) {
      triint = new GaussTriangleIntegrator(1);
    }
  set_integrator(triint);
  triint << keyval->describedclassvalue("fast_integrator");
  set_fast_integrator(triint);
  triint << keyval->describedclassvalue("accurate_integrator");
  set_accurate_integrator(triint);
  clear();
}

TriangulatedSurface::~TriangulatedSurface()
{
  clear();
}

void
TriangulatedSurface::topology_info(ostream&o)
{
  topology_info(nvertex(), nedge(), ntriangle(), o);
}

void
TriangulatedSurface::topology_info(int v, int e, int t, ostream&o)
{
  // Given v vertices i expect 2*v - 4*n_surface triangles
  // and 3*v - 6*n_surface edges
  o << indent
    << scprintf("n_vertex = %d, n_edge = %d, n_triangle = %d:",
                v, e, t)
    << endl;
  int nsurf_e = ((3*v - e)%6 == 0)? (3*v - e)/6 : -1;
  int nsurf_t = ((2*v - t)%4 == 0)? (2*v - t)/4 : -1;
  if ((nsurf_e!=-1) && (nsurf_e == nsurf_t)) {
      o << indent
        << scprintf("  this is consistent with n_closed_surface - n_hole = %d",
                    nsurf_e)
        << endl;
    }
  else {
      o << indent
        << scprintf("  this implies that some surfaces are not closed")
        << endl;
    }
}

void
TriangulatedSurface::set_integrator(const Ref<TriangleIntegrator>& i)
{
  _integrator = i;
}

void
TriangulatedSurface::set_fast_integrator(const Ref<TriangleIntegrator>& i)
{
  _fast_integrator = i;
}

void
TriangulatedSurface::set_accurate_integrator(const Ref<TriangleIntegrator>& i)
{
  _accurate_integrator = i;
}

Ref<TriangleIntegrator>
TriangulatedSurface::integrator(int)
{
  // currently the argument, the integer index of the triangle, is ignored
  return _integrator;
}

Ref<TriangleIntegrator>
TriangulatedSurface::fast_integrator(int)
{
  // currently the argument, the integer index of the triangle, is ignored
  return _fast_integrator == 0?_integrator:_fast_integrator;
}

Ref<TriangleIntegrator>
TriangulatedSurface::accurate_integrator(int)
{
  // currently the argument, the integer index of the triangle, is ignored
  return _accurate_integrator == 0?_integrator:_accurate_integrator;
}

void
TriangulatedSurface::clear_int_arrays()
{
  if (_triangle_vertex) {
      for (int i=0; i<_triangles.size(); i++) {
          delete[] _triangle_vertex[i];
        }
      delete[] _triangle_vertex;
    }
  _triangle_vertex = 0;

  if (_triangle_edge) {
      for (int i=0; i<_triangles.size(); i++) {
          delete[] _triangle_edge[i];
        }
      delete[] _triangle_edge;
    }
  _triangle_edge = 0;

  if (_edge_vertex) {
      for (int i=0; i<_edges.size(); i++) {
          delete[] _edge_vertex[i];
        }
      delete[] _edge_vertex;
    }
  _edge_vertex = 0;

  _completed_surface = 0;
}

void
TriangulatedSurface::clear()
{
  _completed_surface = 0;

  clear_int_arrays();

  _have_values = 0;
  _values.clear();

  _vertices.clear();
  _edges.clear();
  _triangles.clear();

  _tmp_edges.clear();
}

void
TriangulatedSurface::complete_surface()
{
  complete_ref_arrays();
  complete_int_arrays();

  _completed_surface = 1;
}

void
TriangulatedSurface::complete_ref_arrays()
{
  _tmp_edges.clear();
  _index_to_edge.clear();
  _edge_to_index.clear();

  int i;
  int ntri = ntriangle();
  _edges.clear();
  for (i=0; i<ntri; i++) {
      Ref<Triangle> tri = triangle(i);
      add_edge(tri->edge(0));
      add_edge(tri->edge(1));
      add_edge(tri->edge(2));
    }
  int ne = nedge();
  _vertices.clear();
  _index_to_vertex.clear();
  _vertex_to_index.clear();
  for (i=0; i<ne; i++) {
      Ref<Edge> e = edge(i);
      add_vertex(e->vertex(0));
      add_vertex(e->vertex(1));
    }
}

void
TriangulatedSurface::complete_int_arrays()
{
  clear_int_arrays();
  
  int i;
  int ntri = ntriangle();
  int ne = nedge();

  // construct the array that converts the triangle number and vertex
  // number within the triangle to the overall vertex number
  _triangle_vertex = new int*[ntri];
  for (i=0; i<ntri; i++) {
      _triangle_vertex[i] = new int[3];
      for (int j=0; j<3; j++) {
          Ref<Vertex> v = triangle(i)->vertex(j);
          _triangle_vertex[i][j] = _vertex_to_index[v];
        }
    }

  // construct the array that converts the triangle number and edge number
  // within the triangle to the overall edge number
  _triangle_edge = new int*[ntri];
  for (i=0; i<ntri; i++) {
      _triangle_edge[i] = new int[3];
      for (int j=0; j<3; j++) {
          Ref<Edge> e = triangle(i)->edge(j);
          _triangle_edge[i][j] = _edge_to_index[e];
        }
    }

  // construct the array that converts the edge number and vertex number
  // within the edge to the overall vertex number
  _edge_vertex = new int*[ne];
  for (i=0; i<ne; i++) {
      _edge_vertex[i] = new int[2];
      for (int j=0; j<2; j++) {
          Ref<Vertex> v = edge(i)->vertex(j);
          _edge_vertex[i][j] = _vertex_to_index[v];
        }
    }
}

void
TriangulatedSurface::compute_values(Ref<Volume>&vol)
{
  int n = _vertices.size();
  _values.resize(n);

  for (int i=0; i<n; i++) {
      vol->set_x(vertex(i)->point());
      _values[i] = vol->value();
    }
  _have_values = 1;
}

double
TriangulatedSurface::flat_area()
{
  double result = 0.0;
  for (std::set<Ref<Triangle> >::iterator i=_triangles.begin();
       i!=_triangles.end(); i++) {
      result += (*i)->flat_area();
    }
  return result;
}

double
TriangulatedSurface::flat_volume()
{
  double result = 0.0;
  for (int i=0; i<_triangles.size(); i++) {

      // get the vertices of the triangle
      SCVector3 A(vertex(triangle_vertex(i,0))->point());
      SCVector3 B(vertex(triangle_vertex(i,1))->point());
      SCVector3 C(vertex(triangle_vertex(i,2))->point());

      // project the vertices onto the xy plane
      SCVector3 Axy(A); Axy[2] = 0.0;
      SCVector3 Bxy(B); Bxy[2] = 0.0;
      SCVector3 Cxy(C); Cxy[2] = 0.0;

      // construct the legs of the triangle in the xy plane
      SCVector3 BAxy = Bxy - Axy;
      SCVector3 CAxy = Cxy - Axy;

      // find the lengths of the legs of the triangle in the xy plane
      double baxy = sqrt(BAxy.dot(BAxy));
      double caxy = sqrt(CAxy.dot(CAxy));

      // if one of the legs is of length zero, then there is
      // no contribution from this triangle
      if (baxy < 1.e-16 || caxy < 1.e-16) continue;

      // find the sine of the angle between the legs of the triangle
      // in the xy plane
      double costheta = BAxy.dot(CAxy)/(baxy*caxy);
      double sintheta = sqrt(1.0 - costheta*costheta);

      // the area of the triangle in the xy plane
      double areaxy = 0.5 * baxy * caxy * sintheta;

      // the height of the three corners of the triangle
      // (relative to the z plane)
      double hA = A[2];
      double hB = B[2];
      double hC = C[2];

      // the volume of the space under the triangle
      double volume = areaxy * (hA + (hB + hC - 2.0*hA)/3.0);

      // the orientation of the triangle along the projection axis (z)
      SCVector3 BA(B-A);
      SCVector3 CA(C-A);
      double z_orientation = BA.cross(CA)[2];

      if (z_orientation > 0.0) {
          result += volume;
        }
      else {
          result -= volume;
        }

    }

  // If the volume is negative, then the surface gradients were
  // opposite in sign to the direction assumed.  Flip the sign
  // to fix.
  return fabs(result);
}

double
TriangulatedSurface::area()
{
  double area = 0.0;
  TriangulatedSurfaceIntegrator triint(this);
  for (triint = 0; triint.update(); triint++) {
      area += triint.w();
    }
  return area;
}

double
TriangulatedSurface::volume()
{
  double volume = 0.0;
  TriangulatedSurfaceIntegrator triint(this);
  for (triint = 0; triint.update(); triint++) {
      volume += triint.weight()*triint.dA()[2]*triint.current()->point()[2];
    }
  return volume;
}

void
TriangulatedSurface::add_vertex(const Ref<Vertex>&t)
{
  int i = _vertices.size();
  _vertices.insert(t);
  if (i != _vertices.size()) {
      _index_to_vertex.push_back(t);
      _vertex_to_index[t] = i;
      if (_index_to_vertex.size() != _vertex_to_index.size()) {
          ExEnv::errn() << "TriangulatedSurface::add_vertex: length mismatch" << endl;
          abort();
        }
    }
}

void
TriangulatedSurface::add_edge(const Ref<Edge>&t)
{
  int i = _edges.size();
  _edges.insert(t);
  if (i != _edges.size()) {
      _index_to_edge.push_back(t);
      _edge_to_index[t] = i;
      if (_index_to_edge.size() != _edge_to_index.size()) {
          ExEnv::errn() << "TriangulatedSurface::add_edge: length mismatch" << endl;
          abort();
        }
    }
}

void
TriangulatedSurface::add_triangle(const Ref<Triangle>&t)
{
  if (_completed_surface) clear();
  int i = _triangles.size();
  _triangles.insert(t);
  if (i != _triangles.size()) {
      _index_to_triangle.push_back(t);
      _triangle_to_index[t] = i;
      if (_index_to_triangle.size() != _triangle_to_index.size()) {
          ExEnv::errn() << "TriangulatedSurface::add_triangle: length mismatch" << endl;
          abort();
        }
    }
}

void
TriangulatedSurface::add_triangle(const Ref<Vertex>& v1,
                                  const Ref<Vertex>& v2,
                                  const Ref<Vertex>& v3)
{
  // Find this triangle's edges if they have already be created
  // for some other triangle.
  Ref<Edge> e0, e1, e2;

  const std::set<Ref<Edge> > &v1edges = _tmp_edges[v1];

  const std::set<Ref<Edge> > &v2edges = _tmp_edges[v2];

  std::set<Ref<Edge> >::const_iterator ix;
  for (ix = v1edges.begin(); ix != v1edges.end(); ix++) {
      const Ref<Edge>& e = *ix;
      if (e->vertex(0) == v2 || e->vertex(1) == v2) {
          e0 = e;
        }
      else if (e->vertex(0) == v3 || e->vertex(1) == v3) {
          e2 = e;
        }
    }
  for (ix = v2edges.begin(); ix != v2edges.end(); ix++) {
      const Ref<Edge>& e = *ix;
      if (e->vertex(0) == v3 || e->vertex(1) == v3) {
          e1 = e;
        }
    }

  if (e0 == 0) {
      e0 = newEdge(v1,v2);
      _tmp_edges[v1].insert(e0);
      _tmp_edges[v2].insert(e0);
    }
  if (e1 == 0) {
      e1 = newEdge(v2,v3);
      _tmp_edges[v2].insert(e1);
      _tmp_edges[v3].insert(e1);
    }
  if (e2 == 0) {
      e2 = newEdge(v3,v1);
      _tmp_edges[v3].insert(e2);
      _tmp_edges[v1].insert(e2);
    }
  
  int orientation;
  if (e0->vertex(0) == v1) {
      orientation = 0;
    }
  else {
      orientation = 1;
    }
  
  add_triangle(newTriangle(e0,e1,e2,orientation));
}

// If a user isn't keeping track of edges while add_triangle is being
// used to build the surface, then this can be called to see if an edge
// already exists (at a great performance cost).
Ref<Edge>
TriangulatedSurface::find_edge(const Ref<Vertex>& v1, const Ref<Vertex>& v2)
{
  std::set<Ref<Triangle> >::iterator i;

  for (i=_triangles.begin(); i!=_triangles.end(); i++) {
      Ref<Triangle> t = *i;
      Ref<Edge> e1 = t->edge(0);
      Ref<Edge> e2 = t->edge(1);
      Ref<Edge> e3 = t->edge(2);
      if (e1->vertex(0) == v1 && e1->vertex(1) == v2) return e1;
      if (e1->vertex(1) == v1 && e1->vertex(0) == v2) return e1;
      if (e2->vertex(0) == v1 && e2->vertex(1) == v2) return e2;
      if (e2->vertex(1) == v1 && e2->vertex(0) == v2) return e2;
      if (e3->vertex(0) == v1 && e3->vertex(1) == v2) return e3;
      if (e3->vertex(1) == v1 && e3->vertex(0) == v2) return e3;
    }

  return 0;
}

void
TriangulatedSurface::print(ostream&o) const
{
  o << indent << "TriangulatedSurface:" << endl;
  int i;

  int np = nvertex();
  o << indent << scprintf(" %3d Vertices:",np) << endl;
  for (i=0; i<np; i++) {
      Ref<Vertex> p = vertex(i);
      o << indent << scprintf("  %3d:",i);
      for (int j=0; j<3; j++) {
          o << scprintf(" % 15.10f", p->point()[j]);
        }
      o << endl;
    }

  int ne = nedge();
  o << indent << scprintf(" %3d Edges:",ne) << endl;
  for (i=0; i<ne; i++) {
      Ref<Edge> e = edge(i);
      Ref<Vertex> v0 = e->vertex(0);
      Ref<Vertex> v1 = e->vertex(1);
      std::map<Ref<Vertex>,int>::const_iterator v0i=_vertex_to_index.find(v0);
      std::map<Ref<Vertex>,int>::const_iterator v1i=_vertex_to_index.find(v1);
      int v0int = v0i==_vertex_to_index.end()? -1: v0i->second;
      int v1int = v1i==_vertex_to_index.end()? -1: v1i->second;
      o << indent
        << scprintf("  %3d: %3d %3d",i, v0int, v1int)
        << endl;
    }

  int nt = ntriangle();
  o << indent << scprintf(" %3d Triangles:",nt) << endl;
  for (i=0; i<nt; i++) {
      Ref<Triangle> tri = triangle(i);
      Ref<Edge> e0 = tri->edge(0);
      Ref<Edge> e1 = tri->edge(1);
      Ref<Edge> e2 = tri->edge(2);
      std::map<Ref<Edge>,int>::const_iterator e0i = _edge_to_index.find(e0);
      std::map<Ref<Edge>,int>::const_iterator e1i = _edge_to_index.find(e1);
      std::map<Ref<Edge>,int>::const_iterator e2i = _edge_to_index.find(e2);
      int e0int = e0i==_edge_to_index.end()? -1: e0i->second;
      int e1int = e1i==_edge_to_index.end()? -1: e1i->second;
      int e2int = e2i==_edge_to_index.end()? -1: e2i->second;
      o << indent
        << scprintf("  %3d: %3d %3d %3d",i, e0int, e1int, e2int)
        << endl;
    }
}

void
TriangulatedSurface::print_vertices_and_triangles(ostream&o) const
{
  o << indent << "TriangulatedSurface:" << endl;
  int i;

  int np = nvertex();
  o << indent << scprintf(" %3d Vertices:",np) << endl;
  for (i=0; i<np; i++) {
      Ref<Vertex> p = vertex(i);
      o << indent << scprintf("  %3d:",i);
      for (int j=0; j<3; j++) {
          o << scprintf(" % 15.10f", p->point()[j]);
        }
      o << endl;
    }

  int nt = ntriangle();
  o << indent << scprintf(" %3d Triangles:",nt) << endl;
  for (i=0; i<nt; i++) {
      Ref<Triangle> tri = triangle(i);
      o << indent
        << scprintf("  %3d: %3d %3d %3d",i,
                    _triangle_vertex[i][0],
                    _triangle_vertex[i][1],
                    _triangle_vertex[i][2])
        << endl;
    }
}

void
TriangulatedSurface::render(const Ref<Render> &render)
{
  Ref<RenderedPolygons> poly = new RenderedPolygons;
  poly->initialize(_vertices.size(), _triangles.size(),
                   RenderedPolygons::Vertex);
  std::set<Ref<Vertex> >::iterator iv;
  std::set<Ref<Triangle> >::iterator it;
  std::map<Ref<Vertex>, int> vertex_to_index;
  int i = 0;
  for (iv = _vertices.begin(); iv != _vertices.end(); iv++, i++) {
      Ref<Vertex> v = *iv;
      vertex_to_index[v] = i;
      poly->set_vertex(i,
                       v->point()[0],
                       v->point()[1],
                       v->point()[2]);
      poly->set_vertex_rgb(i, 0.3, 0.3, 0.3);
    }
  i = 0;
  for (it = _triangles.begin(); it != _triangles.end(); it++, i++) {
      Ref<Triangle> t = *it;
      poly->set_face(i,
                     vertex_to_index[t->vertex(0)],
                     vertex_to_index[t->vertex(1)],
                     vertex_to_index[t->vertex(2)]);
    }
  render->render(poly.pointer());
}

void
TriangulatedSurface::print_geomview_format(ostream&o) const
{
  o << "OFF" << endl;

  o << nvertex() << " " << ntriangle() << " " << nedge() << endl;
  int i;

  int np = nvertex();
  for (i=0; i<np; i++) {
      Ref<Vertex> p = vertex(i);
      for (int j=0; j<3; j++) {
          o << scprintf(" % 15.10f", p->point()[j]);
        }
      o << endl;
    }

  int nt = ntriangle();
  for (i=0; i<nt; i++) {
      Ref<Triangle> tri = triangle(i);
      o << scprintf(" 3 %3d %3d %3d",
                    _triangle_vertex[i][0],
                    _triangle_vertex[i][1],
                    _triangle_vertex[i][2])
          << endl;
    }
}

void
TriangulatedSurface::recompute_index_maps()
{
  int i;
  std::set<Ref<Vertex> >::iterator iv;
  std::set<Ref<Edge> >::iterator ie;
  std::set<Ref<Triangle> >::iterator it;

  // fix the index maps
  _vertex_to_index.clear();
  _edge_to_index.clear();
  _triangle_to_index.clear();

  _index_to_vertex.clear();
  _index_to_edge.clear();
  _index_to_triangle.clear();

  _index_to_vertex.resize(_vertices.size());
  for (i=0, iv = _vertices.begin(); iv != _vertices.end(); i++, iv++) {
      _vertex_to_index[*iv] = i;
      _index_to_vertex[i] = *iv;
    }

  _index_to_edge.resize(_edges.size());
  for (i=0, ie = _edges.begin(); ie != _edges.end(); i++, ie++) {
      _edge_to_index[*ie] = i;
      _index_to_edge[i] = *ie;
    }

  _index_to_triangle.resize(_triangles.size());
  for (i=0, it = _triangles.begin(); it != _triangles.end(); i++, it++) {
      _triangle_to_index[*it] = i;
      _index_to_triangle[i] = *it;
    }

}

Edge*
TriangulatedSurface::newEdge(const Ref<Vertex>& v0, const Ref<Vertex>& v1) const
{
  return new Edge(v0,v1);
}

Triangle*
TriangulatedSurface::newTriangle(const Ref<Edge>& e0,
                                 const Ref<Edge>& e1,
                                 const Ref<Edge>& e2,
                                 int orientation) const
{
  return new Triangle(e0,e1,e2,orientation);
}

//////////////////////////////////////////////////////////////////////
// TriangulatedSurfaceIntegrator

TriangulatedSurfaceIntegrator::
  TriangulatedSurfaceIntegrator(const Ref<TriangulatedSurface>&ts)
{
  set_surface(ts);
  use_default_integrator();

  _itri = 0;
  _irs = 0;
}

TriangulatedSurfaceIntegrator::
  TriangulatedSurfaceIntegrator()
{
  use_default_integrator();
}

void
TriangulatedSurfaceIntegrator::
  operator =(const TriangulatedSurfaceIntegrator&i)
{
  set_surface(i._ts);

  _integrator = i._integrator;
  _itri = i._itri;
  _irs = i._irs;
}

TriangulatedSurfaceIntegrator::
  ~TriangulatedSurfaceIntegrator()
{
}

void
TriangulatedSurfaceIntegrator::set_surface(const Ref<TriangulatedSurface>&s)
{
  _ts = s;
  _current = new Vertex();
}

int
TriangulatedSurfaceIntegrator::
  vertex_number(int i)
{
  return _ts->triangle_vertex(_itri,i);
}

Ref<Vertex>
TriangulatedSurfaceIntegrator::
  current()
{
  return _current;
}

int
TriangulatedSurfaceIntegrator::n()
{
  int result = 0;
  int ntri = _ts->ntriangle();
  for (int i=0; i<ntri; i++) {
      result += (_ts.pointer()->*_integrator)(i)->n();
    }
  return result;
}

int
TriangulatedSurfaceIntegrator::update()
{
  if (_itri < 0 || _itri >= _ts->ntriangle()) return 0;

  TriangleIntegrator* i = (_ts.pointer()->*_integrator)(_itri).pointer();
  _s = i->s(_irs);
  _r = i->r(_irs);
  _weight = i->w(_irs);
  Ref<Triangle> t = _ts->triangle(_itri);
  Ref<TriInterpCoef> coef = i->coef(t->order(),_irs);
  t->interpolate(coef, _r, _s, _current, _dA);
  _surface_element = _dA.norm();
  static double cum;
  if (_irs == 0) cum = 0.0;
  cum += _surface_element * _weight;
  //ExEnv::outn() << scprintf("%2d dA = %12.8f, w = %12.8f, Sum wdA = %12.8f",
  //                 _irs, _surface_element, _weight, cum)
  //     << endl;

  return (int) 1;
}

void
TriangulatedSurfaceIntegrator::
  operator ++()
{
  int n = (_ts.pointer()->*_integrator)(_itri)->n();
  if (_irs == n-1) {
      _irs = 0;
      if (_grp == 0) _itri++;
      else _itri += _grp->n();
    }
  else {
      _irs++;
    }
}

void
TriangulatedSurfaceIntegrator::distribute(const Ref<MessageGrp> &grp)
{
  _grp = grp;
}

int
TriangulatedSurfaceIntegrator::
  operator = (int i)
{
  _itri = i;
  _irs = 0;
  return i;
}

void
TriangulatedSurfaceIntegrator::use_default_integrator()
{
  _integrator = &TriangulatedSurface::integrator;
}

void
TriangulatedSurfaceIntegrator::use_fast_integrator()
{
  _integrator = &TriangulatedSurface::fast_integrator;
}

void
TriangulatedSurfaceIntegrator::use_accurate_integrator()
{
  _integrator = &TriangulatedSurface::accurate_integrator;
}

/////////////////////////////////////////////////////////////////////////
// TriangulatedImplicitSurface
static ClassDesc TriangulatedImplicitSurface_cd(
  typeid(TriangulatedImplicitSurface),"TriangulatedImplicitSurface",1,"public TriangulatedSurface",
  0, create<TriangulatedImplicitSurface>, 0);

TriangulatedImplicitSurface::
TriangulatedImplicitSurface(const Ref<KeyVal>&keyval):
  TriangulatedSurface(keyval)
{
  inited_ = 0;

  vol_ << keyval->describedclassvalue("volume");
  if (keyval->error() != KeyVal::OK) {
      ExEnv::errn() << "TriangulatedImplicitSurface(const Ref<KeyVal>&keyval): "
           << "requires \"volume\"" << endl;
      abort();
    }

  isovalue_ = keyval->doublevalue("value");
  if (keyval->error() != KeyVal::OK) isovalue_ = 0.0;

  fix_orientation_ = keyval->booleanvalue("fix_orientation");
  if (keyval->error() != KeyVal::OK) fix_orientation_ = 1;

  remove_short_edges_ = keyval->booleanvalue("remove_short_edges");
  if (keyval->error() != KeyVal::OK) remove_short_edges_ = 1;

  remove_slender_triangles_ = keyval->booleanvalue("remove_slender_triangles");
  if (keyval->error() != KeyVal::OK) remove_slender_triangles_ = 0;

  remove_small_triangles_ = keyval->booleanvalue("remove_small_triangles");
  if (keyval->error() != KeyVal::OK) remove_small_triangles_ = 0;

  short_edge_factor_ = keyval->doublevalue("short_edge_factor");
  if (keyval->error() != KeyVal::OK) short_edge_factor_ = 0.4;

  slender_triangle_factor_ = keyval->doublevalue("slender_triangle_factor");
  if (keyval->error() != KeyVal::OK) slender_triangle_factor_ = 0.2;

  small_triangle_factor_ = keyval->doublevalue("small_triangle_factor");
  if (keyval->error() != KeyVal::OK) small_triangle_factor_ = 0.2;

  resolution_ = keyval->doublevalue("resolution");
  if (keyval->error() != KeyVal::OK) resolution_ = 1.0;

  order_ = keyval->intvalue("order");
  if (keyval->error() != KeyVal::OK) order_ = 1;

  int initialize = keyval->booleanvalue("initialize");
  if (initialize) init();
}

void
TriangulatedImplicitSurface::init()
{
  if (_verbose) ExEnv::outn() << "TriangulatedImplicitSurface: init start" << endl;
  ImplicitSurfacePolygonizer isogen(vol_);
  isogen.set_resolution(resolution_);

  if (_verbose) ExEnv::outn() << "TriangulatedImplicitSurface: isosurface" << endl;
  isogen.isosurface(isovalue_,*this);
#if WRITE_OOGL
  if (_debug) {
      render(new OOGLRender("surfiso.oogl"));
    }
#endif
  if (_verbose) ExEnv::outn() << "TriangulatedImplicitSurface: orientation" << endl;
  if (fix_orientation_) fix_orientation();
#if WRITE_OOGL
  if (_debug) {
      render(new OOGLRender("surffix.oogl"));
    }
#endif
  if (remove_short_edges_) {
      if (_verbose) ExEnv::outn() << "TriangulatedImplicitSurface: short edges" << endl;
      remove_short_edges(short_edge_factor_*resolution_,vol_,isovalue_);
      if (_verbose) ExEnv::outn() << "TriangulatedImplicitSurface: orientation" << endl;
      if (fix_orientation_) fix_orientation();
    }
  if (remove_slender_triangles_ || remove_small_triangles_) {
      if (_verbose) ExEnv::outn() << "TriangulatedImplicitSurface: slender" << endl;
      double height_cutoff = slender_triangle_factor_ * resolution_;
      double area_cutoff = small_triangle_factor_*resolution_*resolution_*0.5;
      remove_slender_triangles(remove_slender_triangles_, height_cutoff,
                               remove_small_triangles_, area_cutoff,
                               vol_,isovalue_);
      if (_verbose) ExEnv::outn() << "TriangulatedImplicitSurface: orientation" << endl;
      if (fix_orientation_) fix_orientation();
    }

  // see if a higher order approximation to the surface is required
  if (order_ > 1) {
      if (_verbose) ExEnv::outn() << "TriangulatedImplicitSurface: order" << endl;
      int i;
      for (i=0; i<nedge(); i++) {
          edge(i)->set_order(order_, vol_, isovalue_);
        }
      for (i=0; i<ntriangle(); i++) {
          triangle(i)->set_order(order_, vol_, isovalue_);
        }
    }
  inited_ = 1;
  if (_verbose) ExEnv::outn() << "TriangulatedImplicitSurface: init done" << endl;
}

TriangulatedImplicitSurface::~TriangulatedImplicitSurface()
{
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
