//
// isosurf.cc
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

#include <iostream>

#include <math/isosurf/isosurf.h>
#include <math/isosurf/implicit.h>
#include <math/scmat/vector3.h>

using namespace std;
using namespace sc;

// This static datum is used to interface to the implicit.c routine
// provided in Graphics Gems IV.
ImplicitSurfacePolygonizer* ImplicitSurfacePolygonizer::current = 0;

// These functions are used to interface to the implicit.c routine provided
// in Graphics Gems IV.
extern "C" int
ImplicitSurfacePolygonizer_add_triangle_to_current(int i1, int i2, int i3,
                                                   sc::detail::VERTICES v)
{
  return ImplicitSurfacePolygonizer::add_triangle_to_current(i1, i2, i3, v);
}

extern "C" double
ImplicitSurfacePolygonizer_value_of_current(double x,double y,double z)
{
  return ImplicitSurfacePolygonizer::value_of_current(x,y,z);
}

////////////////////////////////////////////////////////////////////////////
// IsosurfaceGen members

IsosurfaceGen::IsosurfaceGen():
  _resolution(0.25)
{
}

IsosurfaceGen::~IsosurfaceGen()
{
}

void
IsosurfaceGen::set_resolution(double r)
{
  _resolution = r;
}

////////////////////////////////////////////////////////////////////////////
// ImplicitSurfacePolygonizer members

ImplicitSurfacePolygonizer::ImplicitSurfacePolygonizer(const Ref<Volume>&vol):
  _volume(vol)
  ,_tmp_vertices(0)
{
}

ImplicitSurfacePolygonizer::~ImplicitSurfacePolygonizer()
{
}

static SCVector3 current_x;
void
ImplicitSurfacePolygonizer::isosurface(double value,
                                       TriangulatedSurface& surf)
{
  surf.clear();

  // Find the bounding box.
  SCVector3 p0;
  SCVector3 p1;
  _volume->boundingbox(value - 0.001, value + 0.001, p0, p1);
  SCVector3 diag = p1 - p0;
  SCVector3 midpoint = 0.5*diag + p0;
  double biggest_width = diag.maxabs();
  int bounds = (int)(0.5*biggest_width/_resolution) + 2;

  // polygonize will find a starting point and do bounds checking
  // from that point.  To make sure the bounding box is big enough
  // its size must be doubled.  Since polygonization is implicit
  // there is no performance penalty.
  bounds *= 2;

  // Initialize the static pointer to this, so the C polygonizer can find us.
  current = this;
  _surf = &surf;
  _value = value;
  // Find the polygons.
  char *msg = polygonize(ImplicitSurfacePolygonizer_value_of_current, _resolution, bounds,
                         midpoint[0], midpoint[1], midpoint[2],
                         ImplicitSurfacePolygonizer_add_triangle_to_current,
                         NOTET);
  current = 0;
  _surf = 0;
  if (msg) {
      ExEnv::errn() << "ImplicitSurfacePolygonizer::isosurface: failed: "
           << msg << endl;
      abort();
    }

  // Clean up temporaries.
  _tmp_vertices.clear();

  ExEnv::out0() << "about to complete the surface" << endl;

  // finish the surface
  surf.complete_surface();

  ExEnv::out0() << "completed the surface" << endl;
  ExEnv::out0() << "flat area = " << surf.flat_area() << endl;
  ExEnv::out0() << "  ntri = " << setw(10) << surf.ntriangle()
       << " bytes = "
       << setw(10) << surf.ntriangle() * sizeof(Triangle)
       << endl;
  ExEnv::out0() << "  nedg = " << setw(10) << surf.nedge()
       << " bytes = "
       << setw(10) << surf.nedge() * sizeof(Edge)
       << endl;
  ExEnv::out0() << "  nver = " << setw(10) << surf.nvertex()
       << " bytes = "
       << setw(10) << surf.nvertex() * sizeof(Vertex)
       << endl;

  // compute normals if they weren't computed from the gradients
  if (!_volume->gradient_implemented()) {
      int i,j;
      // compute the normals as the average of the normals of
      // all the connected triangles
      for (i=0; i<surf.ntriangle(); i++) {
          Ref<Triangle> t = surf.triangle(i);
          SCVector3 tmp;
          SCVector3 BA = t->vertex(1)->point() - t->vertex(0)->point();
          SCVector3 CA = t->vertex(2)->point() - t->vertex(0)->point();
          SCVector3 N = BA.cross(CA);
          double n = N.norm();
          if (n < 1.0e-8) {
              tmp = 0.0;
            }
          else {
              n = 1.0/n;
              for (int j=0; j<3; j++) {
                  tmp[j] = - N[j]*n;
                }
            }
          for (j=0; j<3; j++) {
              int iv = surf.vertex_index(t->vertex(j));
              if (iv>=0) {
                  Ref<Vertex> v = surf.vertex(iv);
                  if (v->has_normal()) {
                      tmp += v->normal();
                    }
                  tmp.normalize();
                  v->set_normal(tmp);
                }
            }
        }
      // normalize all the normals
      for (i=0; i<surf.nvertex(); i++) {
          Ref<Vertex> v = surf.vertex(i);
          if (v->has_normal()) {
              SCVector3 n = v->normal();
              n.normalize();
              v->set_normal(n);
            }
          else {
              ExEnv::outn() << "ERROR: isosurf has a vertex without a triangle" << endl;
              abort();
            }
        }
    }
}

double
ImplicitSurfacePolygonizer::value_of_current(double x,double y,double z)
{
  current_x[0] = x; current_x[1] = y; current_x[2] = z;
  current->_volume->set_x(current_x);
  return current->_volume->value() - current->_value;
}

int
ImplicitSurfacePolygonizer::add_triangle_to_current(int i1, int i2, int i3,
                                                    sc::detail::VERTICES v)
{
  int oldlength = current->_tmp_vertices.size();
  current->_tmp_vertices.resize(v.count);
  for (int i=oldlength; i<v.count; i++) {
      SCVector3 newpoint;
      newpoint[0] = v.ptr[i].position.x;
      newpoint[1] = v.ptr[i].position.y;
      newpoint[2] = v.ptr[i].position.z;
      current->_volume->set_x(newpoint);
      SCVector3 normal;
      if (current->_volume->gradient_implemented()) {
          current->_volume->get_gradient(normal);
          normal.normalize();
          current->_tmp_vertices[i] = new Vertex(newpoint, normal);
        }
      else {
          current->_tmp_vertices[i] = new Vertex(newpoint);
        }
    }

  Ref<Vertex> v1 = current->_tmp_vertices[i1];
  Ref<Vertex> v2 = current->_tmp_vertices[i2];
  Ref<Vertex> v3 = current->_tmp_vertices[i3];
  
  static int tricnt = 0;
  if (++tricnt%100 == 0) {
      ExEnv::out0() << "adding triangle " << tricnt << endl;
      ExEnv::out0() << "  ntri = " << setw(10) << current->_surf->ntriangle()
           << " bytes = "
           << setw(10) << current->_surf->ntriangle() * sizeof(Triangle)
           << endl;
    }
  current->_surf->add_triangle(v1,v2,v3);

  return 1;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
