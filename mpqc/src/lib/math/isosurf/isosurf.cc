//
// isosurf.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@ca.sandia.gov>
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

#include <iostream.h>

#include <math/scmat/vector3.h>
#include <math/isosurf/isosurf.h>
#include <math/isosurf/implicit.h>
#include <math/isosurf/vtsRAVLMap.h>

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

ImplicitSurfacePolygonizer::ImplicitSurfacePolygonizer(const RefVolume&vol):
  _volume(vol)
{
}

ImplicitSurfacePolygonizer::~ImplicitSurfacePolygonizer()
{
}

ImplicitSurfacePolygonizer* ImplicitSurfacePolygonizer::current = 0;

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
  char *msg = polygonize(value_of_current, _resolution, bounds,
                         midpoint[0], midpoint[1], midpoint[2],
                         add_triangle_to_current, NOTET);
  current = 0;
  _surf = 0;
  if (msg) {
      cerr << "ImplicitSurfacePolygonizer::isosurface: failed: "
           << msg << endl;
      abort();
    }

  // Clean up temporaries.
  _tmp_vertices.clear();

  // finish the surface
  surf.complete_surface();

  // compute normals if they weren't computed from the gradients
  if (!_volume->gradient_implemented()) {
      int i;
      // make a list of what triangles are connected to each vertex
      RefTriangleAVLSet empty;
      RefVertexRefTriangleAVLSetRAVLMap vertex_to_triangles(empty);
      for (i=0; i<surf.ntriangle(); i++) {
          RefTriangle t = surf.triangle(i);
	  RefVertex v0 = t->vertex(0);
	  RefVertex v1 = t->vertex(1);
	  RefVertex v2 = t->vertex(2);
          vertex_to_triangles[v0].add(t);
          vertex_to_triangles[v1].add(t);
          vertex_to_triangles[v2].add(t);
        }
      for (i=0; i<surf.nvertex(); i++) {
          RefVertex v = surf.vertex(i);
          RefTriangleAVLSet triangles;
          triangles |= vertex_to_triangles[v];
          SCVector3 norm(0.0);
          SCVector3 tmp(0.0);
          for (Pix J = triangles.first(); J; triangles.next(J)) {
              RefTriangle t(triangles(J));
              // compute the normal to the surface
              // (using flat triangles)
              SCVector3 BA = t->vertex(1)->point() - t->vertex(0)->point();
              SCVector3 CA = t->vertex(2)->point() - t->vertex(0)->point();
              SCVector3 N = BA.cross(CA);
              double n = N.norm();
              if (n < 1.0e-15) {
                  tmp = 0.0;
                }
              else {
                  n = 1.0/n;
                  for (int j=0; j<3; j++) {
                      tmp[j] = - N[j]*n;
                    }
                }
              norm = norm + tmp;
            }
          norm.normalize();
          v->set_normal(norm);
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
                                                   VERTICES v)
{
  int i;
  for (i=current->_tmp_vertices.length(); i<v.count; i++) {
      SCVector3 newpoint;
      newpoint[0] = v.ptr[i].position.x;
      newpoint[1] = v.ptr[i].position.y;
      newpoint[2] = v.ptr[i].position.z;
      current->_volume->set_x(newpoint);
      SCVector3 normal;
      if (current->_volume->gradient_implemented()) {
          current->_volume->get_gradient(normal);
          normal.normalize();
        }
      current->_tmp_vertices.add(new Vertex(newpoint, normal));
    }

  RefVertex v1 = current->_tmp_vertices[i1];
  RefVertex v2 = current->_tmp_vertices[i2];
  RefVertex v3 = current->_tmp_vertices[i3];
  
  current->_surf->add_triangle(v1,v2,v3);

  return 1;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")
// End:
