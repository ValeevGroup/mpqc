//
// polysphere.cc
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
// This is loosely based on a program by Jon Leech.

#include <stdlib.h>
#include <iostream>
#include <math.h>

#include <util/render/polysphere.h>
#include <util/render/polygons.h>

using namespace std;
using namespace sc;

static inline void
switch2(int& i, int& j)
{
  int tmp = i;
  i = j;
  j = tmp;
}

class edge {
    int vertex[2];
  public:
    void set(int i, int j) {
        if (i == j) {
            ExEnv::errn() << "edge: bad nodes" << endl;
            abort();
          }
        vertex[0] = i; vertex[1] = j;
      }
    int& v(int i) { return vertex[i]; }
};

class triangle {
    int edge_[3];
    int orientation[3];
  public:
    void set(int e0, int o0,
             int e1, int o1,
             int e2, int o2) {
        edge_[0] = e0; orientation[0] = o0;
        edge_[1] = e1; orientation[1] = o1;
        edge_[2] = e2; orientation[2] = o2;
      }
    void set(int e0, int o0,
             int e1, int o1,
             int e2, int o2, edge* edges) {
        edge& E0 = edges[e0];
        edge& E1 = edges[e1];
        edge& E2 = edges[e2];

        if (  ((o0==0? E0.v(1): E0.v(0)) != (o1==0? E1.v(0): E1.v(1)))
            ||((o1==0? E1.v(1): E1.v(0)) != (o2==0? E2.v(0): E2.v(1)))
            ||((o2==0? E2.v(1): E2.v(0)) != (o0==0? E0.v(0): E0.v(1)))) {
            ExEnv::errn() << "triangle: bad edges or orientations" << endl;
            abort();
          }
        edge_[0] = e0; orientation[0] = o0;
        edge_[1] = e1; orientation[1] = o1;
        edge_[2] = e2; orientation[2] = o2;
      }
    int& e(int i) { return edge_[i]; }
    int& o(int i) { return orientation[i]; }
};

static void
subdivide(int level, int maxlevel,
          edge* edges, triangle* triangles,
          int nv, int ne, int nf, const Ref<RenderedPolygons>& poly)
{
  int i;
  
  if (level >= maxlevel) {
      // fill in poly
      for (i=0; i<nf; i++) {
          edge& e0 = edges[triangles[i].e(0)];
          edge& e1 = edges[triangles[i].e(1)];
          int v0 = (triangles[i].o(0)==0? e0.v(0): e0.v(1));
          int v1 = (triangles[i].o(0)==0? e0.v(1): e0.v(0));
          int v2 = (triangles[i].o(1)==0? e1.v(1): e1.v(0));
          if (v0 == v1 || v1 == v2 || v2 == v0) {
              ExEnv::errn() << "bad triangle" << endl;
              abort();
            }
          poly->set_face(i, v0, v1, v2);
        }
      return;
    }

  int nv2 = nv + ne;
  int ne2 = 3*nf + 2*ne;
  int nf2 = 4*nf;

  edge* newedges = new edge[ne2];
  triangle* newtriangles = new triangle[nf2];

  // split the edges and compute the new points
  int ipoint = nv;
  int inewedge = 0;
  for (i=0; i<ne; i++) {
      // split the edges
      newedges[inewedge].v(0) = edges[i].v(0);
      newedges[inewedge].v(1) = ipoint;
      inewedge++;
      newedges[inewedge].v(0) = ipoint;
      newedges[inewedge].v(1) = edges[i].v(1);
      inewedge++;

      // find the midpoint and normalize to unit length
      double v[3];
      v[0] = poly->vertex(edges[i].v(0), 0) + poly->vertex(edges[i].v(1), 0);
      v[1] = poly->vertex(edges[i].v(0), 1) + poly->vertex(edges[i].v(1), 1);
      v[2] = poly->vertex(edges[i].v(0), 2) + poly->vertex(edges[i].v(1), 2);
      double norm = 1.0/sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
      v[0] *= norm;
      v[1] *= norm;
      v[2] *= norm;
      poly->set_vertex(ipoint, v[0], v[1], v[2]);

      ipoint++;
    }

  // form the new triangles and the edges crossing the old triangles
  int inewtriangle = 0;
  for (i=0; i<nf; i++) {
      // the edges of the old triangles
      int e0 = triangles[i].e(0);
      int e1 = triangles[i].e(1);
      int e2 = triangles[i].e(2);

      // the orientations of the old triangles
      int o0 = triangles[i].o(0);
      int o1 = triangles[i].o(1);
      int o2 = triangles[i].o(2);

      // find the points in the centers of the old triangle's edges
      int p0 = nv + e0;
      int p1 = nv + e1;
      int p2 = nv + e2;

      // the edges obtained by dividing the old edges
      int e00 = 2*e0;
      int e01 = 2*e0 + 1;
      int e10 = 2*e1;
      int e11 = 2*e1 + 1;
      int e20 = 2*e2;
      int e21 = 2*e2 + 1;

      // connect the points to form new edges
      int ne0 = inewedge; newedges[inewedge++].set(p0, p1);
      int ne1 = inewedge; newedges[inewedge++].set(p1, p2);
      int ne2 = inewedge; newedges[inewedge++].set(p2, p0);

      // if original edges had reversed orientation, switch
      // ordering of split edges
      if (o0) switch2(e00, e01);
      if (o1) switch2(e10, e11);
      if (o2) switch2(e20, e21);

      // form the new triangles
      newtriangles[inewtriangle++].set(e00, o0,
                                       ne2, 1,
                                       e21, o2,
                                       newedges);
      newtriangles[inewtriangle++].set(e10, o1,
                                       ne0, 1,
                                       e01, o0,
                                       newedges);
      newtriangles[inewtriangle++].set(e20, o2,
                                       ne1, 1,
                                       e11, o1,
                                       newedges);
      newtriangles[inewtriangle++].set(ne0, 0,
                                       ne1, 0,
                                       ne2, 0,
                                       newedges);
    }

  subdivide(level+1,maxlevel,newedges,newtriangles,nv2,ne2,nf2,poly);
  delete[] newedges;
  delete[] newtriangles;
}

void
sc::polysphere(int level, const Ref<RenderedPolygons>& poly)
{
  int i;

  // compute the number of vertices, edges, and faces
  int nf = 8;
  int ne = 12;
  int nv = 6;
  for (i=2; i<=level; i++) {
      nv = nv + ne;
      ne = 3*nf + 2*ne;
      nf = 4*nf;
    }

  poly->initialize(nv,nf);

  // Fill in the first level, an octahedron.

  // the vertices
  poly->set_vertex(0, 1.0, 0.0, 0.0);
  poly->set_vertex(1,-1.0, 0.0, 0.0);
  poly->set_vertex(2, 0.0, 1.0, 0.0);
  poly->set_vertex(3, 0.0,-1.0, 0.0);
  poly->set_vertex(4, 0.0, 0.0, 1.0);
  poly->set_vertex(5, 0.0, 0.0,-1.0);

  edge *edges = new edge[12];
  edges[0].set(0, 5);
  edges[1].set(0, 2);
  edges[2].set(0, 4);
  edges[3].set(0, 3);
  edges[4].set(1, 5);
  edges[5].set(1, 2);
  edges[6].set(1, 4);
  edges[7].set(1, 3);
  edges[8].set(5, 2);
  edges[9].set(2, 4);
  edges[10].set(4, 3);
  edges[11].set(3, 5);

  triangle *triangles = new triangle[8];
  triangles[0].set(0, 0, 8, 0, 1, 1);
  triangles[1].set(1, 0, 9, 0, 2, 1);
  triangles[2].set(2, 0,10, 0, 3, 1);
  triangles[3].set(3, 0,11, 0, 0, 1);
  triangles[4].set(4, 0,11, 1, 7, 1);
  triangles[5].set(5, 0, 8, 1, 4, 1);
  triangles[6].set(6, 0, 9, 1, 5, 1);
  triangles[7].set(7, 0,10, 1, 6, 1);

  nf = 8;
  ne = 12;
  nv = 6;

  subdivide(1,level,edges,triangles,nv,ne,nf,poly);

  delete[] edges;
  delete[] triangles;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
