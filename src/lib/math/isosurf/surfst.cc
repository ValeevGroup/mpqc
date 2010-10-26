//
// surfst.cc
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

#include <stdio.h>
#include <math.h>

#include <util/misc/formio.h>

#include <math/scmat/matrix.h>
#include <math/isosurf/surf.h>

using namespace std;
using namespace sc;

#ifndef WRITE_OOGL // this is useful for debugging this routine
#define WRITE_OOGL 1
#endif

#if WRITE_OOGL
#include <util/render/oogl.h>
#include <util/render/polygons.h>
#endif

void
TriangulatedSurface::remove_slender_triangles(
    int remove_slender, double height_cutoff,
    int remove_small, double area_cutoff,
    const Ref<Volume> &vol, double isoval)
{
  int i,j,k;
  std::set<Ref<Triangle> >::iterator it,jt,kt;
  std::set<Ref<Edge> >::iterator ie,je,ke;
  std::set<Ref<Vertex> >::iterator iv,jv;

  int surface_was_completed = _completed_surface;

  if (!_completed_surface) {
      complete_ref_arrays();
    }
  else {
      clear_int_arrays();
    }

  _have_values = 0;
  _values.clear();
  
  if (_verbose) {
      ExEnv::outn() << "TriangulatedSurface::remove_slender_triangles:" << endl
           << "initial: ";
      topology_info();
    }

#if WRITE_OOGL
  if (_debug) {
      render(new OOGLRender("surfstinit.oogl"));
    }
#endif

  int deleted_edges_length;
  do {
      std::set<Ref<Triangle> > deleted_triangles;
      std::set<Ref<Edge> > deleted_edges;
      std::set<Ref<Vertex> > deleted_vertices;

      std::set<Ref<Vertex> > vertices_of_deleted_triangles;

      std::set<Ref<Triangle> > new_triangles;
      std::set<Ref<Edge> > new_edges;
      std::set<Ref<Vertex> > new_vertices;

      // a vertex to set-of-connected-triangles map
      std::map<Ref<Vertex>,std::set<Ref<Triangle> > > connected_triangle_map;
      for (it = _triangles.begin(); it != _triangles.end(); it++) {
          Ref<Triangle> tri = *it;
          for (j = 0; j<3; j++) {
              Ref<Vertex> v  = tri->vertex(j);
              connected_triangle_map[v].insert(tri);
            }
        }

      // a vertex to set-of-connected-edges map
      std::map<Ref<Vertex>,std::set<Ref<Edge> > > connected_edge_map;
      for (ie = _edges.begin(); ie != _edges.end(); ie++) {
          Ref<Edge> e = *ie;
          for (j = 0; j<2; j++) {
              Ref<Vertex> v = e->vertex(j);
              connected_edge_map[v].insert(e);
            }
        }
  
      for (it = _triangles.begin(); it != _triangles.end(); it++) {
          Ref<Triangle> tri = *it;
          // find the heights of the vertices in tri
          double l[3], l2[3], h[3];
          for (j=0; j<3; j++) {
              l[j] = tri->edge(j)->straight_length();
              if (l[j] <= 0.0) {
                  ExEnv::errn() << "TriangulatedSurface::"
                       << "remove_slender_triangles: bad edge length"
                       << endl;
                  abort();
                }
              l2[j] = l[j]*l[j];
            }
          double y = 2.0*(l2[0]*l2[1]+l2[0]*l2[2]+l2[1]*l2[2])
                     - l2[0]*l2[0] - l2[1]*l2[1] - l2[2]*l2[2];
          if (y < 0.0) y = 0.0;
          double x = 0.5*sqrt(y);
          for (j=0; j<3; j++) h[j] = x/l[j];

          // find the shortest height
          int hmin;
          if (h[0] < h[1]) hmin = 0;
          else hmin = 1;
          if (h[2] < h[hmin]) hmin = 2;

          // see if the shortest height is below the cutoff
          if (remove_slender && h[hmin] < height_cutoff) {
              // find the vertex that gets eliminated
              Ref<Vertex> vertex;
              for (j=0; j<3; j++) {
                  if (tri->vertex(j) != tri->edge(hmin)->vertex(0)
                      &&tri->vertex(j) != tri->edge(hmin)->vertex(1)) {
                      vertex = tri->vertex(j);
                      break;
                    }
                }
              std::set<Ref<Triangle> > connected_triangles;
              connected_triangles = connected_triangle_map[vertex];

              // if one of the connected triangles has a vertex
              // in a deleted triangle, save this one until the
              // next pass
              int skip = 0;
              for (jt = connected_triangles.begin();
                   jt != connected_triangles.end();
                   jt++) {
                  Ref<Triangle> tri = *jt;
                  for (j=0; j<3; j++) {
                      Ref<Vertex> v = tri->vertex(j);
                      if (vertices_of_deleted_triangles.find(v)
                          != vertices_of_deleted_triangles.end()) {
                          skip = 1;
                          break;
                        }
                    }
                  if (skip) break;
                }
              if (skip) continue;

              // find all of the edges contained in the connected triangles
              std::set<Ref<Edge> > all_edges;
              for (jt = connected_triangles.begin();
                   jt != connected_triangles.end();
                   jt++) {
                  Ref<Triangle> ctri = *jt;
                  Ref<Edge> e0 = ctri->edge(0);
                  Ref<Edge> e1 = ctri->edge(1);
                  Ref<Edge> e2 = ctri->edge(2);
                  all_edges.insert(e0);
                  all_edges.insert(e1);
                  all_edges.insert(e2);
                }
              // find all of the edges connected to the deleted vertex
              // (including the short edge)
              std::set<Ref<Edge> > connected_edges;
              connected_edges = connected_edge_map[vertex];
              // find the edges forming the perimeter of the deleted triangles
              // (these are used to form the new triangles)
              std::set<Ref<Edge> > perimeter_edges;
              perimeter_edges = all_edges;
              erase_elements_by_value(perimeter_edges,
                                      connected_edges.begin(),
                                      connected_edges.end());

              // If deleting this point causes a flattened piece of
              // surface, reject it.  This is tested by checking for
              // triangles that have all vertices contained in the set
              // of vertices connected to the deleted vertex.
              std::set<Ref<Vertex> > connected_vertices;
              for (je = perimeter_edges.begin(); je != perimeter_edges.end();
                   je++) {
                  Ref<Edge> e = *je;
                  Ref<Vertex> v0 = e->vertex(0);
                  Ref<Vertex> v1 = e->vertex(1);
                  connected_vertices.insert(v0);
                  connected_vertices.insert(v1);
                }
              std::set<Ref<Triangle> > triangles_connected_to_perimeter;
              for (jv = connected_vertices.begin();
                   jv != connected_vertices.end();
                   jv++) {
                  triangles_connected_to_perimeter.insert(
                      connected_triangle_map[*jv].begin(),
                      connected_triangle_map[*jv].end());
                }
              for (jt = triangles_connected_to_perimeter.begin();
                   jt != triangles_connected_to_perimeter.end();
                   jt++) {
                  Ref<Triangle> t = *jt;
                  Ref<Vertex> v0 = t->vertex(0);
                  Ref<Vertex> v1 = t->vertex(1);
                  Ref<Vertex> v2 = t->vertex(2);
                  if (connected_vertices.find(v0)!=connected_vertices.end()
                      &&connected_vertices.find(v1)!=connected_vertices.end()
                      &&connected_vertices.find(v2)!=connected_vertices.end())
                      {
                      skip = 1;
                      break;
                    }
                }
              if (skip) {
                  continue;
                }

              deleted_triangles.insert(connected_triangles.begin(),
                                       connected_triangles.end());
              deleted_vertices.insert(vertex);
              deleted_edges.insert(connected_edges.begin(),
                                   connected_edges.end());

              for (jt = deleted_triangles.begin();
                   jt != deleted_triangles.end();
                   jt++) {
                  Ref<Triangle> t = *jt;
                  for (j=0; j<2; j++) {
                      Ref<Vertex> v = t->vertex(j);
                      vertices_of_deleted_triangles.insert(v);
                    }
                }

              // find a new point that replaces the deleted vertex
              // (for now use one of the original, since it must lie on the
              // analytic surface)
              Ref<Vertex> replacement_vertex;
              Ref<Edge> short_edge;
              if (hmin==0) {
                  if (l[1] < l[2]) short_edge = tri->edge(1);
                  else short_edge = tri->edge(2);
                }
              else if (hmin==1) {
                  if (l[0] < l[2]) short_edge = tri->edge(0);
                  else short_edge = tri->edge(2);
                }
              else {
                  if (l[0] < l[1]) short_edge = tri->edge(0);
                  else short_edge = tri->edge(1);
                }
              if (short_edge->vertex(0) == vertex) {
                  replacement_vertex = short_edge->vertex(1);
                }
              else {
                  replacement_vertex = short_edge->vertex(0);
                }
              new_vertices.insert(replacement_vertex);
              // for each vertex on the perimeter form a new edge to the
              // replacement vertex (unless the replacement vertex
              // is equal to the perimeter vertex)
              std::map<Ref<Vertex>,Ref<Edge> > new_edge_map;
              for (je = perimeter_edges.begin(); je != perimeter_edges.end();
                   je++) {
                  Ref<Edge> e = *je;
                  for (k = 0; k<2; k++) {
                      Ref<Vertex> v = e->vertex(k);
                      if (v == replacement_vertex) continue;
                      if (new_edge_map.find(v) == new_edge_map.end()) {
                          Ref<Edge> new_e;
                          // if the edge already exists then use the
                          // existing edge
                          for (ke = perimeter_edges.begin();
                               ke != perimeter_edges.end();
                               ke++) {
                              Ref<Edge> tmp = *ke;
                              if ((tmp->vertex(0) == replacement_vertex
                                   &&tmp->vertex(1) == v)
                                  ||(tmp->vertex(1) == replacement_vertex
                                     &&tmp->vertex(0) == v)) {
                                  new_e = tmp;
                                  break;
                                }
                            }
                          if (ke == perimeter_edges.end()) {
                              new_e = newEdge(replacement_vertex,
                                              v);
                            }
                          new_edge_map[v] = new_e;
                          new_edges.insert(new_e);
                        }
                    }
                }
              // for each edge on the perimeter form a new triangle with the
              // replacement vertex (unless the edge contains the replacement
              // vertex)
              for (je = perimeter_edges.begin(); je != perimeter_edges.end();
                   je++) {
                  Ref<Edge> e1 = *je;
                  Ref<Vertex> v0 = e1->vertex(0);
                  Ref<Vertex> v1 = e1->vertex(1);
                  Ref<Edge> e2 = new_edge_map[v0];
                  Ref<Edge> e3 = new_edge_map[v1];
                  if (v0 == replacement_vertex
                      || v1 == replacement_vertex) continue;
                  // Compute the correct orientation of e1 within the new
                  // triangle, by finding the orientation within the old
                  // triangle.
                  int orientation = 0;
                  for (kt = connected_triangles.begin();
                       kt != connected_triangles.end();
                       kt++) {
                      if ((*kt)->contains(e1)) {
                          orientation
                              = (*kt)->orientation(e1);
                          break;
                        }
                    }
                  Ref<Triangle> newtri(newTriangle(e1,e2,e3,orientation));
                  new_triangles.insert(newtri);
                }
            }
        }

#if WRITE_OOGL
      if (_debug) {
          char filename[100];
          static int pass = 0;
          sprintf(filename, "surfst%04d.oogl", pass);
          ExEnv::outn() << scprintf("PASS = %04d\n", pass);
          Ref<Render> render = new OOGLRender(filename);
          Ref<RenderedPolygons> poly = new RenderedPolygons;
          poly->initialize(_vertices.size(), _triangles.size(),
                           RenderedPolygons::Vertex);
          // the number of triangles and edges touching a vertex
          int *n_triangle = new int[_vertices.size()];
          int *n_edge = new int[_vertices.size()];
          memset(n_triangle,0,sizeof(int)*_vertices.size());
          memset(n_edge,0,sizeof(int)*_vertices.size());
          std::set<Ref<Triangle> >::iterator it;
          std::set<Ref<Edge> >::iterator ie;
          std::set<Ref<Vertex> >::iterator iv;
          std::map<Ref<Vertex>, int> vertex_to_index;
          int i = 0;
          for (iv = _vertices.begin(); iv != _vertices.end(); iv++, i++) {
              Ref<Vertex> v = *iv;
              vertex_to_index[v] = i;
              poly->set_vertex(i,
                               v->point()[0],
                               v->point()[1],
                               v->point()[2]);
              if (deleted_vertices.find(v) != deleted_vertices.end()) {
                  poly->set_vertex_rgb(i, 1.0, 0.0, 0.0);
                }
              else {
                  poly->set_vertex_rgb(i, 0.3, 0.3, 0.3);
                }
            }
          i = 0;
          for (it = _triangles.begin(); it != _triangles.end(); it++, i++) {
              Ref<Triangle> t = *it;
              int i0 = vertex_to_index[t->vertex(0)];
              int i1 = vertex_to_index[t->vertex(1)];
              int i2 = vertex_to_index[t->vertex(2)];
              n_triangle[i0]++;
              n_triangle[i1]++;
              n_triangle[i2]++;
              poly->set_face(i,i0,i1,i2);
            }
          for (ie = _edges.begin(); ie != _edges.end(); ie++, i++) {
              Ref<Edge> e = *ie;
              int i0 = vertex_to_index[e->vertex(0)];
              int i1 = vertex_to_index[e->vertex(1)];
              n_edge[i0]++;
              n_edge[i1]++;
            }
          i = 0;
          for (iv = _vertices.begin(); iv != _vertices.end(); iv++, i++) {
              Ref<Vertex> v = *iv;
              if (n_triangle[i] != n_edge[i]) {
                  ExEnv::outn() << "found bad vertex"
                       << " nedge = " << n_edge[i]
                       << " ntriangle = " << n_triangle[i]
                       << endl;
                  if (deleted_vertices.find(v) != deleted_vertices.end()) {
                      poly->set_vertex_rgb(i, 1.0, 1.0, 0.0);
                    }
                  else {
                      poly->set_vertex_rgb(i, 0.0, 1.0, 0.0);
                    }
                }
            }
          render->render(poly.pointer());
          pass++;
          delete[] n_triangle;
          delete[] n_edge;
        }
#endif

      erase_elements_by_value(_triangles,
                              deleted_triangles.begin(),
                              deleted_triangles.end());
      erase_elements_by_value(_edges,
                              deleted_edges.begin(),
                              deleted_edges.end());
      erase_elements_by_value(_vertices,
                              deleted_vertices.begin(),
                              deleted_vertices.end());

      _triangles.insert(new_triangles.begin(), new_triangles.end());
      _edges.insert(new_edges.begin(), new_edges.end());
      _vertices.insert(new_vertices.begin(), new_vertices.end());

      if (_verbose) {
          ExEnv::outn() << "intermediate: ";
          topology_info();
        }

      deleted_edges_length = deleted_edges.size();
    } while(deleted_edges_length != 0);

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

  // restore the int arrays if they were there to begin with
  if (surface_was_completed) {
      complete_int_arrays();
      _completed_surface = 1;
    }

  if (_verbose) {
      ExEnv::outn() << "final: ";
      topology_info();
    }
}


/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
