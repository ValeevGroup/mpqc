//
// surfse.cc
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

#include <algorithm>

#include <util/misc/formio.h>

#include <math/scmat/matrix.h>
#include <math/isosurf/surf.h>

using namespace std;
using namespace sc;

void
TriangulatedSurface::remove_short_edges(double length_cutoff,
                                        const Ref<Volume> &vol, double isoval)
{
  int j,k;
  std::set<Ref<Triangle> >::iterator it,jt,kt;
  std::set<Ref<Edge> >::iterator ie,je;

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
      ExEnv::outn() << "TriangulatedSurface::remove_short_edges:" << endl
           << "initial: ";
      topology_info();
    }

  int deleted_edges_length;
  do {
      std::set<Ref<Triangle> > deleted_triangles;
      std::set<Ref<Edge> > deleted_edges;
      std::set<Ref<Vertex> > deleted_vertices;

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
  
      for (ie = _edges.begin(); ie != _edges.end(); ie++) {
          Ref<Edge> edge = *ie;
          double length = edge->straight_length();
          if (length < length_cutoff) {
              std::set<Ref<Triangle> > connected_triangles;
              Ref<Vertex> v0 = edge->vertex(0);
              Ref<Vertex> v1 = edge->vertex(1);
              connected_triangles.insert(connected_triangle_map[v0].begin(),
                                         connected_triangle_map[v0].end());
              connected_triangles.insert(connected_triangle_map[v1].begin(),
                                         connected_triangle_map[v1].end());
              int skip = 0;
              for (jt = connected_triangles.begin();
                   jt != connected_triangles.end();
                   jt++) {
                  Ref<Triangle> tri = *jt;
                  if (tri->edge(0)->straight_length() < length
                      ||tri->edge(1)->straight_length() < length
                      ||tri->edge(2)->straight_length() < length
                      ||deleted_triangles.find(tri)!=deleted_triangles.end()) {
                      skip = 1;
                      break;
                    }
                }
              if (skip) continue;
              deleted_triangles.insert(connected_triangles.begin(),
                                       connected_triangles.end());
	      v0 = edge->vertex(0);
	      v1 = edge->vertex(1);
              deleted_vertices.insert(v0);
              deleted_vertices.insert(v1);
              // find all of the edges connected to the short edge
              // (including the short edge)
              std::set<Ref<Edge> > connected_edges;
	      v0 = edge->vertex(0);
	      v1 = edge->vertex(1);
              connected_edges.insert(connected_edge_map[v0].begin(),
                                     connected_edge_map[v0].end());
              connected_edges.insert(connected_edge_map[v1].begin(),
                                     connected_edge_map[v1].end());
              deleted_edges.insert(connected_edges.begin(),
                                   connected_edges.end());
              // find the edges forming the perimeter of the deleted triangles
              // (these are used to form the new triangles)
              std::set<Ref<Edge> > perimeter_edges;
              int embedded_triangle = 0;
              for (jt = connected_triangles.begin();
                   jt != connected_triangles.end();
                   jt++) {
                  Ref<Triangle> tri = *jt;
                  for (j=0; j<3; j++) {
                      Ref<Edge> e = tri->edge(j);
                      if (connected_edges.find(e) == connected_edges.end()) {
                          // check to see if another triangle has claimed
                          // that this edge is a perimeter edge.  if so
                          // then this isn't a perimeter edge after all
                          // and it must be deleted.  this also implies
                          // that there is at least one triangle that
                          // has no perimeter edge.  these triangles and
                          // their nonperimeter vertices must be
                          // deleted.
                          Ref<Edge> e = tri->edge(j);
                          if (perimeter_edges.find(e)
                              != perimeter_edges.end()) {
                              perimeter_edges.erase(e);
                              deleted_edges.insert(e);
                              embedded_triangle = 1;
                            }
                          else {
                              perimeter_edges.insert(e);
                            }
                        }
                    }
                }
              // if a triangle is embedded make sure its vertices are
              // all deleted.  (the triangle itself should already be
              // connected and thus deleted).
              if (embedded_triangle) {
                  // make a list of vertices on the perimeter (so i
                  // don't delete them
                  std::set<Ref<Vertex> > perimeter_vertices;
                  for (je = perimeter_edges.begin();
                       je != perimeter_edges.end();
                       je++) {
                      Ref<Edge> e = *je;
                      for (j=0; j<2; j++) {
                          Ref<Vertex> v = e->vertex(j);
                          perimeter_vertices.insert(v);
                        }
                    }
                  // find the embedded_triangle
                  for (jt = connected_triangles.begin();
                       jt != connected_triangles.end();
                       jt++) {
                      Ref<Triangle> tri = *jt;
                      // see if this triangle is embedded
                      for (j=0; j<3; j++) {
                          Ref<Edge> e = tri->edge(j);
                          if (perimeter_edges.find(e) != perimeter_edges.end())
                              break;
                        }
                      // if embedded then delete the triangle's vertices
                      if (j==3) {
                          for (j=0; j<3; j++) {
                              Ref<Vertex> v = tri->vertex(j);
                              if (perimeter_vertices.find(v) == perimeter_vertices.end())
                                  deleted_vertices.insert(v);
                            }
                        }
                    }
                }
              // find a new point that replaces the deleted edge
              // (for now use one of the original, since it must lie on the
              // analytic surface)
              Ref<Vertex> replacement_vertex = edge->vertex(0);
              // however, if we have a volume, find a new vertex on
              // the analytic surface near the center of the edge
              if (vol.nonnull()) {
                  SCVector3 point, norm;
                  int hn = edge->interpolate(0.5,point,norm,vol,isoval);
                  replacement_vertex = new Vertex(point);
                  if (hn) replacement_vertex->set_normal(norm);
                }
              new_vertices.insert(replacement_vertex);
              // for each vertex on the perimeter form a new edge to the
              // replacement vertex
              std::map<Ref<Vertex>,Ref<Edge> > new_edge_map;
              for (je = perimeter_edges.begin(); je!=perimeter_edges.end();
                   je++) {
                  Ref<Edge> e = *je;
                  for (k = 0; k<2; k++) {
                      Ref<Vertex> v = e->vertex(k);
                      if (new_edge_map.find(v) == new_edge_map.end()) {
                          Ref<Edge> new_e = newEdge(
                              replacement_vertex,
                              v
                              );
                          new_edge_map[v] = new_e;
                          new_edges.insert(new_e);
                        }
                    }
                }
              // for each edge on the perimeter form a new triangle with the
              // replacement vertex
              for (je = perimeter_edges.begin(); je != perimeter_edges.end();
                   je++) {
                  Ref<Edge> e1 = *je;
                  Ref<Vertex> v0 = e1->vertex(0);
                  Ref<Vertex> v1 = e1->vertex(1);
                  Ref<Edge> e2 = new_edge_map[v0];
                  Ref<Edge> e3 = new_edge_map[v1];
                  if (e1.null() || e2.null() || e3.null()) {
                      ExEnv::errn() << "TriangulatedSurface::remove_short_edges: "
                           << "building new triangle but edges are null:"
                           << endl;
                      if (e1.null()) ExEnv::errn() << "  e1" << endl;
                      if (e2.null()) ExEnv::errn() << "  e2" << endl;
                      if (e3.null()) ExEnv::errn() << "  e3" << endl;
                      abort();
                    }
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
              //ExEnv::outn() << "WARNING: only one short edge removed" << endl;
              //break;
            }
        }

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
          topology_info();
        }

      deleted_edges_length = deleted_edges.size();
      //ExEnv::outn() << "WARNING: one pass short edge removal" << endl;
      //deleted_edges_length = 0; // do one pass
    } while(deleted_edges_length != 0);

  // fix the index maps
  recompute_index_maps();

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
