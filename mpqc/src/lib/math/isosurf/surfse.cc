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

#include <util/misc/formio.h>

#include <math/scmat/matrix.h>
#include <math/isosurf/surf.h>

void
TriangulatedSurface::remove_short_edges(double length_cutoff,
                                        const RefVolume &vol, double isoval)
{
  int j,k;
  AVLSet<RefTriangle>::iterator it,jt,kt;
  AVLSet<RefEdge>::iterator ie,je;

  int surface_was_completed = _completed_surface;

  if (!_completed_surface) {
      complete_ref_arrays();
    }
  else {
      clear_int_arrays();
    }

  _have_values = 0;
  _values.clear();
  
  int nvertex = _vertices.length();
  int nedge = _edges.length();
  int ntriangle = _triangles.length();

  if (_verbose) {
      cout << "TriangulatedSurface::remove_short_edges:" << endl
           << "initial: ";
      topology_info();
    }

  int deleted_edges_length;
  do {
      AVLSet<RefTriangle> deleted_triangles;
      AVLSet<RefEdge> deleted_edges;
      AVLSet<RefVertex> deleted_vertices;

      AVLSet<RefTriangle> new_triangles;
      AVLSet<RefEdge> new_edges;
      AVLSet<RefVertex> new_vertices;

      // a vertex to set-of-connected-triangles map
      AVLMap<RefVertex,AVLSet<RefTriangle> > connected_triangle_map;
      for (it = _triangles.begin(); it != _triangles.end(); it++) {
          RefTriangle tri = *it;
          for (j = 0; j<3; j++) {
              RefVertex v  = tri->vertex(j);
              connected_triangle_map[v].insert(tri);
            }
        }

      // a vertex to set-of-connected-edges map
      AVLMap<RefVertex,AVLSet<RefEdge> > connected_edge_map;
      for (ie = _edges.begin(); ie != _edges.end(); ie++) {
          RefEdge e = *ie;
          for (j = 0; j<2; j++) {
              RefVertex v = e->vertex(j);
              connected_edge_map[v].insert(e);
            }
        }
  
      for (ie = _edges.begin(); ie != _edges.end(); ie++) {
          RefEdge edge = *ie;
          double length = edge->straight_length();
          if (length < length_cutoff) {
              AVLSet<RefTriangle> connected_triangles;
              RefVertex v0 = edge->vertex(0);
              RefVertex v1 = edge->vertex(1);
              // = operator here causes problems
              connected_triangles |= connected_triangle_map[v0];
              connected_triangles |= connected_triangle_map[v1];
              int skip = 0;
              for (jt = connected_triangles.begin();
                   jt != connected_triangles.end();
                   jt++) {
                  RefTriangle tri = *jt;
                  if (tri->edge(0)->straight_length() < length
                      ||tri->edge(1)->straight_length() < length
                      ||tri->edge(2)->straight_length() < length
                      ||deleted_triangles.contains(tri)) {
                      skip = 1;
                      break;
                    }
                }
              if (skip) continue;
              deleted_triangles |= connected_triangles;
	      v0 = edge->vertex(0);
	      v1 = edge->vertex(1);
              deleted_vertices.insert(v0);
              deleted_vertices.insert(v1);
              // find all of the edges connected to the short edge
              // (including the short edge)
              AVLSet<RefEdge> connected_edges;
	      v0 = edge->vertex(0);
	      v1 = edge->vertex(1);
              connected_edges |= connected_edge_map[v0];
              connected_edges |= connected_edge_map[v1];
              deleted_edges |= connected_edges;
              // find the edges forming the perimeter of the deleted triangles
              // (these are used to form the new triangles)
              AVLSet<RefEdge> perimeter_edges;
              int embedded_triangle = 0;
              for (jt = connected_triangles.begin();
                   jt != connected_triangles.end();
                   jt++) {
                  RefTriangle tri = *jt;
                  for (j=0; j<3; j++) {
                      RefEdge e = tri->edge(j);
                      if (!connected_edges.contains(e)) {
                          // check to see if another triangle has claimed
                          // that this edge is a perimeter edge.  if so
                          // then this isn't a perimeter edge after all
                          // and it must be deleted.  this also implies
                          // that there is at least one triangle that
                          // has no perimeter edge.  these triangles and
                          // their nonperimeter vertices must be
                          // deleted.
                          RefEdge e = tri->edge(j);
                          if (perimeter_edges
                              .contains(e)) {
                              perimeter_edges.remove(e);
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
                  AVLSet<RefVertex> perimeter_vertices;
                  for (je = perimeter_edges.begin();
                       je != perimeter_edges.end();
                       je++) {
                      RefEdge e = *je;
                      for (j=0; j<2; j++) {
                          RefVertex v = e->vertex(j);
                          perimeter_vertices.insert(v);
                        }
                    }
                  // find the embedded_triangle
                  for (jt = connected_triangles.begin();
                       jt != connected_triangles.end();
                       jt++) {
                      RefTriangle tri = *jt;
                      // see if this triangle is embedded
                      for (j=0; j<3; j++) {
                          RefEdge e = tri->edge(j);
                          if (perimeter_edges.contains(e))
                              break;
                        }
                      // if embedded then delete the triangle's vertices
                      if (j==3) {
                          for (j=0; j<3; j++) {
                              RefVertex v = tri->vertex(j);
                              if (!perimeter_vertices.contains(v))
                                  deleted_vertices.insert(v);
                            }
                        }
                    }
                }
              // find a new point that replaces the deleted edge
              // (for now use one of the original, since it must lie on the
              // analytic surface)
              RefVertex replacement_vertex = edge->vertex(0);
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
              AVLMap<RefVertex,RefEdge> new_edge_map;
              for (je = perimeter_edges.begin(); je!=perimeter_edges.end();
                   je++) {
                  RefEdge e = *je;
                  for (k = 0; k<2; k++) {
                      RefVertex v = e->vertex(k);
                      if (!new_edge_map.contains(v)) {
                          RefEdge new_e = newEdge(
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
                  RefEdge e1 = *je;
                  RefVertex v0 = e1->vertex(0);
                  RefVertex v1 = e1->vertex(1);
                  RefEdge e2 = new_edge_map[v0];
                  RefEdge e3 = new_edge_map[v1];
                  if (e1.null() || e2.null() || e3.null()) {
                      cerr << "TriangulatedSurface::remove_short_edges: "
                           << "building new triangle but edges are null:"
                           << endl;
                      if (e1.null()) cerr << "  e1" << endl;
                      if (e2.null()) cerr << "  e2" << endl;
                      if (e3.null()) cerr << "  e3" << endl;
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
                  RefTriangle newtri(newTriangle(e1,e2,e3,orientation));
                  new_triangles.insert(newtri);
                }
              //cout << "WARNING: only one short edge removed" << endl;
              //break;
            }
        }

      _triangles -= deleted_triangles;
      _edges -= deleted_edges;
      _vertices -= deleted_vertices;

      _triangles |= new_triangles;
      _edges |= new_edges;
      _vertices |= new_vertices;

      if (_verbose) {
          topology_info();
        }

      deleted_edges_length = deleted_edges.length();
      //cout << "WARNING: one pass short edge removal" << endl;
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
      cout << "final: ";
      topology_info();
    }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
