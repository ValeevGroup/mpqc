
#include <math/scmat/matrix.h>
#include <math/isosurf/surf.h>
#include <math/isosurf/triRAVLMap.h>
#include <math/isosurf/veRAVLMap.h>
#include <math/isosurf/vtsRAVLMap.h>
#include <math/isosurf/edgeRAVLMap.h>
#include <math/isosurf/vertexAVLSet.h>

void
TriangulatedSurface::remove_short_edges(double length_cutoff)
{
  int j,k;
  Pix I,J,K;

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
      printf("TriangulatedSurface::remove_short_edges:\n");
      printf("initial: ");
      topology_info();
    }

  int deleted_edges_length;
  do {
      RefTriangleAVLSet deleted_triangles;
      RefEdgeAVLSet deleted_edges;
      RefVertexAVLSet deleted_vertices;

      RefTriangleAVLSet new_triangles;
      RefEdgeAVLSet new_edges;
      RefVertexAVLSet new_vertices;

      // a vertex to set-of-connected-triangles map
      RefTriangleAVLSet empty_triangle_set;
      RefVertexRefTriangleAVLSetRAVLMap connected_triangle_map(
          empty_triangle_set);
      for (I = _triangles.first(); I; _triangles.next(I)) {
          RefTriangle tri = _triangles(I);
          for (j = 0; j<3; j++) {
              RefVertex v  = tri->vertex(j);
              connected_triangle_map[v].add(tri);
            }
        }

      // a vertex to set-of-connected-edges map
      RefEdgeAVLSet empty_edge_set;
      RefVertexRefEdgeAVLSetRAVLMap connected_edge_map(empty_edge_set);
      for (I = _edges.first(); I; _edges.next(I)) {
          RefEdge e = _edges(I);
          for (j = 0; j<2; j++) {
              RefVertex v = e->vertex(j);
              connected_edge_map[v].add(e);
            }
        }
  
      for (I = _edges.first(); I; _edges.next(I)) {
          RefEdge edge = _edges(I);
          double length = edge->straight_length();
          if (length < length_cutoff) {
              RefTriangleAVLSet connected_triangles;
              // = operator here causes problems
              connected_triangles |= connected_triangle_map[edge->vertex(0)];
              connected_triangles |= connected_triangle_map[edge->vertex(1)];
              int skip = 0;
              for (J = connected_triangles.first();
                   J;
                   connected_triangles.next(J)) {
                  RefTriangle tri = connected_triangles(J);
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
              deleted_vertices.add(edge->vertex(0));
              deleted_vertices.add(edge->vertex(1));
              // find all of the edges connected to the short edge
              // (including the short edge)
              RefEdgeAVLSet connected_edges;
              connected_edges |= connected_edge_map[edge->vertex(0)];
              connected_edges |= connected_edge_map[edge->vertex(1)];
              deleted_edges |= connected_edges;
              // find the edges forming the perimeter of the deleted triangles
              // (these are used to form the new triangles)
              RefEdgeAVLSet perimeter_edges;
              int embedded_triangle = 0;
              for (J = connected_triangles.first();
                   J;
                   connected_triangles.next(J)) {
                  RefTriangle tri = connected_triangles(J);
                  for (j=0; j<3; j++) {
                      if (!connected_edges.contains(tri->edge(j))) {
                          // check to see if another triangle has claimed
                          // that this edge is a perimeter edge.  if so
                          // then this isn't a perimeter edge after all
                          // and it must be deleted.  this also implies
                          // that there is at least one triangle that
                          // has no perimeter edge.  these triangles and
                          // their nonperimeter vertices must be
                          // deleted.
                          if (perimeter_edges
                              .contains(tri->edge(j))) {
                              perimeter_edges.del(tri->edge(j));
                              deleted_edges.add(tri->edge(j));
                              embedded_triangle = 1;
                            }
                          else {
                              perimeter_edges.add(tri->edge(j));
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
                  RefVertexAVLSet perimeter_vertices;
                  for (J = perimeter_edges.first();
                       J;
                       perimeter_edges.next(J)) {
                      RefEdge e = perimeter_edges(J);
                      for (j=0; j<2; j++)
                          perimeter_vertices.add(e->vertex(j));
                    }
                  // find the embedded_triangle
                  for (J = connected_triangles.first();
                       J;
                       connected_triangles.next(J)) {
                      RefTriangle tri = connected_triangles(J);
                      // see if this triangle is embedded
                      for (j=0; j<3; j++) {
                          if (perimeter_edges.contains(tri->edge(j)))
                              break;
                        }
                      // if embedded then delete the triangle's vertices
                      if (j==3) {
                          for (j=0; j<3; j++) {
                              if (!perimeter_vertices.contains(tri->vertex(j)))
                                  deleted_vertices.add(tri->vertex(j));
                            }
                        }
                    }
                }
              // find a new point that replaces the deleted edge
              // (for now use one of the original, since it must lie on the
              // analytic surface)
              RefVertex replacement_vertex = edge->vertex(0);
              new_vertices.add(replacement_vertex);
              // for each vertex on the perimeter form a new edge to the
              // replacement vertex
              RefEdge empty_edge;
              RefVertexRefEdgeRAVLMap new_edge_map(empty_edge);
              for (J = perimeter_edges.first(); J; perimeter_edges.next(J)) {
                  RefEdge e = perimeter_edges(J);
                  for (k = 0; k<2; k++) {
                      if (!new_edge_map.contains(e->vertex(k))) {
                          RefEdge new_e = newEdge(
                              replacement_vertex,
                              e->vertex(k)
                              );
                          new_edge_map[e->vertex(k)] = new_e;
                          new_edges.add(new_e);
                        }
                    }
                }
              // for each edge on the perimeter form a new triangle with the
              // replacement vertex
              for (J = perimeter_edges.first(); J; perimeter_edges.next(J)) {
                  RefEdge e1 = perimeter_edges(J);
                  RefEdge e2 = new_edge_map[e1->vertex(0)];
                  RefEdge e3 = new_edge_map[e1->vertex(1)];
                  if (e1.null() || e2.null() || e3.null()) {
                      fprintf(stderr,
                              "TriangulatedSurface::remove_short_edges: "
                              "building new triangle but edges are null:\n");
                      if (e1.null()) fprintf(stderr,"  e1\n");
                      if (e2.null()) fprintf(stderr,"  e2\n");
                      if (e3.null()) fprintf(stderr,"  e3\n");
                      abort();
                    }
                  // Compute the correct orientation of e1 within the new
                  // triangle, by finding the orientation within the old
                  // triangle.
                  int orientation = 0;
                  for (K = connected_triangles.first();
                       K;
                       connected_triangles.next(K)) {
                      if (connected_triangles(K)->contains(e1)) {
                          orientation
                              = connected_triangles(K)->orientation(e1);
                          break;
                        }
                    }
                  RefTriangle newtri(newTriangle(e1,e2,e3,orientation));
                  new_triangles.add(newtri);
                }
              //printf("WARNING: only one short edge removed\n");
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
      //printf("WARNING: one pass short edge removal\n");
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
      printf("final: ");
      topology_info();
    }
}
