
#include <math.h>
#include <math/scmat/matrix.h>
#include <math/isosurf/surf.h>
#include <math/isosurf/triRAVLMap.h>
#include <math/isosurf/veRAVLMap.h>
#include <math/isosurf/vtsRAVLMap.h>
#include <math/isosurf/edgeRAVLMap.h>
#include <math/isosurf/vertexAVLSet.h>

#define VERBOSE 0

void
TriangulatedSurface::remove_slender_triangles(double height_cutoff)
{
  int i,j,k;
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

  if (VERBOSE) {
      printf("TriangulatedSurface::remove_slender_triangles:\n");
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
  
      for (I = _triangles.first(); I; _triangles.next(I)) {
          RefTriangle tri = _triangles(I);
          // find the heights of the vertices in tri
          double l[3], l2[3], h[3];
          for (j=0; j<3; j++) {
              l[j] = tri->edge(j)->length();
              if (l[j] <= 0.0) {
                  fprintf(stderr,"TriangulatedSurface::"
                          "remove_slender_triangles: bad edge length\n");
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
          if (h[hmin] < height_cutoff) {
              // find the vertex that gets eliminated
              RefVertex vertex;
              for (j=0; j<3; j++) {
                  if (tri->vertex(j) != tri->edge(hmin)->vertex(0)
                      &&tri->vertex(j) != tri->edge(hmin)->vertex(1)) {
                      vertex = tri->vertex(j);
                      break;
                    }
                }
              RefTriangleAVLSet connected_triangles;
              connected_triangles |= connected_triangle_map[vertex];

              // if one of the connected triangles is already being
              // deleted, save this one until the next pass
              int skip = 0;
              for (J = connected_triangles.first();
                   J;
                   connected_triangles.next(J)) {
                  RefTriangle tri = connected_triangles(J);
                  if (deleted_triangles.contains(tri)) {
                      skip = 1;
                      break;
                    }
                }
              if (skip) continue;
              
              deleted_triangles |= connected_triangles;
              deleted_vertices.add(vertex);
              // find all of the edges contained in the connected triangles
              RefEdgeAVLSet all_edges;
              for (J = connected_triangles.first();
                   J;
                   connected_triangles.next(J)) {
                  RefTriangle ctri = connected_triangles(J);
                  all_edges.add(ctri->edge(0));
                  all_edges.add(ctri->edge(1));
                  all_edges.add(ctri->edge(2));
                }
              // find all of the edges connected to the deleted vertex
              // (including the short edge)
              RefEdgeAVLSet connected_edges;
              connected_edges |= connected_edge_map[vertex];
              deleted_edges |= connected_edges;
              // find the edges forming the perimeter of the deleted triangles
              // (these are used to form the new triangles)
              RefEdgeAVLSet perimeter_edges;
              perimeter_edges |= all_edges;
              perimeter_edges -= connected_edges;
              // find a new point that replaces the deleted vertex
              // (for now use one of the original, since it must lie on the
              // analytic surface)
              RefVertex replacement_vertex;
              RefEdge short_edge;
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
              new_vertices.add(replacement_vertex);
              // for each vertex on the perimeter form a new edge to the
              // replacement vertex (unless the replacement vertex
              // is equal to the perimeter vertex)
              RefEdge empty_edge;
              RefVertexRefEdgeRAVLMap new_edge_map(empty_edge);
              for (J = perimeter_edges.first(); J; perimeter_edges.next(J)) {
                  RefEdge e = perimeter_edges(J);
                  for (k = 0; k<2; k++) {
                      if (e->vertex(k) == replacement_vertex) continue;
                      if (!new_edge_map.contains(e->vertex(k))) {
                          RefEdge new_e;
                          // if the edge already exists then use the
                          // existing edge
                          for (K = perimeter_edges.first();
                               K;
                               perimeter_edges.next(K)) {
                              RefEdge tmp = perimeter_edges(K);
                              if ((tmp->vertex(0) == replacement_vertex
                                   &&tmp->vertex(1) == e->vertex(k))
                                  ||(tmp->vertex(1) == replacement_vertex
                                     &&tmp->vertex(0) == e->vertex(k))) {
                                  new_e = tmp;
                                  break;
                                }
                            }
                          if (!K) {
                              new_e = newEdge(replacement_vertex,
                                              e->vertex(k));
                            }
                          new_edge_map[e->vertex(k)] = new_e;
                          new_edges.add(new_e);
                        }
                    }
                }
              // for each edge on the perimeter form a new triangle with the
              // replacement vertex (unless the edge contains the replacement
              // vertex)
              for (J = perimeter_edges.first(); J; perimeter_edges.next(J)) {
                  RefEdge e1 = perimeter_edges(J);
                  RefEdge e2 = new_edge_map[e1->vertex(0)];
                  RefEdge e3 = new_edge_map[e1->vertex(1)];
                  if (e1->vertex(0) == replacement_vertex
                      ||e1->vertex(1) == replacement_vertex) continue;
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
            }
        }

      _triangles -= deleted_triangles;
      _edges -= deleted_edges;
      _vertices -= deleted_vertices;

      _triangles |= new_triangles;
      _edges |= new_edges;
      _vertices |= new_vertices;

      deleted_edges_length = deleted_edges.length();
      //printf("WARNING: one pass short edge removal\n");
      //deleted_edges_length = 0; // do one pass
    } while(deleted_edges_length != 0);

  // fix the index maps
  _vertex_to_index.clear();
  _edge_to_index.clear();
  _triangle_to_index.clear();

  _index_to_vertex.clear();
  _index_to_edge.clear();
  _index_to_triangle.clear();

  int ne = _edges.length();
  int nv = _vertices.length();
  int nt = _triangles.length();

  for (i=0, I = _vertices.first(); I; i++, _vertices.next(I)) {
      _vertex_to_index[I] = i;
      _index_to_vertex[i] = I;
    }

  for (i=0, I = _edges.first(); I; i++, _edges.next(I)) {
      _edge_to_index[I] = i;
      _index_to_edge[i] = I;
    }

  for (i=0, I = _triangles.first(); I; i++, _triangles.next(I)) {
      _triangle_to_index[I] = i;
      _index_to_triangle[i] = I;
    }

  // restore the int arrays if they were there to begin with
  if (surface_was_completed) {
      complete_int_arrays();
      _completed_surface = 1;
    }

  if (VERBOSE) {
      printf("final: ");
      topology_info();
    }
}

