
#include <math.h>
#include <math/scmat/matrix.h>
#include <math/isosurf/surf.h>
#include <math/isosurf/triRAVLMap.h>
#include <math/isosurf/veRAVLMap.h>
#include <math/isosurf/vtsRAVLMap.h>
#include <math/isosurf/edgeRAVLMap.h>
#include <math/isosurf/vertexAVLSet.h>

#ifndef WRITE_OOGL // this is useful for debugging this routine
#define WRITE_OOGL 0
#endif

#if WRITE_OOGL
#include <util/container/pixintRAVLMap.h>
#include <util/render/oogl.h>
#include <util/render/polygons.h>
#endif

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

  if (_verbose) {
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
              l[j] = tri->edge(j)->straight_length();
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
              // find the edges forming the perimeter of the deleted triangles
              // (these are used to form the new triangles)
              RefEdgeAVLSet perimeter_edges;
              perimeter_edges |= all_edges;
              perimeter_edges -= connected_edges;

              // If deleting this point causes a flattened piece of
              // surface, reject it.  This is tested by checking for
              // triangles that have all vertices contained in the set
              // of vertices connected to the deleted vertex.
              RefVertexAVLSet connected_vertices;
              for (J = perimeter_edges.first(); J; perimeter_edges.next(J)) {
                  RefEdge e = perimeter_edges(J);
                  connected_vertices.add(e->vertex(0));
                  connected_vertices.add(e->vertex(1));
                }
              RefTriangleAVLSet triangles_connected_to_perimeter;
              for (J = connected_vertices.first();
                   J;
                   connected_vertices.next(J)) {
                  triangles_connected_to_perimeter
                      |= connected_triangle_map[connected_vertices(J)];
                }
              for (J = triangles_connected_to_perimeter.first();
                   J;
                   triangles_connected_to_perimeter.next(J)) {
                  RefTriangle t = triangles_connected_to_perimeter(J);
                  if (connected_vertices.seek(t->vertex(0))
                      &&connected_vertices.seek(t->vertex(1))
                      &&connected_vertices.seek(t->vertex(2))) {
                      skip = 1;
                      break;
                    }
                }
              if (skip) {
                  continue;
                }

              deleted_triangles |= connected_triangles;
              deleted_vertices.add(vertex);
              deleted_edges |= connected_edges;

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

#if WRITE_OOGL
        {
          char filename[100];
          static int pass = 0;
          sprintf(filename, "surfst%04d.oogl", pass);
          RefRender render = new OOGLRender(filename);
          RefRenderedPolygons poly = new RenderedPolygons;
          poly->initialize(_vertices.length(), _triangles.length(),
                           RenderedPolygons::Vertex);
          Pix I;
          PixintRAVLMap pix_to_index(0);
          int i = 0;
          for (I = _vertices.first(); I; _vertices.next(I), i++) {
              RefVertex v = _vertices(I);
              pix_to_index[(Pix)v.pointer()] = i;
              poly->set_vertex(i,
                               v->point()->get_element(0),
                               v->point()->get_element(1),
                               v->point()->get_element(2));
              if (deleted_vertices.seek(v)) {
                  poly->set_vertex_rgb(i, 1.0, 0.0, 0.0);
                }
              else {
                  poly->set_vertex_rgb(i, 0.3, 0.3, 0.3);
                }
            }
          i = 0;
          for (I = _triangles.first(); I; _triangles.next(I), i++) {
              RefTriangle t = _triangles(I);
              poly->set_face(i,
                             pix_to_index[(Pix)t->vertex(0).pointer()],
                             pix_to_index[(Pix)t->vertex(1).pointer()],
                             pix_to_index[(Pix)t->vertex(2).pointer()]);
            }
          if (_verbose) {
              printf("PASS = %04d: ", pass);
            }
          else {
              printf("PASS = %04d\n", pass);
            }
          render->render(poly);
          pass++;
        }
#endif

      _triangles -= deleted_triangles;
      _edges -= deleted_edges;
      _vertices -= deleted_vertices;

      _triangles |= new_triangles;
      _edges |= new_edges;
      _vertices |= new_vertices;

      if (_verbose) {
          printf("intermediate: ");
          topology_info();
        }

      deleted_edges_length = deleted_edges.length();
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

  if (_verbose) {
      printf("final: ");
      topology_info();
    }
}

