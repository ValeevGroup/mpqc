
#include <math/scmat/matrix.h>
#include <math/isosurf/surf.h>
#include <math/isosurf/triRAVLMap.h>
#include <math/isosurf/veRAVLMap.h>
#include <math/isosurf/vtsRAVLMap.h>
#include <math/isosurf/edgeRAVLMap.h>
#include <math/isosurf/vertexAVLSet.h>

void
TriangulatedSurface::fix_orientation()
{
  int i,j;
  Pix I;
  int nflip = 0;
  
  int ne = nedge();
  int ntri = ntriangle();
  RefTriangle *edge_to_triangle0;
  RefTriangle *edge_to_triangle1;
  edge_to_triangle0 = new RefTriangle[ne];
  edge_to_triangle1 = new RefTriangle[ne];

  for (I = _triangles.first(); I; _triangles.next(I)) {
      RefTriangle tri = _triangles(I);
      for (j=0; j<3; j++) {
          RefEdge e = tri->edge(j);
          Pix ix = _edges.seek(e);
          int e_index = _edge_to_index[ix];
          if (edge_to_triangle0[e_index].null()) {
              edge_to_triangle0[e_index] = tri;
            }
          else if (edge_to_triangle1[e_index].null()) {
              edge_to_triangle1[e_index] = tri;
            }
          else {
              fprintf(stderr,"TriangulatedSurface::fix_orientation:"
                      " more than two triangles to an edge\n");
              abort();
            }
        }
    }

  RefTriangleAVLSet unfixed;
  RefTriangleAVLSet fixed;
  RefTriangleAVLSet finished;

  unfixed |= _triangles;

  while (unfixed.length()) {
      // define unfixed.first()'s orientation to be the fixed orientation
      Pix first = unfixed.first();
      fixed.add(unfixed(first));
      unfixed.del(unfixed(first));
      while (fixed.length()) {
          RefTriangle tri = fixed(fixed.first());
          // make all neighbors of tri oriented the same as tri
          for (i=0; i<3; i++) {
              RefEdge e = tri->edge(i);
              int e_index = _edge_to_index[_edges.seek(e)];
              RefTriangle othertri;
              if (edge_to_triangle0[e_index] == tri) {
                  othertri = edge_to_triangle1[e_index];
                }
              else {
                  othertri = edge_to_triangle0[e_index];
                }
              for (j=0; j<3; j++) {
                  if (othertri->edge(j) == e) break;
                }
              if (j == 3) {
                  fprintf(stderr,"TriangulatedSurface::fix_orientation: "
                          " edge_to_triangle wrong\n");
                  abort();
                }
              if (tri->orientation(i) == othertri->orientation(j)) {
                  if (unfixed.contains(othertri)) {
                      unfixed.del(othertri);
                      fixed.add(othertri);
                      othertri->flip();
                      nflip++;
                    }
                  else {
                      fprintf(stderr,"TriangulatedSurface::fix_orientation:"
                              " tried to flip a fixed triangle\n");
                      abort();
                    }
                }
              else if (unfixed.contains(othertri)) {
                  unfixed.del(othertri);
                  fixed.add(othertri);
                }
            }
          // by this point all of tri's neighbors have had their
          // orientation fixed to match that of tri
          fixed.del(tri);
          finished.add(tri);
        }
    }

  printf("%d out of %d triangles were flipped\n", nflip, ntri);

  // just in case
  for (i=0; i<ne; i++) {
      edge_to_triangle0[i] = 0;
      edge_to_triangle1[i] = 0;
    }

  delete[] edge_to_triangle0;
  delete[] edge_to_triangle1;
}
