//
// surfor.cc
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

using namespace std;
using namespace sc;

void
TriangulatedSurface::fix_orientation()
{
  int i,j;
  std::set<Ref<Triangle> >::iterator I;
  int nflip = 0;
  
  int ne = nedge();
  int ntri = ntriangle();
  Ref<Triangle> *edge_to_triangle0;
  Ref<Triangle> *edge_to_triangle1;
  edge_to_triangle0 = new Ref<Triangle>[ne];
  edge_to_triangle1 = new Ref<Triangle>[ne];

  for (I = _triangles.begin(); I != _triangles.end(); I++) {
      Ref<Triangle> tri = *I;
      for (j=0; j<3; j++) {
          Ref<Edge> e = tri->edge(j);
          int e_index = _edge_to_index[e];
          if (edge_to_triangle0[e_index] == 0) {
              edge_to_triangle0[e_index] = tri;
            }
          else if (edge_to_triangle1[e_index] == 0) {
              edge_to_triangle1[e_index] = tri;
            }
          else {
              ExEnv::errn() << "TriangulatedSurface::fix_orientation:"
                   << " more than two triangles to an edge" << endl;
              abort();
            }
        }
    }

  std::set<Ref<Triangle> > unfixed;
  std::set<Ref<Triangle> > fixed;
  std::set<Ref<Triangle> > finished;

  unfixed.insert(_triangles.begin(), _triangles.end());

  while (unfixed.size()) {
      // define unfixed.first()'s orientation to be the fixed orientation
      std::set<Ref<Triangle> >::iterator first = unfixed.begin();
      fixed.insert(*first);
      unfixed.erase(*first);
      while (fixed.size()) {
          Ref<Triangle> tri = *fixed.begin();
          // make all neighbors of tri oriented the same as tri
          for (i=0; i<3; i++) {
              Ref<Edge> e = tri->edge(i);
              int e_index = _edge_to_index[e];
              Ref<Triangle> othertri;
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
                  ExEnv::errn() << "TriangulatedSurface::fix_orientation: "
                       << " edge_to_triangle wrong" << endl;
                  abort();
                }
              if (tri->orientation(i) == othertri->orientation(j)) {
                  if (unfixed.find(othertri) != unfixed.end()) {
                      unfixed.erase(othertri);
                      fixed.insert(othertri);
                      othertri->flip();
                      nflip++;
                    }
                  else {
                      ExEnv::errn() << "TriangulatedSurface::fix_orientation:"
                           << " tried to flip a fixed triangle" << endl;
                      abort();
                    }
                }
              else if (unfixed.find(othertri) != unfixed.end()) {
                  unfixed.erase(othertri);
                  fixed.insert(othertri);
                }
            }
          // by this point all of tri's neighbors have had their
          // orientation fixed to match that of tri
          fixed.erase(tri);
          finished.insert(tri);
        }
    }

  if (_verbose) {
      ExEnv::outn()
          << scprintf("%d out of %d triangles were flipped\n", nflip, ntri);
    }

  // just in case
  for (i=0; i<ne; i++) {
      edge_to_triangle0[i] = 0;
      edge_to_triangle1[i] = 0;
    }

  delete[] edge_to_triangle0;
  delete[] edge_to_triangle1;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
