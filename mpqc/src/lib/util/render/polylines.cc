//
// polylines.cc
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

#include <stdlib.h>

#include <util/render/render.h>
#include <util/render/object.h>
#include <util/render/polylines.h>
#include <util/render/color.h>

using namespace sc;

static ClassDesc RenderedPolylines_cd(
  typeid(RenderedPolylines),"RenderedPolylines",1,"public RenderedObject",
  0, create<RenderedPolylines>, 0);

RenderedPolylines::RenderedPolylines()
{
  nvertex_ = 0;
  vertices_ = 0;
  vertex_rgb_ = 0;
  polylines_ = 0;
  nvertex_in_polyline_ = 0;
}

RenderedPolylines::RenderedPolylines(const Ref<KeyVal>& keyval):
  RenderedObject(keyval)
{
  int nvertex = keyval->count("vertices");
  int nline = keyval->count("lines");
  Coloring coloring = None;
  if (keyval->count("vertex_color_list")) {
      coloring = Vertex;
    }
  initialize(nvertex, nline, coloring);

  int i;
  for (i=0; i<nvertex; i++) {
      set_vertex(i,
                 keyval->doublevalue("vertices", i, 0),
                 keyval->doublevalue("vertices", i, 1),
                 keyval->doublevalue("vertices", i, 2));
    }

  if (coloring == Vertex) {
      for (i=0; i<nvertex; i++) {
          set_vertex_rgb(i,
                         keyval->doublevalue("vertex_color_list", i, 0),
                         keyval->doublevalue("vertex_color_list", i, 1),
                         keyval->doublevalue("vertex_color_list", i, 2));
        }
    }

  for (i=0; i<nline; i++) {
      nvertex_in_polyline_[i] = keyval->count("lines", i);
      polylines_[i] = new int[nvertex_in_polyline_[i]];
      for (int j=0; j<nvertex_in_polyline_[i]; j++) {
          polylines_[i][j] = keyval->intvalue("lines", i, j);
        }
    }
}

RenderedPolylines::~RenderedPolylines()
{
  if (vertices_ && vertices_[0]) delete[] vertices_[0];
  if (vertices_) delete[] vertices_;
  if (vertex_rgb_ && vertex_rgb_[0]) delete[] vertex_rgb_[0];
  if (vertex_rgb_) delete[] vertex_rgb_;
  if (polylines_) {
      for (int i=0; i<npolyline_; i++) {
          if (polylines_[i]) delete[] polylines_[i];
        }
      delete[] polylines_;
    }
  if (nvertex_in_polyline_) {
      delete[] nvertex_in_polyline_;
    }
}

void
RenderedPolylines::render(const Ref<Render>& render)
{
  render->polylines(this);
}

void
RenderedPolylines::initialize(int nvertex, int nline,
                             RenderedPolylines::Coloring coloring)
{
  coloring_ = coloring;
  nvertex_ = nvertex;
  npolyline_ = nline;
  
  vertices_ = new double*[nvertex];
  double* tmp = vertices_[0] = new double[3*nvertex];
  int i;
  for (i=1; i<nvertex; i++) {
      tmp += 3;
      vertices_[i] = tmp;
    }

  if (coloring == Vertex) {
      vertex_rgb_ = new double*[nvertex];
      double*tmp = vertex_rgb_[0] = new double[3*nvertex];
      for (i=1; i<nvertex; i++) {
          tmp += 3;
          vertex_rgb_[i] = tmp;
        }
    }
  else {
      vertex_rgb_ = 0;
    }

  polylines_ = new int*[nline];
  nvertex_in_polyline_ = new int[nline];
  for (i=0; i<nline; i++) {
      polylines_[i] = 0;
      nvertex_in_polyline_[i] = 0;
    }
}

void
RenderedPolylines::set_vertex(int i, double x, double y, double z)
{
  vertices_[i][0] = x;
  vertices_[i][1] = y;
  vertices_[i][2] = z;
}

void
RenderedPolylines::set_vertex_rgb(int i, double r, double g, double b)
{
  vertex_rgb_[i][0] = r;
  vertex_rgb_[i][1] = g;
  vertex_rgb_[i][2] = b;
}

void
RenderedPolylines::set_polyline(int i, int v1, int v2, int v3)
{
  if (polylines_[i]) delete[] polylines_[i];
  polylines_[i] = new int[3];
  polylines_[i][0] = v1;
  polylines_[i][1] = v2;
  polylines_[i][2] = v3;

  nvertex_in_polyline_[i] = 3;
}

void
RenderedPolylines::set_polyline(int i, int v1, int v2)
{
  if (polylines_[i]) delete[] polylines_[i];
  polylines_[i] = new int[2];
  polylines_[i][0] = v1;
  polylines_[i][1] = v2;

  nvertex_in_polyline_[i] = 2;
}

void
RenderedPolylines::set_polyline(int i, int v1, int v2, int v3, int v4)
{
  if (polylines_[i]) delete[] polylines_[i];
  polylines_[i] = new int[4];
  polylines_[i][0] = v1;
  polylines_[i][1] = v2;
  polylines_[i][2] = v3;
  polylines_[i][3] = v4;

  nvertex_in_polyline_[i] = 4;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
