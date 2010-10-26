//
// polylines.h
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

#ifndef _util_render_polylines_h
#define _util_render_polylines_h

#include <util/keyval/keyval.h>
#include <util/render/object.h>

namespace sc {

class RenderedPolylines: public RenderedObject {
  public:
    enum Coloring { None, Vertex };
  private:
    RenderedPolylines::Coloring coloring_;
    int nvertex_;
    int npolyline_;
    int *nvertex_in_polyline_;
    int **polylines_;
    double **vertices_;
    double **vertex_rgb_;
  public:
    RenderedPolylines();
    RenderedPolylines(const Ref<KeyVal>&);
    ~RenderedPolylines();

    void initialize(int nvertex, int npolylines,
                    RenderedPolylines::Coloring c = RenderedPolylines::None);
    int npolyline() { return npolyline_; }
    int nvertex() { return nvertex_; }
    int nvertex_in_polyline(int i) const { return nvertex_in_polyline_[i]; }
    double vertex(int i, int j) const { return vertices_[i][j]; }
    double vertex_rgb(int i, int j) const { return vertex_rgb_[i][j]; }
    int polyline(int i, int j) const { return polylines_[i][j]; }
    int have_vertex_rgb() const { return coloring_ == Vertex; }

    void set_vertex(int, double x, double y, double z);
    void set_vertex_rgb(int, double r, double g, double b);
    void set_polyline(int i, int v1, int v2);
    void set_polyline(int i, int v1, int v2, int v3);
    void set_polyline(int i, int v1, int v2, int v3, int v4);

    void render(const Ref<Render>&);
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
