//
// polygons.h
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

#ifndef _util_render_polygons_h
#define _util_render_polygons_h

#include <util/keyval/keyval.h>
#include <util/render/object.h>

namespace sc {

class RenderedPolygons: public RenderedObject {
  public:
    enum Coloring { None, Vertex /*, Face*/ };
  private:
    RenderedPolygons::Coloring coloring_;
    int nvertex_;
    int nface_;
    double** vertices_;
    double** vertex_rgb_;
    int** faces_;
    int* nvertex_in_face_;
  protected:
    void render(const Ref<Render>&);
  public:
    RenderedPolygons();
    RenderedPolygons(const Ref<KeyVal>&);
    ~RenderedPolygons();

    void initialize(int nvertex, int nface,
                    RenderedPolygons::Coloring c = RenderedPolygons::None);
    void set_vertex(int, double x, double y, double z);
    void set_vertex_rgb(int, double r, double g, double b);
    void set_vertex_color(int i, const Color&c) {
        set_vertex_rgb(i, c.red(), c.green(), c.blue());
      }
    void set_face(int iface, int v1, int v2, int v3);
    void set_face(int iface, int v1, int v2, int v3, int v4);

    int nvertex() const { return nvertex_; }
    int nface() const { return nface_; }
    int nvertex_in_face(int iface) const { return nvertex_in_face_[iface]; }
    const double *vertex(int i) const { return vertices_[i]; }
    double vertex(int i, int j) const { return vertices_[i][j]; }
    int face(int i,int j) const { return faces_[i][j]; }
    double vertex_rgb(int i, int j) const { return vertex_rgb_[i][j]; }
    int have_vertex_rgb() const { return coloring_ == Vertex; }
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
