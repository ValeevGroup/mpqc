//
// edge.h
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

#ifndef _math_isosurf_edge_h
#define _math_isosurf_edge_h

#include <set>

#include <math/isosurf/vertex.h>

namespace sc {

class Edge: public RefCount {
  private:
    int _order;
    Ref<Vertex> *_vertices; // nvertices = _order + 1
  public:
    Edge(const Ref<Vertex> &p1,
         const Ref<Vertex> &p2)
    {
      _order = 1;
      _vertices = new Ref<Vertex>[2];
      _vertices[0]=p1; _vertices[1]=p2;
    }
    Edge(const Ref<Vertex> &p1,
         const Ref<Vertex> &p2,
         const Ref<Vertex> &p3);
    Edge(const Ref<Vertex> &p1,
         const Ref<Vertex> &p2,
         const Ref<Vertex> &p3,
         const Ref<Vertex> &p4);
    ~Edge();
    int order() const { return _order; }
    double straight_length();
    // return the endpoints
    Ref<Vertex> vertex(int i) const
    {
      return i?_vertices[_order]:_vertices[0];
    }
    // returns endpoints or interior vertex 0 <= i <= order
    Ref<Vertex> interior_vertex(int i) const
    {
      return _vertices[i];
    }
    // add the endpoints to the set
    void add_vertices(std::set<Ref<Vertex> >&);
    void set_order(int order, const Ref<Volume>&vol,double isovalue);
    // find the position of a point on the edge
    int interpolate(double location, SCVector3&point, SCVector3&norm);
    // find the true position of a point using the isosurface
    int interpolate(double location, SCVector3&point, SCVector3&norm,
                     const Ref<Volume> &vol, double isovalue);
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
