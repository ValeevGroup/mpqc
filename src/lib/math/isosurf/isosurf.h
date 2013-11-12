//
// isosurf.h
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

#ifndef _math_isosurf_isosurf_h
#define _math_isosurf_isosurf_h

#include <vector>

#include <math/isosurf/surf.h>

namespace sc {

class IsosurfaceGen {
  protected:
    double _resolution;
  public:
    IsosurfaceGen();
    virtual ~IsosurfaceGen();
    virtual void isosurface(double value,
                            TriangulatedSurface& surf) = 0;
    virtual void set_resolution(double);
};

class ImplicitSurfacePolygonizer: public IsosurfaceGen {
    // These static data and members are used to interface to the
    // implicit.c routine provided in Graphics Gems IV.
    static ImplicitSurfacePolygonizer* current;
    // The following should not really be used publically.
    // they are public to permit access through internal
    // C-language functions.
  public:
    /// For internal use only.
    static int add_triangle_to_current(int,int,int,sc::detail::VERTICES);
    /// For internal use only.
    static double value_of_current(double x, double y, double z);
  protected:
    Ref<Volume> _volume;

    std::vector<Ref<Vertex> >  _tmp_vertices;

    TriangulatedSurface* _surf;
    double _value;
  public:
    ImplicitSurfacePolygonizer(const Ref<Volume>&);
    virtual ~ImplicitSurfacePolygonizer();
    virtual void isosurface(double value,
                            TriangulatedSurface& surf);
};  

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
