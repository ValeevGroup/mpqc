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

#ifdef __GNUC__
#pragma interface
#endif

#include <math/isosurf/surf.h>
#include  <math/isosurf/edgeRAVLMap.h>

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
  private:
    // These static data and members are used to interface to the
    // implicit.c routine provided in Graphics Gems IV.
    static ImplicitSurfacePolygonizer* current;
    static int add_triangle_to_current(int,int,int,VERTICES);
    static double value_of_current(double x, double y, double z);
  protected:
    RefVolume _volume;

    ArraysetRefVertex  _tmp_vertices;
    TriangulatedSurface* _surf;
    double _value;
  public:
    ImplicitSurfacePolygonizer(const RefVolume&);
    virtual ~ImplicitSurfacePolygonizer();
    virtual void isosurface(double value,
                            TriangulatedSurface& surf);
};  

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
