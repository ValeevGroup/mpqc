//
// oogl.h
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@ca.sandia.gov>
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

#ifndef _util_render_oogl_h
#define _util_render_oogl_h

#include <stdio.h>
#include <util/render/render.h>

class OOGLRender: public Render {
#   define CLASSNAME OOGLRender
#   define HAVE_KEYVAL_CTOR
#   include <util/class/classd.h>
  private:
    char* filename_;
    FILE* fp_;
    int oogl_spheres_;
  public:
    OOGLRender(const char * filename);
    OOGLRender(FILE * fp = stdout);
    OOGLRender(const RefKeyVal&);
    virtual ~OOGLRender();

    void render(const RefRenderedObject&);
    void clear();

    void set(const RefRenderedObjectSet&);
    void sphere(const RefRenderedSphere&);
    void polygons(const RefRenderedPolygons&);
    void polylines(const RefRenderedPolylines&);
};
DescribedClass_REF_dec(OOGLRender);

#endif

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")
// End:
