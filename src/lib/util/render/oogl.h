//
// oogl.h
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

#ifndef _util_render_oogl_h
#define _util_render_oogl_h

#include <iostream>
#include <util/render/render.h>

namespace sc {

class OOGLRender: public FileRender {
  private:
    int oogl_spheres_;
  public:
    OOGLRender(const char * filename);
    OOGLRender(std::ostream &o = ExEnv::out0());
    OOGLRender(const Ref<KeyVal>&);
    virtual ~OOGLRender();

    void render(const Ref<RenderedObject>&);
    void animate(const Ref<AnimatedObject>&);

    const char *file_extension();

    void set(const Ref<RenderedObjectSet>&);
    void sphere(const Ref<RenderedSphere>&);
    void polygons(const Ref<RenderedPolygons>&);
    void polylines(const Ref<RenderedPolylines>&);
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
