//
// vertex.h
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

#ifndef _math_isosurf_vertex_h
#define _math_isosurf_vertex_h

#ifdef __GNUC__
#pragma interface
#endif

#include <stdlib.h>
#include <math.h>
#include <util/ref/ref.h>
#include <math/scmat/matrix.h>
#include <math/isosurf/volume.h>
#include <math/isosurf/implicit.h>

namespace sc {

class Vertex: public RefCount {
  private:
    SCVector3 _point;
    SCVector3 *_normal; // _normal is optional
  public:
    Vertex();
    Vertex(const SCVector3& point,const SCVector3& normal);
    Vertex(const SCVector3& point);
    ~Vertex();
    const SCVector3& point() const { return _point; }
    int has_normal() const { return _normal != 0; }
    const SCVector3& normal() const { return *_normal; }
    void set_point(const SCVector3&p);
    void set_normal(const SCVector3&p);
    operator SCVector3&();

    void print(std::ostream&o=ExEnv::out0());
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
