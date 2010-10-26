//
// transform.h
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

#ifndef _util_render_transform_h
#define _util_render_transform_h

#include <iostream>

#include <util/class/class.h>
#include <util/keyval/keyval.h>
#include <util/render/algebra3.h>

namespace sc {

class Transform: public DescribedClass {
  private:
    mat4 transform_;
  public:
    Transform() { transform_ = identity3D(); }
    Transform(const Ref<KeyVal>&);
    ~Transform();
    mat4& transform() { return transform_; }
    void translate(double, double, double);
    void translate(const vec3&);
    void rotate(const vec3&, double angle_degrees);
    void scale(double);
    void print(std::ostream& = ExEnv::out0()) const;
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
