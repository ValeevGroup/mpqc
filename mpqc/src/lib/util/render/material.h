//
// material.h
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

#ifndef _util_render_material_h
#define _util_render_material_h

#include <iostream>

#include <util/class/class.h>
#include <util/keyval/keyval.h>
#include <util/render/parameter.h>
#include <util/render/color.h>

namespace sc {

class Material: public DescribedClass {
  private:
    Parameter<Color> diffuse_;
    Parameter<Color> ambient_;
  public:
    Material();
    Material(const Ref<KeyVal>&);
    ~Material();
    Parameter<Color>& diffuse() { return diffuse_; }
    Parameter<Color>& ambient() { return ambient_; }
    void print(std::ostream& = ExEnv::out0()) const;
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
