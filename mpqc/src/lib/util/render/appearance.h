//
// appearance.h
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

#ifndef _util_render_appearance_h
#define _util_render_appearance_h

#include <iostream.h>

#include <util/class/class.h>
#include <util/keyval/keyval.h>
#include <util/render/parameter.h>

class Appearance: public DescribedClass {
#   define CLASSNAME Appearance
#   define HAVE_KEYVAL_CTOR
#   include <util/class/classd.h>
  private:
    Parameter<int> level_; // level of accuracy used to generate spheres, etc
  public:
    Appearance();
    Appearance(const RefKeyVal&);
    ~Appearance();
    Parameter<int>& level() { return level_; }

    void print(ostream& = cout);
};
DescribedClass_REF_dec(Appearance);

#endif

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")
// End:
