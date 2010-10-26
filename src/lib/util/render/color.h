//
// color.h
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

#ifndef _util_render_color_h
#define _util_render_color_h

#include <util/keyval/keyval.h>

namespace sc {

class Color {
  private:
    double red_;
    double green_;
    double blue_;
  public:
    Color() {}
    Color(double r, double g, double b): red_(r), green_(g), blue_(b) {}
    Color(const Ref<KeyVal>&);
    double red() const { return red_; }
    double green() const { return green_; }
    double blue() const { return blue_; }
    void set_rgb(double r, double g, double b) {
        red_ = r;
        green_ = g;
        blue_ = b;
      }
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
