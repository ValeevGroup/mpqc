//
// parameter.h
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

#ifndef _util_render_parameter_h
#define _util_render_parameter_h

#include <util/keyval/keyval.h>

namespace sc {

template <class T>
class Parameter {
    T parameter_;
    int is_set_;
    int overrides_;
  public:
    Parameter(): is_set_(0), overrides_(0) {}
    void set(const T& a) { parameter_ = a; is_set_ = 1; }
    void override(int overrides = 1) { overrides_ = overrides; }
    const T& value() const { return parameter_; }
    int overrides() const { return overrides_; }
    int is_set() const { return is_set_; }
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
