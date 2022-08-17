//
// tformv3.h
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

#ifndef _chemistry_qc_intv3_tranform_h
#define _chemistry_qc_intv3_tranform_h

#include <chemistry/qc/basis/gaussshell.h>
#include <chemistry/qc/basis/transform.h>
#include <chemistry/qc/intv3/macros.h>
#include <chemistry/qc/intv3/int1e.h>
#include <chemistry/qc/intv3/int2e.h>

namespace sc {

class Integral;

class SphericalTransformComponentV3 : public SphericalTransformComponent {
  public:
    void init(int a, int b, int c, double coef, int pureindex) {
      a_ = a;
      b_ = b;
      c_ = c;
      coef_ = coef;
      pureindex_ = pureindex;
      cartindex_ = INT_CARTINDEX(a+b+c,a,b);
    }
};

class SphericalTransformV3 : public SphericalTransform {
  public:
    SphericalTransformV3(int l, int subl=-1):SphericalTransform(l,subl) {
      init();
    }

    SphericalTransformComponent * new_component() {
      return new SphericalTransformComponentV3;
    }
};

class ISphericalTransformV3 : public ISphericalTransform {
  public:
    ISphericalTransformV3(int l, int subl=-1):ISphericalTransform(l,subl) {
      init();
    }

    SphericalTransformComponent * new_component() {
      return new SphericalTransformComponentV3;
    }
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
